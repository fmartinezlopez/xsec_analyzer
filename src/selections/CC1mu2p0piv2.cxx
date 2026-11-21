// Standard library includes
#include <cmath>

// XSecAnalyzer includes
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/KICalculator.hh"
#include "XSecAnalyzer/Selections/CC1mu2p0piv2.hh"

namespace {
  constexpr double MUON_P_MIN = 0.10; // GeV/c

  // Proton momentum window, used both in the signal definition and in the
  // reco selection. NOTE: The proton multiplicity analysis (CC1muNp0pi) uses
  // 0.25 GeV/c. The default here was raised to 0.30 GeV/c, which roughly
  // halves the CC0pi >= 3p background (missed third proton just above
  // threshold). Set this back to 0.25 to recover the exact correspondence
  // with the multiplicity 2p bin.
  constexpr double PROTON_P_MIN = 0.30; // GeV/c
  constexpr double PROTON_P_MAX = 1.00; // GeV/c

  constexpr double PI_PLUS_P_MIN = 0.00; // GeV/c

  // Proton candidate attachment to the reco neutrino vertex: the start of the
  // closer proton candidate must lie within PROTON_VTX_DIST_NEAR of the
  // vertex, and the start of the farther one within PROTON_VTX_DIST_FAR.
  // (The 2022 requirement of 4 cm for both loses collinear pairs.)
  constexpr double PROTON_VTX_DIST_NEAR = 4.; // cm
  constexpr double PROTON_VTX_DIST_FAR = 8.; // cm

  // NuGraph2 particle ID, which replaces the LLR proton PID cut
  // (DEFAULT_PROTON_PID_CUT is still evaluated and stored, but not applied).
  // Both proton candidates must have at least this fraction of their hits
  // labelled HIP by NuGraph2. Unlike the LLR score, this works for short
  // (< 5 cm) tracks, which is where most subleading protons are.
  constexpr double PROTON_NG2_HIP_FRAC_MIN = 0.6;
  // The muon candidate must carry the NuGraph2 MIP semantic label
  // (pfng2semlabel: 0 MIP, 1 HIP, 2 shower, 3 Michel, 4 diffuse)
  constexpr int NG2_MIP_LABEL = 0;

  // Proton candidates shorter than this use the vertex-to-track-end
  // direction instead of the track-fit direction (see get_reco_p3_proton)
  constexpr double PROTON_SHORT_TRACK_LENGTH = 12.; // cm

}

CC1mu2p0piv2::CC1mu2p0piv2( const std::string& sel_name )
  : SelectionBase( sel_name )
{
  ki_calc_type_ = KICalculator::kOpt1;
  this->define_fv( 21.5, 234.85, -95.0, 95.0, 21.5, 966.8 );
}

bool CC1mu2p0piv2::is_selected( AnalysisEvent& ev ) {

  // Get access to the input PeLEE ntuple branches. The tree is requested by
  // name because derived selections add other input trees (ev.in() without
  // a name returns the first one in alphabetical order).
  const auto& in = ev.in( "nuselection/NeutrinoSelectionFilter" );

  // Get access to the output TTree branches
  auto& out = ev.out();

  FiducialVolume pcv( 10., 246.35, -106.5, 106.5, 10., 1026.8 );

  // Require event to pass original cosmic rejection (removed for NuGraph2 pass)
  int slice_orig_pass_id;
  in.at( "slice_orig_pass_id" ) >> slice_orig_pass_id;
  bool sel_slice_orig_pass = ( slice_orig_pass_id == 1 );

  // Check whether the space-charge-corrected position of the reconstructed
  // neutrino vertex lies within the fiducial volume
  float nu_vx, nu_vy, nu_vz;
  in.at( "reco_nu_vtx_sce_x" ) >> nu_vx;
  in.at( "reco_nu_vtx_sce_y" ) >> nu_vy;
  in.at( "reco_nu_vtx_sce_z" ) >> nu_vz;

  bool sel_reco_vertex_in_FV = this->get_fv().is_inside( nu_vx, nu_vy, nu_vz );

  // Evaluate cuts on the topological score cut and the cosmic impact parameter
  float topological_score, cosmic_ip;
  in.at( "topological_score" ) >> topological_score;
  in.at( "CosmicIP" ) >> cosmic_ip;
  bool sel_topo_cut_passed = topological_score > TOPO_SCORE_CUT;
  bool sel_cosmic_ip_cut_passed = cosmic_ip > COSMIC_IP_CUT;

  // Retrieve information about the PFParticles in the event
  int num_pfps;
  std::vector< unsigned int >* gen_vec;
  std::vector< float > *track_start_x, *track_start_y, *track_start_z;
  in.at( "n_pfps" ) >> num_pfps;
  in.at( "pfp_generation_v" ) >> gen_vec;
  in.at( "trk_sce_start_x_v" ) >> track_start_x;
  in.at( "trk_sce_start_y_v" ) >> track_start_y;
  in.at( "trk_sce_start_z_v" ) >> track_start_z;

  // Evaluate the containment cut by looping over all PFParticles and checking
  // their starting positions. Pass the cut by default. Also count the direct
  // neutrino daughters while we're at it.
  bool sel_pfp_starts_in_PCV = true;
  int num_gen2_pfps = 0;
  for ( int p = 0; p < num_pfps; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = gen_vec->at( p );
    if ( generation != 2u ) continue;

    ++num_gen2_pfps;

    // Use the track reconstruction results to get the start point for
    // every PFParticle for the purpose of verifying containment. We could
    // in principle differentiate between tracks and showers here, but
    // (1) we cut out all showers later on in the selection anyway, and
    // (2) the blinded PeLEE data ntuples do not include shower information.
    // We therefore apply the track reconstruction here unconditionally.
    float x = track_start_x->at( p );
    float y = track_start_y->at( p );
    float z = track_start_z->at( p );

    // Verify that the start of the PFParticle lies within the containment
    // volume.
    //
    // See https://stackoverflow.com/a/2488507 for an explanation of the use
    // of &= here. Don't worry, it's type-safe since both operands are bool.
    sel_pfp_starts_in_PCV &= pcv.is_inside( x, y, z );
  }

  // Search for a muon candidate among the PFParticles in the event. A cut on
  // the muon momentum is handled later.
  std::vector< float > *track_score, *track_start_distance, *track_length,
    *track_llr_pid_score;
  in.at( "trk_score_v" ) >> track_score;
  in.at( "trk_distance_v" ) >> track_start_distance;
  in.at( "trk_len_v" ) >> track_length;
  in.at( "trk_llr_pid_score_v" ) >> track_llr_pid_score;

  std::vector< int > muon_candidate_indices;
  // NOTE: CC1muNp0pi stores these scores in a std::vector< int >, which
  // truncates them. A std::vector< float > is used here so that the
  // comparison below works as intended.
  std::vector< float > muon_pid_scores;

  bool sel_has_muon_candidate = false;

  for ( int p = 0; p < num_pfps; ++p ) {
    // Only direct neutrino daughters (generation == 2) will be considered as
    // possible muon candidates
    unsigned int generation = gen_vec->at( p );
    if ( generation != 2u ) continue;

    float trk_score = track_score->at( p );
    float start_dist = track_start_distance->at( p );
    float trk_length = track_length->at( p );
    float pid_score = track_llr_pid_score->at( p );

    if ( trk_score > MUON_TRACK_SCORE_CUT
      && start_dist < MUON_VTX_DISTANCE_CUT
      && trk_length > MUON_LENGTH_CUT
      && pid_score > MUON_PID_CUT )
    {
      muon_candidate_indices.push_back( p );
      muon_pid_scores.push_back( pid_score );
    }
  }

  size_t num_candidates = muon_candidate_indices.size();
  if ( num_candidates > 0u ) sel_has_muon_candidate = true;

  // By default, set the muon candidate index to the bogus value
  int muon_candidate_idx = BOGUS_INDEX;
  // If exactly one muon candidate was found, then assign the index of its
  // corresponding PFParticle
  if ( num_candidates == 1u ) {
    muon_candidate_idx = muon_candidate_indices.front();
  }
  // In the case of multiple muon candidates, choose the one with the highest
  // PID score (most muon-like) as the one to use
  else if ( num_candidates > 1u ) {
    float highest_score = LOW_FLOAT;
    int chosen_index = BOGUS_INDEX;
    for ( size_t c = 0; c < num_candidates; ++c ) {
      float score = muon_pid_scores.at( c );
      if ( highest_score < score ) {
        highest_score = score;
        chosen_index = muon_candidate_indices.at( c );
      }
    }
    muon_candidate_idx = chosen_index;
  }

  // Combine some of the prior cuts into a numu CC preselection
  bool sel_nu_mu_cc = sel_slice_orig_pass && sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV
    && sel_has_muon_candidate && sel_topo_cut_passed;

  // Fail the shower cut if any showers were reconstructed
  // NOTE: We could do this quicker like this,
  //   sel_no_reco_showers_ = ( num_showers_ > 0 );
  // but it might be nice to be able to adjust the track score for this cut.
  // Therefore, we do it the hard way.
  int reco_shower_count = 0;
  for ( int p = 0; p < num_pfps; ++p ) {

    // Only check direct neutrino daughters (generation == 2)
    unsigned int generation = gen_vec->at( p );
    if ( generation != 2u ) continue;

    float tscore = track_score->at( p );
    if ( tscore <= TRACK_SCORE_CUT ) ++reco_shower_count;
  }
  // Check the shower cut
  bool sel_no_reco_showers = ( reco_shower_count == 0 );

  // Set default values of the next few cuts before checking them
  bool sel_passed_proton_pid_cut = true;
  bool sel_protons_contained = true;
  bool sel_protons_passed_mom_cuts = true;
  bool sel_muon_contained = false;
  bool sel_muon_passed_mom_cuts = false;
  bool sel_muon_quality_ok = false;
  bool sel_has_p_candidate = false;

  std::vector< float > *track_end_x, *track_end_y, *track_end_z,
    *muon_mom_range, *muon_mom_mcs, *track_KE_proton;
  in.at( "trk_sce_end_x_v" ) >> track_end_x;
  in.at( "trk_sce_end_y_v" ) >> track_end_y;
  in.at( "trk_sce_end_z_v" ) >> track_end_z;
  in.at( "trk_range_muon_mom_v" ) >> muon_mom_range;
  in.at( "trk_mcs_muon_mom_v" ) >> muon_mom_mcs;
  in.at( "trk_energy_proton_v" ) >> track_KE_proton;

  // Both muon momentum estimators are stored in the output tree, together
  // with their fractional difference
  double p_mu_range = BOGUS, p_mu_mcs = BOGUS, frac_diff_range_mcs = BOGUS;

  // Indices and range-based momenta of the proton candidates
  std::vector< int > p_candidate_indices;
  std::vector< float > p_candidate_moms;

  for ( int p = 0; p < num_pfps; ++p ) {
    // Only worry about direct neutrino daughters (PFParticles considered
    // daughters of the reconstructed neutrino)
    unsigned int generation = gen_vec->at( p );
    if ( generation != 2u ) continue;

    // Check that we can find a muon candidate in the event. If more than
    // one is found, also fail the cut.
    if ( p == muon_candidate_idx ) {

      // Check whether the muon candidate is contained. Use the same
      // containment volume as the protons. TODO: revisit this as needed.
      float endx = track_end_x->at( p );
      float endy = track_end_y->at( p );
      float endz = track_end_z->at( p );
      bool end_contained = pcv.is_inside( endx, endy, endz );

      if ( end_contained ) sel_muon_contained = true;

      // Check that the muon candidate is above threshold. Use the best
      // momentum based on whether it was contained or not. No upper limit
      // is applied.
      float muon_mom = LOW_FLOAT;
      float range_muon_mom = muon_mom_range->at( p );
      float mcs_muon_mom = muon_mom_mcs->at( p );

      p_mu_range = range_muon_mom;
      p_mu_mcs = mcs_muon_mom;

      if ( sel_muon_contained ) muon_mom = range_muon_mom;
      else muon_mom = mcs_muon_mom;

      if ( muon_mom >= MUON_P_MIN ) {
        sel_muon_passed_mom_cuts = true;
      }

      // Evaluate the muon candidate quality cut by comparing MCS and
      // range-based momentum estimators. The flag is stored but is not
      // (for now) part of the selection.
      if ( range_muon_mom > 0. ) {
        frac_diff_range_mcs = std::abs( range_muon_mom - mcs_muon_mom )
          / range_muon_mom;
        if ( frac_diff_range_mcs < MUON_MOM_QUALITY_CUT ) {
          sel_muon_quality_ok = true;
        }
      }
    }
    else {

      float trk_score = track_score->at( p );
      if ( trk_score <= TRACK_SCORE_CUT ) continue;

      // Bad tracks in the searchingfornues TTree can have
      // bogus track lengths. This skips those.
      float trk_length = track_length->at( p );
      if ( trk_length <= 0. ) continue;

      // We found a reco track that is not the muon candidate. All such
      // tracks are considered proton candidates.
      sel_has_p_candidate = true;

      float llr_pid_score = track_llr_pid_score->at( p );

      // Check whether the current proton candidate fails the proton PID cut
      if ( llr_pid_score > DEFAULT_PROTON_PID_CUT ) {
        sel_passed_proton_pid_cut = false;
      }

      // Check whether the current proton candidate fails the containment cut
      float endx = track_end_x->at( p );
      float endy = track_end_y->at( p );
      float endz = track_end_z->at( p );
      bool end_contained = pcv.is_inside( endx, endy, endz );
      if ( !end_contained ) sel_protons_contained = false;

      // Check whether the range-based reco momentum of the current proton
      // candidate lies outside of the allowed range
      float KEp = track_KE_proton->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );
      if ( p_mom < PROTON_P_MIN || p_mom > PROTON_P_MAX ) {
        sel_protons_passed_mom_cuts = false;
      }

      p_candidate_indices.push_back( p );
      p_candidate_moms.push_back( p_mom );
    }

  }

  // Find the leading and subleading proton candidates, i.e., the two with
  // the highest range-based reco momenta
  int lead_p_candidate_idx = BOGUS_INDEX;
  int sub_p_candidate_idx = BOGUS_INDEX;

  float lead_p_mom = LOW_FLOAT;
  float sub_p_mom = LOW_FLOAT;

  size_t num_p_candidates = p_candidate_indices.size();
  for ( size_t c = 0u; c < num_p_candidates; ++c ) {
    float p_mom = p_candidate_moms.at( c );
    if ( p_mom > lead_p_mom ) {
      sub_p_mom = lead_p_mom;
      sub_p_candidate_idx = lead_p_candidate_idx;

      lead_p_mom = p_mom;
      lead_p_candidate_idx = p_candidate_indices.at( c );
    }
    else if ( p_mom > sub_p_mom ) {
      sub_p_mom = p_mom;
      sub_p_candidate_idx = p_candidate_indices.at( c );
    }
  }

  // Require exactly two proton candidates and no other direct neutrino
  // daughters (e.g., ones skipped above due to a bogus track length)
  bool sel_two_p_candidates = ( num_p_candidates == 2u );
  bool sel_no_extra_pfps = ( num_gen2_pfps == 3 );

  // Additional requirements from the original CC1mu2p0pi selection. These
  // are stored but not (yet) applied.

  // No PFParticles other than the direct neutrino daughters
  bool sel_no_daughter_pfps = ( num_pfps == num_gen2_pfps );

  // Both proton candidates start or end near the reco vertex, and both
  // pass the tighter track score cut
  TVector3 x3_reco_vtx( nu_vx, nu_vy, nu_vz );
  auto track_end_dist = [&]( int idx ) -> float {
    TVector3 x3_end( track_end_x->at( idx ), track_end_y->at( idx ),
      track_end_z->at( idx ) );
    return ( x3_end - x3_reco_vtx ).Mag();
  };

  bool sel_protons_attached = true;
  bool sel_protons_track_like = true;
  for ( size_t c = 0u; c < num_p_candidates; ++c ) {
    int p = p_candidate_indices.at( c );

    float start_dist = track_start_distance->at( p );
    float end_dist = track_end_dist( p );

    if ( start_dist > MUON_VTX_DISTANCE_CUT
      && end_dist > MUON_VTX_DISTANCE_CUT ) sel_protons_attached = false;

    if ( track_score->at( p ) < MUON_TRACK_SCORE_CUT ) {
      sel_protons_track_like = false;
    }
  }

  out[ "no_daughter_pfps" ] = sel_no_daughter_pfps;
  out[ "protons_attached" ] = sel_protons_attached;
  out[ "protons_track_like" ] = sel_protons_track_like;
  out[ "num_pfps" ] = num_pfps;

  // Requirements of the current working point (applied below)
  std::vector< unsigned int >* pfp_n_desc;
  in.at( "pfp_n_descendents_v" ) >> pfp_n_desc;

  // (1) Proton candidates attached to the vertex (see PROTON_VTX_DIST_*)
  // (2) No proton candidate has descendant PFParticles (targets pions that
  //     decay or interact, and protons that reinteract)
  float p_start_dist_min = BOGUS, p_start_dist_max = BOGUS;
  bool sel_protons_no_descendents = true;
  for ( size_t c = 0u; c < num_p_candidates; ++c ) {
    int p = p_candidate_indices.at( c );
    float sd = track_start_distance->at( p );
    if ( c == 0u || sd < p_start_dist_min ) p_start_dist_min = sd;
    if ( c == 0u || sd > p_start_dist_max ) p_start_dist_max = sd;
    if ( pfp_n_desc->at( p ) > 0u ) sel_protons_no_descendents = false;
  }
  bool sel_protons_near_vertex = ( num_p_candidates > 0u )
    && ( p_start_dist_min < PROTON_VTX_DIST_NEAR )
    && ( p_start_dist_max < PROTON_VTX_DIST_FAR );

  // (3) Range/MCS momentum agreement for contained muon candidates. Exiting
  // muons (MCS only) are not affected. Targets NC events in which a
  // contained pion or proton is taken as the muon.
  bool sel_muon_quality_if_contained = !sel_muon_contained
    || sel_muon_quality_ok;

  out[ "protons_near_vertex" ] = sel_protons_near_vertex;
  out[ "protons_no_descendents" ] = sel_protons_no_descendents;
  out[ "muon_quality_if_contained" ] = sel_muon_quality_if_contained;

  // NuGraph2 particle ID (applied below). If the NuGraph2 per-PFP vectors
  // are missing or shorter than expected (e.g., a sample without the NG2
  // pass), the requirements fail.
  std::vector< int >* ng2_semlabel;
  std::vector< float > *ng2_hipfrac, *ng2_mipfrac, *ng2_mclfrac,
    *ng2_bkgfrac;
  in.at( "pfng2semlabel" ) >> ng2_semlabel;
  in.at( "pfng2hipfrac" ) >> ng2_hipfrac;
  in.at( "pfng2mipfrac" ) >> ng2_mipfrac;
  in.at( "pfng2mclfrac" ) >> ng2_mclfrac;
  in.at( "pfng2bkgfrac" ) >> ng2_bkgfrac;

  auto ng2_ok = [&]( int idx ) -> bool {
    if ( idx == BOGUS_INDEX ) return false;
    size_t uidx = static_cast< size_t >( idx );
    return ng2_semlabel && ng2_hipfrac && ng2_semlabel->size() > uidx
      && ng2_hipfrac->size() > uidx;
  };

  // (4) Both proton candidates are HIP-like according to NuGraph2
  bool sel_protons_ng2_hip = ( num_p_candidates > 0u );
  for ( size_t c = 0u; c < num_p_candidates; ++c ) {
    int p = p_candidate_indices.at( c );
    if ( !ng2_ok( p ) || !( ng2_hipfrac->at( p ) >= PROTON_NG2_HIP_FRAC_MIN ) )
    {
      sel_protons_ng2_hip = false;
    }
  }

  // (5) The muon candidate is MIP-like according to NuGraph2
  bool sel_muon_ng2_mip = ng2_ok( muon_candidate_idx )
    && ( ng2_semlabel->at( muon_candidate_idx ) == NG2_MIP_LABEL );

  out[ "protons_ng2_hip" ] = sel_protons_ng2_hip;
  out[ "muon_ng2_mip" ] = sel_muon_ng2_mip;

  // Event-level quantities useful for offline cut studies
  int nslice;
  in.at( "nslice" ) >> nslice;
  out[ "nslice" ] = nslice;
  out[ "topological_score" ] = static_cast< double >( topological_score );
  out[ "cosmic_ip" ] = static_cast< double >( cosmic_ip );

  // Per-candidate quantities for the muon, leading and subleading proton
  // candidates. These allow PID, attachment and truth-matching studies to be
  // done offline.
  std::vector< float > *track_pid_chipr, *track_bragg_p, *track_bragg_mu;
  std::vector< unsigned int > *pfp_trk_daughters, *pfp_shr_daughters;
  in.at( "trk_pid_chipr_v" ) >> track_pid_chipr;
  in.at( "trk_bragg_p_v" ) >> track_bragg_p;
  in.at( "trk_bragg_mu_v" ) >> track_bragg_mu;
  in.at( "pfp_trk_daughters_v" ) >> pfp_trk_daughters;
  in.at( "pfp_shr_daughters_v" ) >> pfp_shr_daughters;

  // Backtracked truth information is MC-only. In data these branches may be
  // present but empty, so the sizes are checked before use.
  std::vector< int >* bt_pdg;
  std::vector< float > *bt_purity, *bt_completeness, *bt_px, *bt_py, *bt_pz;
  in.at( "backtracked_pdg" ) >> bt_pdg;
  in.at( "backtracked_purity" ) >> bt_purity;
  in.at( "backtracked_completeness" ) >> bt_completeness;
  in.at( "backtracked_px" ) >> bt_px;
  in.at( "backtracked_py" ) >> bt_py;
  in.at( "backtracked_pz" ) >> bt_pz;

  auto store_candidate = [&]( const std::string& prefix, int idx ) {
    double llr = BOGUS, tscore = BOGUS, length = BOGUS, start_dist = BOGUS,
      end_dist = BOGUS, chipr = BOGUS, bragg_p = BOGUS, bragg_mu = BOGUS,
      bt_pur = BOGUS, bt_comp = BOGUS, bt_mom = BOGUS;
    int n_desc = BOGUS_INDEX, n_trk_d = BOGUS_INDEX, n_shr_d = BOGUS_INDEX,
      bt_code = 0;
    bool contained = false;

    if ( idx != BOGUS_INDEX ) {
      llr = track_llr_pid_score->at( idx );
      tscore = track_score->at( idx );
      length = track_length->at( idx );
      start_dist = track_start_distance->at( idx );
      end_dist = track_end_dist( idx );
      chipr = track_pid_chipr->at( idx );
      bragg_p = track_bragg_p->at( idx );
      bragg_mu = track_bragg_mu->at( idx );
      n_desc = pfp_n_desc->at( idx );
      n_trk_d = pfp_trk_daughters->at( idx );
      n_shr_d = pfp_shr_daughters->at( idx );
      contained = pcv.is_inside( track_end_x->at( idx ),
        track_end_y->at( idx ), track_end_z->at( idx ) );

      size_t uidx = static_cast< size_t >( idx );
      if ( bt_pdg && bt_pdg->size() > uidx ) {
        bt_code = bt_pdg->at( idx );
        bt_pur = bt_purity->at( idx );
        bt_comp = bt_completeness->at( idx );
        bt_mom = TVector3( bt_px->at( idx ), bt_py->at( idx ),
          bt_pz->at( idx ) ).Mag();
      }
    }

    out[ prefix + "_llr_pid" ] = llr;
    out[ prefix + "_trk_score" ] = tscore;
    out[ prefix + "_length" ] = length;
    out[ prefix + "_start_dist" ] = start_dist;
    out[ prefix + "_end_dist" ] = end_dist;
    out[ prefix + "_chi2_p" ] = chipr;
    out[ prefix + "_bragg_p" ] = bragg_p;
    out[ prefix + "_bragg_mu" ] = bragg_mu;
    out[ prefix + "_n_descendents" ] = n_desc;
    out[ prefix + "_n_trk_daughters" ] = n_trk_d;
    out[ prefix + "_n_shr_daughters" ] = n_shr_d;
    out[ prefix + "_contained" ] = contained;
    out[ prefix + "_bt_pdg" ] = bt_code;
    out[ prefix + "_bt_purity" ] = bt_pur;
    out[ prefix + "_bt_completeness" ] = bt_comp;
    out[ prefix + "_bt_p" ] = bt_mom;

    // NuGraph2 per-PFP quantities
    int ng2_label = BOGUS_INDEX;
    double hipf = BOGUS, mipf = BOGUS, mclf = BOGUS, bkgf = BOGUS;
    if ( ng2_ok( idx ) ) {
      size_t uidx = static_cast< size_t >( idx );
      ng2_label = ng2_semlabel->at( idx );
      hipf = ng2_hipfrac->at( idx );
      if ( ng2_mipfrac && ng2_mipfrac->size() > uidx ) {
        mipf = ng2_mipfrac->at( idx );
      }
      if ( ng2_mclfrac && ng2_mclfrac->size() > uidx ) {
        mclf = ng2_mclfrac->at( idx );
      }
      if ( ng2_bkgfrac && ng2_bkgfrac->size() > uidx ) {
        bkgf = ng2_bkgfrac->at( idx );
      }
    }
    out[ prefix + "_ng2_semlabel" ] = ng2_label;
    out[ prefix + "_ng2_hipfrac" ] = hipf;
    out[ prefix + "_ng2_mipfrac" ] = mipf;
    out[ prefix + "_ng2_mclfrac" ] = mclf;
    out[ prefix + "_ng2_bkgfrac" ] = bkgf;
  };

  store_candidate( "mu", muon_candidate_idx );
  store_candidate( "lead_p", lead_p_candidate_idx );
  store_candidate( "sub_p", sub_p_candidate_idx );

  // All right, we've applied all selection cuts. Set the flag that indicates
  // whether all were passed (and thus the event is selected as a CC2p0pi
  // candidate). Note that the muon is not required to be contained; the
  // range/MCS quality cut is applied only to contained muon candidates.
  // The proton PID is the NuGraph2 HIP requirement; the LLR proton PID flag
  // (sel_passed_proton_pid_cut) is stored but no longer applied.
  bool sel_CC2p0pi = sel_nu_mu_cc && sel_no_reco_showers
    && sel_muon_passed_mom_cuts && sel_muon_quality_if_contained
    && sel_muon_ng2_mip
    && sel_has_p_candidate && sel_two_p_candidates && sel_no_extra_pfps
    && sel_protons_ng2_hip && sel_protons_contained
    && sel_protons_passed_mom_cuts
    && sel_protons_near_vertex && sel_protons_no_descendents;

  // Store the muon and proton candidate indices in the output tree (they may
  // be BOGUS_INDEX if such a candidate was not found)
  out[ "muon_candidate_idx" ] = muon_candidate_idx;
  out[ "lead_p_candidate_idx" ] = lead_p_candidate_idx;
  out[ "sub_p_candidate_idx" ] = sub_p_candidate_idx;

  std::vector< float > *track_dir_x, *track_dir_y, *track_dir_z;
  in.at( "trk_dir_x_v" ) >> track_dir_x;
  in.at( "trk_dir_y_v" ) >> track_dir_y;
  in.at( "trk_dir_z_v" ) >> track_dir_z;

  // Helper function that builds a reco 3-momentum from the direction of the
  // track with the given index and the given momentum magnitude. As in
  // CC1mu2p0pi, the direction is reversed if the end of the track is closer
  // to the reco neutrino vertex than its start.
  TVector3 x3_nu_vertex( nu_vx, nu_vy, nu_vz );

  auto get_reco_p3 = [&]( int idx, float mom ) -> TVector3 {
    TVector3 p3( track_dir_x->at( idx ), track_dir_y->at( idx ),
      track_dir_z->at( idx ) );
    p3 = p3.Unit() * mom;

    float start_dist = track_start_distance->at( idx );

    TVector3 x3_end( track_end_x->at( idx ), track_end_y->at( idx ),
      track_end_z->at( idx ) );
    x3_end -= x3_nu_vertex;
    float end_dist = x3_end.Mag();

    if ( start_dist > end_dist ) p3 *= -1.0;

    return p3;
  };

  // Proton version: for short tracks (length < PROTON_SHORT_TRACK_LENGTH) the
  // track-fit direction is unreliable, so use the direction from the reco
  // neutrino vertex to the far end of the track (the end, or the start if
  // the track is flipped). For selected signal this halves the 90th
  // percentile of the subleading-proton angular error (43 -> 21 deg) and
  // improves the pair-variable resolutions by 10-20%. Longer tracks keep the
  // (flipped) track-fit direction, which is slightly better there.
  auto get_reco_p3_proton = [&]( int idx, float mom ) -> TVector3 {
    TVector3 p3 = get_reco_p3( idx, mom );
    if ( track_length->at( idx ) >= PROTON_SHORT_TRACK_LENGTH ) return p3;

    TVector3 x3_start( track_start_x->at( idx ), track_start_y->at( idx ),
      track_start_z->at( idx ) );
    TVector3 x3_end( track_end_x->at( idx ), track_end_y->at( idx ),
      track_end_z->at( idx ) );

    float start_dist = track_start_distance->at( idx );
    float end_dist = ( x3_end - x3_nu_vertex ).Mag();
    TVector3 x3_far = ( start_dist > end_dist ) ? x3_start : x3_end;

    TVector3 dir = x3_far - x3_nu_vertex;
    if ( dir.Mag() <= 0. ) return p3;

    return dir.Unit() * mom;
  };

  // Set the reco 3-momentum of the muon candidate if we found one
  TVector3 p3mu;
  bool muon = muon_candidate_idx != BOGUS_INDEX;
  if ( muon ) {

    // The selection flag indicating whether the muon candidate is contained
    // was already set when the selection was applied. Use it to choose the
    // best momentum estimator to use.
    float muon_mom = LOW_FLOAT;
    if ( sel_muon_contained ) {
      muon_mom = muon_mom_range->at( muon_candidate_idx );
    }
    else {
      muon_mom = muon_mom_mcs->at( muon_candidate_idx );
    }

    p3mu = get_reco_p3( muon_candidate_idx, muon_mom );
  }

  // Set the reco 3-momentum of the leading proton candidate if we found one
  TVector3 p3_lead_p;
  bool lead_p = lead_p_candidate_idx != BOGUS_INDEX;
  if ( lead_p ) {
    float KEp = track_KE_proton->at( lead_p_candidate_idx );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    p3_lead_p = get_reco_p3_proton( lead_p_candidate_idx, p_mom );
  }

  // Set the reco 3-momentum of the subleading proton candidate if we found one
  TVector3 p3_sub_p;
  bool sub_p = sub_p_candidate_idx != BOGUS_INDEX;
  if ( sub_p ) {
    float KEp = track_KE_proton->at( sub_p_candidate_idx );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

    p3_sub_p = get_reco_p3_proton( sub_p_candidate_idx, p_mom );
  }

  // Set the reco 3-momenta of all proton candidates (i.e., all generation == 2
  // tracks except the muon candidate) assuming we found both a muon candidate
  // and at least one proton candidate.
  std::vector< TVector3 > p3_p_vec;
  if ( muon && lead_p ) {
    for ( int p = 0; p < num_pfps; ++p ) {
      // Skip the muon candidate
      if ( p == muon_candidate_idx ) continue;

      // Only include direct neutrino daughters (generation == 2)
      unsigned int generation = gen_vec->at( p );
      if ( generation != 2u ) continue;

      float KEp = track_KE_proton->at( p );
      float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );

      p3_p_vec.push_back( get_reco_p3_proton( p, p_mom ) );
    }

    // Sort the reco proton 3-momenta in order from highest to lowest magnitude
    std::sort( p3_p_vec.begin(), p3_p_vec.end(), [](const TVector3& a,
      const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );
  }

  // Distance between the reco neutrino vertex and the start of each of the
  // two proton candidate tracks (useful for rejecting secondary protons)
  double lead_p_start_dist = BOGUS, sub_p_start_dist = BOGUS;
  if ( lead_p ) {
    lead_p_start_dist = track_start_distance->at( lead_p_candidate_idx );
  }
  if ( sub_p ) {
    sub_p_start_dist = track_start_distance->at( sub_p_candidate_idx );
  }


  // Sideband topologies, for background constraints. Every direct neutrino
  // daughter other than the muon candidate is classified with the cuts of
  // this selection: a shower (track score <= TRACK_SCORE_CUT), a bad track
  // (no valid length), or a track; a track is a good proton if contained,
  // NuGraph2 HIP fraction >= PROTON_NG2_HIP_FRAC_MIN, without descendants
  // and with range-based p in the proton window; it is pion-like if it
  // fails the HIP requirement, carries the NuGraph2 MIP label and starts
  // within PROTON_VTX_DIST_NEAR of the vertex. All sidebands need the
  // selection's muon, no shower or bad track, and good protons near the
  // vertex as in the selection. Their kinematics (sb_* branches) use the
  // muon and the two "protons" of the topology: the two good protons (3p:
  // the hardest two; one good proton and a pion: both, the pion taken as a
  // proton), or the single good proton.
  bool sb_muon_ok = sel_nu_mu_cc && sel_muon_passed_mom_cuts
    && sel_muon_quality_if_contained && sel_muon_ng2_mip;
  int sb_n_trk = 0, sb_n_bad = 0, sb_n_shw = 0;
  std::vector< std::pair< float, int > > sb_good, sb_pion; // (p, PFP index)
  float sb_sd_min = 1e9, sb_sd_max = -1e9;
  for ( int p = 0; p < num_pfps; ++p ) {
    if ( gen_vec->at( p ) != 2u || p == muon_candidate_idx ) continue;
    float ts = track_score->at( p );
    if ( ts <= TRACK_SCORE_CUT ) { ++sb_n_shw; continue; }
    if ( !( ts > TRACK_SCORE_CUT ) ) continue;
    if ( !( track_length->at( p ) > 0. ) ) { ++sb_n_bad; continue; }
    ++sb_n_trk;
    bool cont = pcv.is_inside( track_end_x->at( p ), track_end_y->at( p ),
      track_end_z->at( p ) );
    float hip = ng2_ok( p ) ? ng2_hipfrac->at( p ) : -1.f;
    int sem = ng2_ok( p ) ? ng2_semlabel->at( p ) : BOGUS_INDEX;
    float KEp = track_KE_proton->at( p );
    float p_mom = real_sqrt( KEp*KEp + 2.*PROTON_MASS*KEp );
    float sd = track_start_distance->at( p );
    bool proton_like = cont && hip >= PROTON_NG2_HIP_FRAC_MIN
      && pfp_n_desc->at( p ) == 0u;
    if ( proton_like && p_mom >= PROTON_P_MIN && p_mom <= PROTON_P_MAX ) {
      sb_good.emplace_back( p_mom, p );
      if ( std::isfinite( sd ) ) {
        sb_sd_min = std::min( sb_sd_min, sd );
        sb_sd_max = std::max( sb_sd_max, sd );
      }
    }
    if ( hip < PROTON_NG2_HIP_FRAC_MIN && sem == NG2_MIP_LABEL
      && sd < PROTON_VTX_DIST_NEAR ) sb_pion.emplace_back( p_mom, p );
  }
  int sb_n_good = sb_good.size(), sb_n_pion = sb_pion.size();
  bool sb_clean = sb_muon_ok && sb_n_bad == 0 && sb_n_shw == 0
    && sb_sd_min < PROTON_VTX_DIST_NEAR && sb_sd_max < PROTON_VTX_DIST_FAR;
  bool sb_3p = sb_clean && sb_n_trk == 3 && sb_n_good == 3;
  bool sb_npi = sb_clean && sb_n_pion == 1
    && ( ( sb_n_trk == 2 && sb_n_good == 1 ) || ( sb_n_trk == 3 && sb_n_good == 2 ) );
  bool sb_1p = sb_clean && sb_n_trk == 1 && sb_n_good == 1;

  auto by_mom = []( const auto& a, const auto& b ) { return a.first > b.first; };
  std::vector< std::pair< float, int > > sb_prot = sb_good;
  if ( sb_npi && sb_n_good == 1 ) sb_prot.push_back( sb_pion.front() );
  std::stable_sort( sb_prot.begin(), sb_prot.end(), by_mom );
  TVector3 sb_p3[ 2 ];
  for ( size_t r = 0; r < 2u && r < sb_prot.size(); ++r ) {
    sb_p3[ r ] = get_reco_p3_proton( sb_prot.at( r ).second, sb_prot.at( r ).first );
  }
  bool sb_two = muon && ( sb_3p || sb_npi );
  fill_observables( out, sb_two, p3mu, sb_p3[ 0 ], sb_p3[ 1 ], "sb_",
    muon && sb_1p );
  out[ "sb_topo_3p" ] = sb_3p;
  out[ "sb_topo_npi" ] = sb_npi;
  out[ "sb_topo_1p" ] = sb_1p;
  out[ "sb_n_good_p" ] = sb_n_good;
  out[ "sb_n_pion" ] = sb_n_pion;

  // We're done, store the results in the output tree and return
  // the boolean flag indicating whether the full selection was passed or not
  out[ "reco_vertex_in_FV" ] = sel_reco_vertex_in_FV;
  out[ "pfp_starts_in_PCV" ] = sel_pfp_starts_in_PCV;
  out[ "has_muon_candidate" ] = sel_has_muon_candidate;
  out[ "topo_cut_passed" ] = sel_topo_cut_passed;
  out[ "nu_mu_cc" ] = sel_nu_mu_cc;
  out[ "muon_contained" ] = sel_muon_contained;
  out[ "muon_passed_mom_cuts" ] = sel_muon_passed_mom_cuts;
  out[ "no_reco_showers" ] = sel_no_reco_showers;
  out[ "has_p_candidate" ] = sel_has_p_candidate;
  out[ "two_p_candidates" ] = sel_two_p_candidates;
  out[ "no_extra_pfps" ] = sel_no_extra_pfps;
  out[ "muon_quality_ok" ] = sel_muon_quality_ok;
  out[ "protons_contained" ] = sel_protons_contained;
  out[ "passed_proton_pid_cut" ] = sel_passed_proton_pid_cut;
  out[ "protons_passed_mom_cuts" ] = sel_protons_passed_mom_cuts;
  out[ "cosmic_ip_cut_passed" ] = sel_cosmic_ip_cut_passed;

  out[ "num_gen2_pfps" ] = num_gen2_pfps;
  out[ "num_p_candidates" ] = static_cast< int >( num_p_candidates );

  out[ "reco_p_mu_range" ] = p_mu_range;
  out[ "reco_p_mu_mcs" ] = p_mu_mcs;
  out[ "reco_mu_frac_diff_range_mcs" ] = frac_diff_range_mcs;

  out[ "reco_lead_p_start_dist" ] = lead_p_start_dist;
  out[ "reco_sub_p_start_dist" ] = sub_p_start_dist;

  fill_observables( out, muon && lead_p && sub_p, p3mu, p3_lead_p, p3_sub_p );
  out[ "reco_p3_p_vec" ] = p3_p_vec;

  // Keep the Pandora results for derived selections
  float raw_vx;
  in.at( "reco_nu_vtx_x" ) >> raw_vx;
  pd_ = PandoraResult();
  pd_.selected = sel_CC2p0pi;
  pd_.nu_mu_cc = sel_nu_mu_cc;
  pd_.reco_vertex_in_FV = sel_reco_vertex_in_FV;
  pd_.has_vertex = std::isfinite( nu_vx ) && std::isfinite( nu_vy )
    && std::isfinite( nu_vz ) && nu_vx > -900. && raw_vx > -900.;
  pd_.has_muon_candidate = sel_has_muon_candidate;
  pd_.muon_contained = sel_muon_contained;
  pd_.muon_passed_mom_cuts = sel_muon_passed_mom_cuts;
  pd_.muon_quality_if_contained = sel_muon_quality_if_contained;
  pd_.muon_ng2_mip = sel_muon_ng2_mip;
  pd_.no_reco_showers = sel_no_reco_showers;
  pd_.two_p_candidates = sel_two_p_candidates;
  pd_.no_extra_pfps = sel_no_extra_pfps;
  pd_.protons_contained = sel_protons_contained;
  pd_.protons_near_vertex = sel_protons_near_vertex;
  pd_.protons_no_descendents = sel_protons_no_descendents;
  if ( ng2_ok( lead_p_candidate_idx ) ) {
    pd_.lead_p_ng2_hipfrac = ng2_hipfrac->at( lead_p_candidate_idx );
  }
  if ( ng2_ok( sub_p_candidate_idx ) ) {
    pd_.sub_p_ng2_hipfrac = ng2_hipfrac->at( sub_p_candidate_idx );
  }
  pd_.p_mu_range = p_mu_range;
  pd_.p_mu_mcs = p_mu_mcs;
  pd_.has_p3 = muon && lead_p && sub_p;
  pd_.has_lead_p = lead_p;
  pd_.has_sub_p = sub_p;
  pd_.p3mu = p3mu;
  pd_.p3_lead_p = p3_lead_p;
  pd_.p3_sub_p = p3_sub_p;
  if ( pd_.has_p3 ) pd_.cos_gamma_pp = std::cos( p3_lead_p.Angle( p3_sub_p ) );
  pd_.sb_3p = sb_3p;
  pd_.sb_npi = sb_npi;
  pd_.sb_1p = sb_1p;

  return sel_CC2p0pi;
}

void CC1mu2p0piv2::fill_observables( TreeMap::Wrapper& out, bool ok,
  const TVector3& p3mu, const TVector3& p3_lead_p, const TVector3& p3_sub_p,
  const std::string& prefix, bool lead_only_ok ) const
{
  // Compute reco observables if we have a muon candidate and both a leading
  // and a subleading proton candidate in the event

  // General kinematics
  double p_mu = BOGUS, cos_theta_mu = BOGUS, p_lead_p = BOGUS,
    cos_theta_lead_p = BOGUS, p_sub_p = BOGUS, cos_theta_sub_p = BOGUS,
    cos_theta_mu_lead_p = BOGUS, cos_theta_mu_sub_p = BOGUS;

  // Proton pair
  double cos_gamma_pp = BOGUS, p_pp = BOGUS, cos_theta_pp = BOGUS,
    theta_mu_pp = BOGUS, p_rel = BOGUS, p_rel_cm = BOGUS, m_pp = BOGUS,
    T_pp = BOGUS, T_asym_pp = BOGUS;

  // Kinematic imbalance (muon + proton pair)
  double delta_pT = BOGUS, delta_phiT = BOGUS, delta_alphaT = BOGUS,
    delta_pTx = BOGUS, delta_pTy = BOGUS, delta_pL = BOGUS, pn = BOGUS,
    p_miss_minus = BOGUS;

  // Kinematic imbalance (muon + leading proton only)
  double lead_delta_pT = BOGUS, lead_delta_phiT = BOGUS,
    lead_delta_alphaT = BOGUS, lead_delta_pTx = BOGUS, lead_delta_pTy = BOGUS,
    lead_delta_pL = BOGUS, lead_pn = BOGUS;

  // Calorimetric and momentum-transfer observables
  double Ecal = BOGUS, E_QE = BOGUS, q0 = BOGUS, q3 = BOGUS, Q2 = BOGUS,
    p_init_lead = BOGUS, cos_theta_pinit_sub = BOGUS,
    cos_theta_lead_q = BOGUS, cos_theta_sub_q = BOGUS,
    cos_theta_pp_q = BOGUS;

  if ( ok ) {

    double E_mu = real_sqrt( p3mu.Mag2() + MUON_MASS*MUON_MASS );
    double E_lead_p = real_sqrt( p3_lead_p.Mag2() + PROTON_MASS*PROTON_MASS );
    double E_sub_p = real_sqrt( p3_sub_p.Mag2() + PROTON_MASS*PROTON_MASS );

    double T_lead_p = E_lead_p - PROTON_MASS;
    double T_sub_p = E_sub_p - PROTON_MASS;

    TVector3 p3_pp = p3_lead_p + p3_sub_p;
    TLorentzVector p4_pp( p3_pp, E_lead_p + E_sub_p );

    // Hadronic invariant mass of the proton pair
    double W = p4_pp.M();

    // General kinematics
    p_mu = p3mu.Mag();
    cos_theta_mu = p3mu.CosTheta();
    p_lead_p = p3_lead_p.Mag();
    cos_theta_lead_p = p3_lead_p.CosTheta();
    p_sub_p = p3_sub_p.Mag();
    cos_theta_sub_p = p3_sub_p.CosTheta();
    cos_theta_mu_lead_p = std::cos( p3mu.Angle( p3_lead_p ) );
    cos_theta_mu_sub_p = std::cos( p3mu.Angle( p3_sub_p ) );

    // Proton pair
    cos_gamma_pp = std::cos( p3_lead_p.Angle( p3_sub_p ) );
    p_pp = p3_pp.Mag();
    cos_theta_pp = p3_pp.CosTheta();
    p_rel = 0.5 * ( p3_lead_p - p3_sub_p ).Mag();
    p_rel_cm = real_sqrt( 0.25*W*W - PROTON_MASS*PROTON_MASS );
    m_pp = W;
    T_pp = T_lead_p + T_sub_p;
    if ( T_pp > 0. ) T_asym_pp = ( T_lead_p - T_sub_p ) / T_pp;

    // Kinematic imbalance using the proton pair as the hadronic system.
    // Angles from the KICalculator are in degrees.
    KICalculator ki_calc_pp( p3mu, p3_pp, ki_calc_type_, MUON_MASS, W );

    delta_pT = ki_calc_pp.pT();
    delta_phiT = ki_calc_pp.delta_phiT();
    delta_alphaT = ki_calc_pp.delta_alphaT();
    delta_pTx = ki_calc_pp.pTx();
    delta_pTy = ki_calc_pp.pTy();
    theta_mu_pp = ki_calc_pp.theta_lep_had();
    p_miss_minus = ki_calc_pp.p_miss_minus();
    E_QE = ki_calc_pp.E_QE();

    // NOTE: The KICalculator computes the hadronic kinetic energy as
    // E_had - m_had. For the proton pair this is E_lead + E_sub - W, which is
    // smaller than T_lead + T_sub by W - 2*m_p. Correct for this here to
    // obtain the calorimetric neutrino energy.
    Ecal = ki_calc_pp.Ecal() + W - 2.*PROTON_MASS;

    // NOTE: KICalculator::pL() and the quantities derived from it (pn,
    // pn_perp, pn_par, delta_alpha3D, ...) are not used for the proton pair.
    // The kOpt1-kOpt3 expression assumes an A-1 residual nucleus, which
    // biases pL by roughly -m_p for two-proton knockout, and the kOpt4
    // expression relies on the uncorrected calorimetric energy. The kOpt4
    // definition is evaluated here with the corrected energy instead. With
    // this definition, pn is the magnitude of the proton pair's
    // center-of-mass momentum, p3_pp - q.
    delta_pL = p3mu.Z() + p3_pp.Z() - Ecal;
    pn = real_sqrt( delta_pT*delta_pT + delta_pL*delta_pL );

    // Kinematic imbalance using only the leading proton (i.e., what
    // CC1muNp0pi computes for the same event)
    KICalculator ki_calc_lead( p3mu, p3_lead_p, ki_calc_type_ );

    lead_delta_pT = ki_calc_lead.pT();
    lead_delta_phiT = ki_calc_lead.delta_phiT();
    lead_delta_alphaT = ki_calc_lead.delta_alphaT();
    lead_delta_pTx = ki_calc_lead.pTx();
    lead_delta_pTy = ki_calc_lead.pTy();
    lead_delta_pL = ki_calc_lead.pL();
    lead_pn = ki_calc_lead.pn();

    // Energy and momentum transfer, assuming a neutrino along +z with the
    // calorimetric energy computed above
    TVector3 q3_vec( -p3mu.X(), -p3mu.Y(), Ecal - p3mu.Z() );

    q0 = Ecal - E_mu;
    q3 = q3_vec.Mag();
    Q2 = q3*q3 - q0*q0;

    // Initial momentum of the struck nucleon, assuming that it is the
    // leading proton, and its angle with respect to the subleading proton
    // (back-to-back for a correlated pair)
    TVector3 p3_init_lead = p3_lead_p - q3_vec;

    p_init_lead = p3_init_lead.Mag();
    cos_theta_pinit_sub = std::cos( p3_init_lead.Angle( p3_sub_p ) );

    cos_theta_lead_q = std::cos( p3_lead_p.Angle( q3_vec ) );
    cos_theta_sub_q = std::cos( p3_sub_p.Angle( q3_vec ) );
    cos_theta_pp_q = std::cos( p3_pp.Angle( q3_vec ) );
  }
  else if ( lead_only_ok ) {
    // Muon and leading proton only (e.g., a one-proton sideband)
    p_mu = p3mu.Mag();
    cos_theta_mu = p3mu.CosTheta();
    p_lead_p = p3_lead_p.Mag();
    cos_theta_lead_p = p3_lead_p.CosTheta();
    cos_theta_mu_lead_p = std::cos( p3mu.Angle( p3_lead_p ) );

    KICalculator ki_calc_lead( p3mu, p3_lead_p, ki_calc_type_ );
    lead_delta_pT = ki_calc_lead.pT();
    lead_delta_phiT = ki_calc_lead.delta_phiT();
    lead_delta_alphaT = ki_calc_lead.delta_alphaT();
    lead_delta_pTx = ki_calc_lead.pTx();
    lead_delta_pTy = ki_calc_lead.pTy();
    lead_delta_pL = ki_calc_lead.pL();
    lead_pn = ki_calc_lead.pn();
  }

  out[ prefix + "p_mu" ] = p_mu;
  out[ prefix + "cos_theta_mu" ] = cos_theta_mu;
  out[ prefix + "p_lead_p" ] = p_lead_p;
  out[ prefix + "cos_theta_lead_p" ] = cos_theta_lead_p;
  out[ prefix + "p_sub_p" ] = p_sub_p;
  out[ prefix + "cos_theta_sub_p" ] = cos_theta_sub_p;
  out[ prefix + "cos_theta_mu_lead_p" ] = cos_theta_mu_lead_p;
  out[ prefix + "cos_theta_mu_sub_p" ] = cos_theta_mu_sub_p;

  out[ prefix + "cos_gamma_pp" ] = cos_gamma_pp;
  out[ prefix + "p_pp" ] = p_pp;
  out[ prefix + "cos_theta_pp" ] = cos_theta_pp;
  out[ prefix + "theta_mu_pp" ] = theta_mu_pp;
  out[ prefix + "p_rel" ] = p_rel;
  out[ prefix + "p_rel_cm" ] = p_rel_cm;
  out[ prefix + "m_pp" ] = m_pp;
  out[ prefix + "T_pp" ] = T_pp;
  out[ prefix + "T_asym_pp" ] = T_asym_pp;

  out[ prefix + "delta_pT" ] = delta_pT;
  out[ prefix + "delta_phiT" ] = delta_phiT;
  out[ prefix + "delta_alphaT" ] = delta_alphaT;
  out[ prefix + "delta_pTx" ] = delta_pTx;
  out[ prefix + "delta_pTy" ] = delta_pTy;
  out[ prefix + "delta_pL" ] = delta_pL;
  out[ prefix + "pn" ] = pn;
  out[ prefix + "p_miss_minus" ] = p_miss_minus;

  out[ prefix + "lead_delta_pT" ] = lead_delta_pT;
  out[ prefix + "lead_delta_phiT" ] = lead_delta_phiT;
  out[ prefix + "lead_delta_alphaT" ] = lead_delta_alphaT;
  out[ prefix + "lead_delta_pTx" ] = lead_delta_pTx;
  out[ prefix + "lead_delta_pTy" ] = lead_delta_pTy;
  out[ prefix + "lead_delta_pL" ] = lead_delta_pL;
  out[ prefix + "lead_pn" ] = lead_pn;

  out[ prefix + "Ecal" ] = Ecal;
  out[ prefix + "E_QE" ] = E_QE;
  out[ prefix + "q0" ] = q0;
  out[ prefix + "q3" ] = q3;
  out[ prefix + "Q2" ] = Q2;
  out[ prefix + "p_init_lead" ] = p_init_lead;
  out[ prefix + "cos_theta_pinit_sub" ] = cos_theta_pinit_sub;
  out[ prefix + "cos_theta_lead_q" ] = cos_theta_lead_q;
  out[ prefix + "cos_theta_sub_q" ] = cos_theta_sub_q;
  out[ prefix + "cos_theta_pp_q" ] = cos_theta_pp_q;

  out[ prefix + "p3_mu" ] = p3mu;
  out[ prefix + "p3_lead_p" ] = p3_lead_p;
  out[ prefix + "p3_sub_p" ] = p3_sub_p;
}

std::string CC1mu2p0piv2::categorize_event( AnalysisEvent& ev ) {

  // Get access to the input and output trees
  const auto& in = ev.in( "nuselection/NeutrinoSelectionFilter" );
  auto& out = ev.out();

  // Retrieve some truth information from the PeLEE ntuple
  float mc_nu_vx, mc_nu_vy, mc_nu_vz;
  int mc_nu_pdg, mc_ccnc, mc_interaction_type;
  in.at( "true_nu_vtx_x" ) >> mc_nu_vx;
  in.at( "true_nu_vtx_y" ) >> mc_nu_vy;
  in.at( "true_nu_vtx_z" ) >> mc_nu_vz;
  in.at( "nu_pdg" ) >> mc_nu_pdg;
  in.at( "ccnc" ) >> mc_ccnc;
  in.at( "interaction" ) >> mc_interaction_type;

  // Check whether the true neutrino vertex lies within the fiducial volume
  bool sig_in_fv = this->get_fv().is_inside( mc_nu_vx, mc_nu_vy, mc_nu_vz );

  bool sig_is_numu = ( mc_nu_pdg == MUON_NEUTRINO );
  bool is_nc = ( mc_ccnc == NEUTRAL_CURRENT );

  // Set some flags to their default values before checking them
  bool sig_no_fs_mesons = true;
  bool sig_mc_no_fs_pi0 = true;
  bool sig_mc_no_charged_pi_above_threshold = true;
  bool sig_mu_in_mom_range = false;
  bool sig_lead_p_mom_in_range = false;
  bool sig_one_p_in_mom_range = false;
  bool sig_two_p_in_mom_range = false;

  std::vector< int >* nu_daughter_pdg;
  std::vector< float > *nu_daughter_energy, *nu_daughter_px,
    *nu_daughter_py, *nu_daughter_pz;

  in.at( "mc_pdg" ) >> nu_daughter_pdg;
  in.at( "mc_E" ) >> nu_daughter_energy;
  in.at( "mc_px" ) >> nu_daughter_px;
  in.at( "mc_py" ) >> nu_daughter_py;
  in.at( "mc_pz" ) >> nu_daughter_pz;

  int num_p_in_mom_range = 0;
  double mom_lead_p = 0.;

  for ( size_t p = 0u; p < nu_daughter_pdg->size(); ++p ) {
    int pdg = nu_daughter_pdg->at( p );
    float energy = nu_daughter_energy->at( p );

    // Do the general check for (anti)mesons first before considering
    // any individual PDG codes
    if ( is_meson_or_antimeson( pdg ) ) {
      sig_no_fs_mesons = false;
    }

    if ( pdg == MUON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(MUON_MASS, 2) );
      if ( mom >= MUON_P_MIN ) {
        sig_mu_in_mom_range = true;
      }
    }
    else if ( pdg == PROTON ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PROTON_MASS, 2) );
      if ( mom > mom_lead_p ) {
        mom_lead_p = mom;
      }
      if ( mom >= PROTON_P_MIN && mom <= PROTON_P_MAX ) {
        ++num_p_in_mom_range;
      }
    }
    // Not part of the signal definition, but the flag is recorded anyway
    else if ( pdg == PI_ZERO ) {
      sig_mc_no_fs_pi0 = false;
    }
    // Also not part of the signal definition
    else if ( std::abs(pdg) == PI_PLUS ) {
      double mom = real_sqrt( std::pow(energy, 2) - std::pow(PI_PLUS_MASS, 2) );
      if ( mom > PI_PLUS_P_MIN ) {
        sig_mc_no_charged_pi_above_threshold = false;
      }
    }
  }

  // Check that the leading proton has a momentum within the allowed range.
  // Since it is the most energetic proton, this also ensures that there are
  // no protons above the upper momentum limit.
  if ( mom_lead_p >= PROTON_P_MIN && mom_lead_p <= PROTON_P_MAX ) {
    sig_lead_p_mom_in_range = true;
  }

  // Check number of protons in momentum range
  sig_one_p_in_mom_range = ( num_p_in_mom_range == 1 );
  sig_two_p_in_mom_range = ( num_p_in_mom_range == 2 );

  // Mesonless numu CC event in the fiducial volume with a muon above
  // threshold and any number of protons
  bool is_cc0pi = sig_in_fv && !is_nc && sig_is_numu && sig_mu_in_mom_range
    && sig_no_fs_mesons;

  bool is_signal = is_cc0pi && sig_lead_p_mom_in_range
    && sig_two_p_in_mom_range;

  // Set the true 3-momentum of the final-state muon if there is one
  size_t num_mc_daughters = nu_daughter_pdg->size();

  TVector3 mc_p3mu;
  bool true_muon = ( sig_is_numu && mc_ccnc == CHARGED_CURRENT );
  if ( true_muon ) {
    // Loop over the MC neutrino daughters, find the muon, and get its
    // true 3-momentum. Note that we assume there is only one muon in
    // this loop.
    bool found_muon = false;
    for ( size_t d = 0u; d < num_mc_daughters; ++d ) {
      int pdg = nu_daughter_pdg->at( d );
      if ( pdg == MUON ) {
        found_muon = true;
        float px = nu_daughter_px->at( d );
        float py = nu_daughter_py->at( d );
        float pz = nu_daughter_pz->at( d );
        mc_p3mu = TVector3( px, py, pz );
        break;
      }
    }

    if ( !found_muon ) {
      std::cout << "WARNING: Missing muon in true numu CC event!\n";
    }
  }

  // Set the true 3-momentum of the final-state protons (if there are any)
  std::vector< TVector3 > mc_p3_p_vec;

  for ( size_t p = 0u; p < num_mc_daughters; ++p ) {
    int pdg = nu_daughter_pdg->at( p );
    if ( pdg == PROTON )
    {
      float px = nu_daughter_px->at( p );
      float py = nu_daughter_py->at( p );
      float pz = nu_daughter_pz->at( p );
      TVector3 temp_p3 = TVector3( px, py, pz );

      mc_p3_p_vec.push_back( temp_p3 );
    }
  }

  // Sort the true proton 3-momenta in order from highest to lowest magnitude
  std::sort( mc_p3_p_vec.begin(), mc_p3_p_vec.end(), [](const TVector3& a,
    const TVector3& b) -> bool { return a.Mag() > b.Mag(); } );

  // Get the leading and subleading proton specifically. For signal events,
  // these are the two protons within the allowed momentum range (any other
  // protons lie below the lower momentum limit).
  TVector3 mc_p3_lead_p;
  TVector3 mc_p3_sub_p;

  bool true_lead_p = ( mc_p3_p_vec.size() > 0u );
  bool true_sub_p = ( mc_p3_p_vec.size() > 1u );

  if ( true_lead_p ) mc_p3_lead_p = mc_p3_p_vec.at( 0 );
  if ( true_sub_p ) mc_p3_sub_p = mc_p3_p_vec.at( 1 );

  // Warn the user if two protons were not found for a signal event
  if ( !true_sub_p && is_signal ) {
    std::cout << "WARNING: Missing proton pair in MC signal event!\n";
  }

  // Compute true observables if the event contains a muon and both a leading
  // and a subleading proton

  // General kinematics
  double mc_p_mu = BOGUS, mc_cos_theta_mu = BOGUS, mc_p_lead_p = BOGUS,
    mc_cos_theta_lead_p = BOGUS, mc_p_sub_p = BOGUS,
    mc_cos_theta_sub_p = BOGUS, mc_cos_theta_mu_lead_p = BOGUS,
    mc_cos_theta_mu_sub_p = BOGUS;

  // Proton pair
  double mc_cos_gamma_pp = BOGUS, mc_p_pp = BOGUS, mc_cos_theta_pp = BOGUS,
    mc_theta_mu_pp = BOGUS, mc_p_rel = BOGUS, mc_p_rel_cm = BOGUS,
    mc_m_pp = BOGUS, mc_T_pp = BOGUS, mc_T_asym_pp = BOGUS;

  // Kinematic imbalance (muon + proton pair)
  double mc_delta_pT = BOGUS, mc_delta_phiT = BOGUS, mc_delta_alphaT = BOGUS,
    mc_delta_pTx = BOGUS, mc_delta_pTy = BOGUS, mc_delta_pL = BOGUS,
    mc_pn = BOGUS, mc_p_miss_minus = BOGUS;

  // Kinematic imbalance (muon + leading proton only)
  double mc_lead_delta_pT = BOGUS, mc_lead_delta_phiT = BOGUS,
    mc_lead_delta_alphaT = BOGUS, mc_lead_delta_pTx = BOGUS,
    mc_lead_delta_pTy = BOGUS, mc_lead_delta_pL = BOGUS, mc_lead_pn = BOGUS;

  // Calorimetric and momentum-transfer observables
  double mc_Ecal = BOGUS, mc_E_QE = BOGUS, mc_q0 = BOGUS, mc_q3 = BOGUS,
    mc_Q2 = BOGUS, mc_p_init_lead = BOGUS, mc_cos_theta_pinit_sub = BOGUS,
    mc_cos_theta_lead_q = BOGUS, mc_cos_theta_sub_q = BOGUS,
    mc_cos_theta_pp_q = BOGUS;

  if ( true_muon && true_lead_p && true_sub_p ) {

    double E_mu = real_sqrt( mc_p3mu.Mag2() + MUON_MASS*MUON_MASS );
    double E_lead_p = real_sqrt( mc_p3_lead_p.Mag2()
      + PROTON_MASS*PROTON_MASS );
    double E_sub_p = real_sqrt( mc_p3_sub_p.Mag2()
      + PROTON_MASS*PROTON_MASS );

    double T_lead_p = E_lead_p - PROTON_MASS;
    double T_sub_p = E_sub_p - PROTON_MASS;

    TVector3 mc_p3_pp = mc_p3_lead_p + mc_p3_sub_p;
    TLorentzVector mc_p4_pp( mc_p3_pp, E_lead_p + E_sub_p );

    // Hadronic invariant mass of the proton pair
    double W = mc_p4_pp.M();

    // General kinematics
    mc_p_mu = mc_p3mu.Mag();
    mc_cos_theta_mu = mc_p3mu.CosTheta();
    mc_p_lead_p = mc_p3_lead_p.Mag();
    mc_cos_theta_lead_p = mc_p3_lead_p.CosTheta();
    mc_p_sub_p = mc_p3_sub_p.Mag();
    mc_cos_theta_sub_p = mc_p3_sub_p.CosTheta();
    mc_cos_theta_mu_lead_p = std::cos( mc_p3mu.Angle( mc_p3_lead_p ) );
    mc_cos_theta_mu_sub_p = std::cos( mc_p3mu.Angle( mc_p3_sub_p ) );

    // Proton pair
    mc_cos_gamma_pp = std::cos( mc_p3_lead_p.Angle( mc_p3_sub_p ) );
    mc_p_pp = mc_p3_pp.Mag();
    mc_cos_theta_pp = mc_p3_pp.CosTheta();
    mc_p_rel = 0.5 * ( mc_p3_lead_p - mc_p3_sub_p ).Mag();
    mc_p_rel_cm = real_sqrt( 0.25*W*W - PROTON_MASS*PROTON_MASS );
    mc_m_pp = W;
    mc_T_pp = T_lead_p + T_sub_p;
    if ( mc_T_pp > 0. ) mc_T_asym_pp = ( T_lead_p - T_sub_p ) / mc_T_pp;

    // Kinematic imbalance using the proton pair as the hadronic system.
    // Angles from the KICalculator are in degrees.
    KICalculator ki_calc_pp( mc_p3mu, mc_p3_pp, ki_calc_type_, MUON_MASS, W );

    mc_delta_pT = ki_calc_pp.pT();
    mc_delta_phiT = ki_calc_pp.delta_phiT();
    mc_delta_alphaT = ki_calc_pp.delta_alphaT();
    mc_delta_pTx = ki_calc_pp.pTx();
    mc_delta_pTy = ki_calc_pp.pTy();
    mc_theta_mu_pp = ki_calc_pp.theta_lep_had();
    mc_p_miss_minus = ki_calc_pp.p_miss_minus();
    mc_E_QE = ki_calc_pp.E_QE();

    // See the notes in is_selected() about the calorimetric energy and the
    // longitudinal variables for the proton pair
    mc_Ecal = ki_calc_pp.Ecal() + W - 2.*PROTON_MASS;

    mc_delta_pL = mc_p3mu.Z() + mc_p3_pp.Z() - mc_Ecal;
    mc_pn = real_sqrt( mc_delta_pT*mc_delta_pT + mc_delta_pL*mc_delta_pL );

    // Kinematic imbalance using only the leading proton (i.e., what
    // CC1muNp0pi computes for the same event)
    KICalculator ki_calc_lead( mc_p3mu, mc_p3_lead_p, ki_calc_type_ );

    mc_lead_delta_pT = ki_calc_lead.pT();
    mc_lead_delta_phiT = ki_calc_lead.delta_phiT();
    mc_lead_delta_alphaT = ki_calc_lead.delta_alphaT();
    mc_lead_delta_pTx = ki_calc_lead.pTx();
    mc_lead_delta_pTy = ki_calc_lead.pTy();
    mc_lead_delta_pL = ki_calc_lead.pL();
    mc_lead_pn = ki_calc_lead.pn();

    // Energy and momentum transfer, assuming a neutrino along +z with the
    // calorimetric energy computed above
    TVector3 q3_vec( -mc_p3mu.X(), -mc_p3mu.Y(), mc_Ecal - mc_p3mu.Z() );

    mc_q0 = mc_Ecal - E_mu;
    mc_q3 = q3_vec.Mag();
    mc_Q2 = mc_q3*mc_q3 - mc_q0*mc_q0;

    // Initial momentum of the struck nucleon, assuming that it is the
    // leading proton, and its angle with respect to the subleading proton
    // (back-to-back for a correlated pair)
    TVector3 mc_p3_init_lead = mc_p3_lead_p - q3_vec;

    mc_p_init_lead = mc_p3_init_lead.Mag();
    mc_cos_theta_pinit_sub = std::cos( mc_p3_init_lead.Angle( mc_p3_sub_p ) );

    mc_cos_theta_lead_q = std::cos( mc_p3_lead_p.Angle( q3_vec ) );
    mc_cos_theta_sub_q = std::cos( mc_p3_sub_p.Angle( q3_vec ) );
    mc_cos_theta_pp_q = std::cos( mc_p3_pp.Angle( q3_vec ) );
  }

  // We're done. Store MC truth information for this event and return the
  // string indicating the category to which it belongs.
  out[ "mc_is_numu" ] = sig_is_numu;
  out[ "mc_vertex_in_FV" ] = sig_in_fv;
  out[ "mc_lead_p_in_range" ] = sig_lead_p_mom_in_range;
  out[ "mc_one_p_in_mom_range" ] = sig_one_p_in_mom_range;
  out[ "mc_two_p_in_mom_range" ] = sig_two_p_in_mom_range;
  out[ "mc_muon_in_mom_range" ] = sig_mu_in_mom_range;
  out[ "mc_no_FS_mesons" ] = sig_no_fs_mesons;
  out[ "mc_no_charged_pions_above_threshold" ]
    = sig_mc_no_charged_pi_above_threshold;
  out[ "mc_no_pi0s" ] = sig_mc_no_fs_pi0;
  out[ "mc_is_signal" ] = is_signal;
  out[ "num_p_in_mom_range" ] = num_p_in_mom_range;

  out[ "true_p_mu" ] = mc_p_mu;
  out[ "true_cos_theta_mu" ] = mc_cos_theta_mu;
  out[ "true_p_lead_p" ] = mc_p_lead_p;
  out[ "true_cos_theta_lead_p" ] = mc_cos_theta_lead_p;
  out[ "true_p_sub_p" ] = mc_p_sub_p;
  out[ "true_cos_theta_sub_p" ] = mc_cos_theta_sub_p;
  out[ "true_cos_theta_mu_lead_p" ] = mc_cos_theta_mu_lead_p;
  out[ "true_cos_theta_mu_sub_p" ] = mc_cos_theta_mu_sub_p;

  out[ "true_cos_gamma_pp" ] = mc_cos_gamma_pp;
  out[ "true_p_pp" ] = mc_p_pp;
  out[ "true_cos_theta_pp" ] = mc_cos_theta_pp;
  out[ "true_theta_mu_pp" ] = mc_theta_mu_pp;
  out[ "true_p_rel" ] = mc_p_rel;
  out[ "true_p_rel_cm" ] = mc_p_rel_cm;
  out[ "true_m_pp" ] = mc_m_pp;
  out[ "true_T_pp" ] = mc_T_pp;
  out[ "true_T_asym_pp" ] = mc_T_asym_pp;

  out[ "true_delta_pT" ] = mc_delta_pT;
  out[ "true_delta_phiT" ] = mc_delta_phiT;
  out[ "true_delta_alphaT" ] = mc_delta_alphaT;
  out[ "true_delta_pTx" ] = mc_delta_pTx;
  out[ "true_delta_pTy" ] = mc_delta_pTy;
  out[ "true_delta_pL" ] = mc_delta_pL;
  out[ "true_pn" ] = mc_pn;
  out[ "true_p_miss_minus" ] = mc_p_miss_minus;

  out[ "true_lead_delta_pT" ] = mc_lead_delta_pT;
  out[ "true_lead_delta_phiT" ] = mc_lead_delta_phiT;
  out[ "true_lead_delta_alphaT" ] = mc_lead_delta_alphaT;
  out[ "true_lead_delta_pTx" ] = mc_lead_delta_pTx;
  out[ "true_lead_delta_pTy" ] = mc_lead_delta_pTy;
  out[ "true_lead_delta_pL" ] = mc_lead_delta_pL;
  out[ "true_lead_pn" ] = mc_lead_pn;

  out[ "true_Ecal" ] = mc_Ecal;
  out[ "true_E_QE" ] = mc_E_QE;
  out[ "true_q0" ] = mc_q0;
  out[ "true_q3" ] = mc_q3;
  out[ "true_Q2" ] = mc_Q2;
  out[ "true_p_init_lead" ] = mc_p_init_lead;
  out[ "true_cos_theta_pinit_sub" ] = mc_cos_theta_pinit_sub;
  out[ "true_cos_theta_lead_q" ] = mc_cos_theta_lead_q;
  out[ "true_cos_theta_sub_q" ] = mc_cos_theta_sub_q;
  out[ "true_cos_theta_pp_q" ] = mc_cos_theta_pp_q;

  out[ "true_p3_mu" ] = mc_p3mu;
  out[ "true_p3_lead_p" ] = mc_p3_lead_p;
  out[ "true_p3_sub_p" ] = mc_p3_sub_p;
  out[ "true_p3_p_vec" ] = mc_p3_p_vec;

  // The categories blend those of CC1muNp0pi (signal split by interaction
  // mode, pion and "other" backgrounds) with those of CC1mu2p0pi (mesonless
  // backgrounds split by proton multiplicity)
  if ( !sig_in_fv ) return "Out FV";
  else if ( is_nc ) return "NC";
  else if ( mc_nu_pdg == ELECTRON_NEUTRINO ) return "#nu_{e} CC";
  else if ( mc_nu_pdg != MUON_NEUTRINO ) return "Other";
  else if ( is_signal ) {
    switch ( mc_interaction_type ) {
      case QE_INTERACTION: return "Signal (CCQE)"; break;
      case MEC_INTERACTION: return "Signal (CCMEC)"; break;
      case RES_INTERACTION: return "Signal (CCRES)"; break;
      default: return "Signal (CC other)";
    }
  }
  // Mesonless numu CC events with the wrong proton multiplicity. Events
  // with a proton above the upper momentum limit are left for the
  // "Other #nu_{#mu} CC" category below.
  else if ( is_cc0pi && mom_lead_p <= PROTON_P_MAX ) {
    if ( num_p_in_mom_range == 0 ) return "CC0#pi0p";
    else if ( sig_one_p_in_mom_range ) return "CC0#pi1p";
    else return "CC0#pi(M>2)p";
  }
  else if ( !sig_mc_no_fs_pi0 || !sig_mc_no_charged_pi_above_threshold ) {
    return "#nu_{#mu} CCN#pi";
  }
  else return "Other #nu_{#mu} CC";
}