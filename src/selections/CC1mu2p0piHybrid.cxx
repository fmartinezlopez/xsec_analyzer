// Standard library includes
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// XSecAnalyzer includes
#include "XSecAnalyzer/Array.hh"
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/Selections/CC1mu2p0piHybrid.hh"

namespace {

  const std::string PELEE_TREE = "nuselection/NeutrinoSelectionFilter";
  const std::string WC_PF_TREE = "wcpselection/T_PFeval";
  const std::string WC_BDT_TREE = "wcpselection/T_BDTvars";
  const std::string WC_EVAL_TREE = "wcpselection/T_eval";
  const std::string LANTERN_TREE = "lantern/EventTree";

  constexpr int ELECTRON_PDG = 11;
  constexpr int PHOTON_PDG = 22;

  // Phase space of the measurement (also used by CC1mu2p0piv2)
  constexpr double MUON_P_MIN = 0.10; // GeV/c
  constexpr double PROTON_P_MIN = 0.30; // GeV/c
  constexpr double PROTON_P_MAX = 1.00; // GeV/c

  // Wire-Cell and LANTERN positions are not corrected for space charge. Their
  // vertices are compared with the fiducial volume after this average shift
  // in x (Pandora raw -> SCE-corrected vertex for selected events).
  constexpr double RAW_TO_SCE_DX = -1.82; // cm

  // ---- Pandora vote: CC1mu2p0piv2 without the attachment, descendant and
  // NuGraph2 muon requirements, and a looser NuGraph2 proton HIP fraction
  constexpr double P_VOTE_HIP_FRAC_MIN = 0.4;

  // ---- Wire-Cell vote: numu_score above this value, a primary muon with
  // p >= MUON_P_MIN, exactly two primary protons in the proton window (the
  // two hardest in [PROTON_P_MIN, PROTON_P_MAX], any third below
  // PROTON_P_MIN), no primary charged pion, Wire-Cell vertex in the
  // fiducial volume
  constexpr double WC_NUMU_SCORE_MIN = 0.0;

  // ---- LANTERN vote: exactly one muon-like and two proton-like primary
  // tracks within L_ATTACH_DIST of the vertex (after merging duplicate
  // prongs), no pion-like track, no primary e/gamma shower above
  // L_SHOWER_E_MAX, vertex in the fiducial volume. Protons are counted in
  // the proton window on the LANTERN (trackRecoE) momentum, as are the
  // assigned kinematics
  constexpr double L_ATTACH_DIST = 4.0; // cm, nearer track end to the vertex
  constexpr double L_MERGE_ANGLE = 25.; // deg, duplicate prongs
  constexpr double L_MU_PROB_MIN = 0.5; // LArPID muon probability
  constexpr double L_P_PROB_MIN = 0.0; // LArPID proton probability
  constexpr double L_PI_KE_MIN = 0.0; // MeV, pion-like tracks above this veto
  constexpr double L_PI_PROB_MIN = 0.5; // LArPID pion probability
  constexpr double L_SHOWER_E_MAX = 50.; // MeV

  // ---- Kinematic assignment: a particle of another reconstruction is the
  // same particle if its direction is within this angle
  constexpr double ALIGN_MAX_DEG = 10.;

  // ---- Pandora rescue: CC1mu2p0piv2 events with nearly collinear protons
  // (which neither Wire-Cell nor LANTERN can confirm)
  constexpr double RESCUE_COS_GAMMA_MIN = 0.9;

  double ke_of_p( double p, double m ) { return std::sqrt( p*p + m*m ) - m; }
  double p_of_ke( double ke, double m ) {
    ke = std::max( ke, 0. );
    return std::sqrt( ke*ke + 2.*ke*m );
  }

  bool aligned( const TVector3& a, const TVector3& b ) {
    if ( a.Mag() <= 0. || b.Mag() <= 0. ) return false;
    return a.Angle( b ) * 180. / M_PI < ALIGN_MAX_DEG;
  }

  // Element of a C-style array branch. Array< T >::vector() is not used: its
  // stored pointer does not follow the Array when it is moved into the
  // TreeMap, so it can refer to an empty vector.
  template < typename T > T elem( Array< T >* a, size_t i ) {
    if ( i >= a->size() ) throw std::out_of_range( "array index out of range" );
    return a->data()[ i ];
  }

  // Scalar of an input TTree
  template < typename T > T get( const TreeMap::Wrapper& tm, const std::string& br ) {
    T v;
    tm.at( br ) >> v;
    return v;
  }
}

CC1mu2p0piHybrid::CC1mu2p0piHybrid() : CC1mu2p0piv2( "CC1mu2p0piHybrid" ) {}

bool CC1mu2p0piHybrid::is_selected( AnalysisEvent& ev ) {

  // Pandora: the CC1mu2p0piv2 selection (also fills all of its output
  // branches, including the Pandora reco_* observables, and pd_)
  bool pandora_default = CC1mu2p0piv2::is_selected( ev );

  auto& out = ev.out();
  const auto& fv = this->get_fv();

  // ------------------------------------------------------------ Pandora vote
  double p_lead = pd_.has_p3 ? pd_.p3_lead_p.Mag() : BOGUS;
  double p_sub = pd_.has_p3 ? pd_.p3_sub_p.Mag() : BOGUS;
  bool vote_P = pd_.nu_mu_cc && pd_.muon_passed_mom_cuts
    && pd_.two_p_candidates && pd_.no_extra_pfps && pd_.protons_contained
    && p_lead <= PROTON_P_MAX && p_sub >= PROTON_P_MIN && pd_.no_reco_showers
    && pd_.muon_quality_if_contained
    && pd_.lead_p_ng2_hipfrac >= P_VOTE_HIP_FRAC_MIN
    && pd_.sub_p_ng2_hipfrac >= P_VOTE_HIP_FRAC_MIN;
  bool P_ok = pd_.has_muon_candidate && pd_.two_p_candidates && pd_.has_p3
    && pd_.p3mu.Mag() > 0. && p_lead > 0. && p_sub > 0.;

  // ---------------------------------------------------------- Wire-Cell vote
  const auto& pf = ev.in( WC_PF_TREE );
  Array< int > *wc_mother = nullptr, *wc_pdg = nullptr;
  Array< float >* wc_mom = nullptr;
  pf.at( "reco_mother" ) >> wc_mother;
  pf.at( "reco_pdg" ) >> wc_pdg;
  pf.at( "reco_startMomentum" ) >> wc_mom;
  size_t n_wc = wc_pdg->size();

  std::vector< std::pair< double, size_t > > wc_protons;
  int wc_n_pi = 0;
  int wc_mu_idx = -1;
  double wc_mu_E = -1.;
  for ( size_t k = 0; k < n_wc; ++k ) {
    if ( elem( wc_mother, k ) != 0 ) continue; // primary particles only
    int pdg = elem( wc_pdg, k );
    TVector3 p3( ( *wc_mom )( k, 0 ), ( *wc_mom )( k, 1 ), ( *wc_mom )( k, 2 ) );
    if ( pdg == PROTON ) wc_protons.emplace_back( p3.Mag(), k );
    if ( std::abs( pdg ) == PI_PLUS ) ++wc_n_pi;
    if ( std::abs( pdg ) == MUON ) {
      double E = ( *wc_mom )( k, 3 );
      if ( E > wc_mu_E ) { wc_mu_E = E; wc_mu_idx = static_cast< int >( k ); }
    }
  }
  std::stable_sort( wc_protons.begin(), wc_protons.end(),
    []( const auto& a, const auto& b ) { return a.first > b.first; } );
  auto wc_pmag = [&]( size_t r ) {
    return r < wc_protons.size() ? wc_protons.at( r ).first : 0.;
  };
  auto wc_p3 = [&]( size_t r ) {
    size_t k = wc_protons.at( r ).second;
    return TVector3( ( *wc_mom )( k, 0 ), ( *wc_mom )( k, 1 ), ( *wc_mom )( k, 2 ) );
  };
  bool wc_has_mu = wc_mu_idx >= 0 && wc_mu_E > 0.;
  TVector3 wc_p3mu;
  if ( wc_has_mu ) {
    size_t k = static_cast< size_t >( wc_mu_idx );
    wc_p3mu.SetXYZ( ( *wc_mom )( k, 0 ), ( *wc_mom )( k, 1 ), ( *wc_mom )( k, 2 ) );
  }
  double wc_E_mcs = get< double >( pf, "mcs_emu_MCS" );
  double wc_p_mcs = wc_E_mcs > 0.
    ? std::sqrt( std::max( wc_E_mcs*wc_E_mcs - MUON_MASS*MUON_MASS, 0. ) ) : -1.;
  float wc_vx = get< float >( pf, "reco_nuvtxX" );
  float wc_vy = get< float >( pf, "reco_nuvtxY" );
  float wc_vz = get< float >( pf, "reco_nuvtxZ" );
  bool wc_in_fv = fv.is_inside( double( wc_vx ) + RAW_TO_SCE_DX, double( wc_vy ),
    double( wc_vz ) );
  float wc_score = get< float >( ev.in( WC_BDT_TREE ), "numu_score" );
  bool wc_is_fc = get< bool >( ev.in( WC_EVAL_TREE ), "match_isFC" );

  bool vote_W = wc_score > WC_NUMU_SCORE_MIN && wc_has_mu
    && wc_p3mu.Mag() >= MUON_P_MIN && wc_pmag( 1 ) >= PROTON_P_MIN
    && wc_pmag( 0 ) <= PROTON_P_MAX && wc_pmag( 2 ) < PROTON_P_MIN && wc_in_fv
    && wc_n_pi == 0;

  // ------------------------------------------------------------ LANTERN vote
  const auto& lt = ev.in( LANTERN_TREE );
  Array< int > *t_sec = nullptr, *t_pid = nullptr, *s_sec = nullptr, *s_pid = nullptr;
  Array< float > *t_E = nullptr, *t_prs = nullptr, *t_mus = nullptr, *t_pis = nullptr,
    *t_sx = nullptr, *t_sy = nullptr, *t_sz = nullptr, *t_ex = nullptr,
    *t_ey = nullptr, *t_ez = nullptr, *t_dx = nullptr, *t_dy = nullptr,
    *t_dz = nullptr, *s_E = nullptr;
  lt.at( "trackIsSecondary" ) >> t_sec;
  lt.at( "trackPID" ) >> t_pid;
  lt.at( "trackRecoE" ) >> t_E;
  lt.at( "trackPrScore" ) >> t_prs;
  lt.at( "trackMuScore" ) >> t_mus;
  lt.at( "trackPiScore" ) >> t_pis;
  lt.at( "trackStartPosX" ) >> t_sx;
  lt.at( "trackStartPosY" ) >> t_sy;
  lt.at( "trackStartPosZ" ) >> t_sz;
  lt.at( "trackEndPosX" ) >> t_ex;
  lt.at( "trackEndPosY" ) >> t_ey;
  lt.at( "trackEndPosZ" ) >> t_ez;
  lt.at( "trackStartDirX" ) >> t_dx;
  lt.at( "trackStartDirY" ) >> t_dy;
  lt.at( "trackStartDirZ" ) >> t_dz;
  lt.at( "showerIsSecondary" ) >> s_sec;
  lt.at( "showerPID" ) >> s_pid;
  lt.at( "showerRecoE" ) >> s_E;

  int lt_found_vtx = get< int >( lt, "foundVertex" );
  TVector3 lt_vtx( get< float >( lt, "vtxX" ), get< float >( lt, "vtxY" ),
    get< float >( lt, "vtxZ" ) );
  float lt_vscore = get< float >( lt, "vtxScore" );
  float lt_cosfrac = get< float >( lt, "vtxFracHitsOnCosmic" );
  bool lt_in_fv = lt_found_vtx > 0
    && fv.is_inside( lt_vtx.X() + RAW_TO_SCE_DX, lt_vtx.Y(), lt_vtx.Z() );

  size_t n_lt = t_pid->size();
  auto vec = []( Array< float >* x, Array< float >* y, Array< float >* z, size_t i ) {
    return TVector3( elem( x, i ), elem( y, i ), elem( z, i ) );
  };
  std::vector< bool > is_mu( n_lt ), is_p( n_lt ), pi_like( n_lt ), drop( n_lt, false );
  std::vector< double > len( n_lt );
  std::vector< TVector3 > out_dir( n_lt ), kin_dir( n_lt );
  for ( size_t i = 0; i < n_lt; ++i ) {
    TVector3 S = vec( t_sx, t_sy, t_sz, i ), E = vec( t_ex, t_ey, t_ez, i );
    double dS = ( S - lt_vtx ).Mag(), dE = ( E - lt_vtx ).Mag();
    bool start_near = dS <= dE;
    TVector3 u = ( start_near ? E : S ) - lt_vtx; // vertex -> far end
    out_dir[ i ] = u.Mag() > 0. ? u.Unit() : TVector3();
    TVector3 sd = vec( t_dx, t_dy, t_dz, i );
    bool sd_ok = std::abs( sd.X() ) <= 1. && std::abs( sd.Y() ) <= 1.
      && std::abs( sd.Z() ) <= 1.;
    kin_dir[ i ] = ( start_near && sd_ok ) ? sd : out_dir[ i ];
    len[ i ] = ( E - S ).Mag();
    int pid = std::abs( elem( t_pid, i ) );
    bool base = elem( t_sec, i ) == 0 && std::min( dS, dE ) < L_ATTACH_DIST;
    is_mu[ i ] = base && pid == MUON && std::exp( elem( t_mus, i ) ) >= L_MU_PROB_MIN;
    is_p[ i ] = base && pid == PROTON && std::exp( elem( t_prs, i ) ) >= L_P_PROB_MIN;
    pi_like[ i ] = base && pid == PI_PLUS && elem( t_E, i ) > L_PI_KE_MIN
      && std::exp( elem( t_pis, i ) ) >= L_PI_PROB_MIN;
  }
  // Merge duplicate prongs: of two collinear candidate tracks, drop the
  // shorter if it is proton-like (or if both are muon-like), and drop the
  // proton-like one of a muon-proton pair
  const double cos_merge = std::cos( L_MERGE_ANGLE * M_PI / 180. );
  for ( size_t i = 0; i < n_lt; ++i ) {
    for ( size_t j = i + 1; j < n_lt; ++j ) {
      if ( !( ( is_mu[ i ] || is_p[ i ] ) && ( is_mu[ j ] || is_p[ j ] ) ) ) continue;
      if ( out_dir[ i ].Mag() <= 0. || out_dir[ j ].Mag() <= 0. ) continue;
      if ( !( out_dir[ i ].Dot( out_dir[ j ] ) > cos_merge ) ) continue;
      size_t sh = len[ i ] < len[ j ] ? i : j;
      size_t lg = len[ i ] < len[ j ] ? j : i;
      if ( is_p[ sh ] ) drop[ sh ] = true;
      if ( is_mu[ sh ] && is_mu[ lg ] ) drop[ sh ] = true;
      if ( is_mu[ i ] && is_p[ j ] ) drop[ j ] = true;
      if ( is_p[ i ] && is_mu[ j ] ) drop[ i ] = true;
    }
  }
  const double ke_mu_min = 1e3 * ke_of_p( MUON_P_MIN, MUON_MASS );
  const double ke_p_lo = 1e3 * ke_of_p( PROTON_P_MIN, PROTON_MASS );
  const double ke_p_hi = 1e3 * ke_of_p( PROTON_P_MAX, PROTON_MASS );
  int lt_n_mu = 0, lt_n_p = 0, lt_n_p_hi = 0, lt_n_pi = 0, lt_n_shw = 0;
  int lt_mu_idx = -1;
  std::vector< std::pair< double, size_t > > lt_protons; // (p, index)
  for ( size_t i = 0; i < n_lt; ++i ) {
    double ke = elem( t_E, i );
    if ( pi_like[ i ] ) ++lt_n_pi;
    if ( drop[ i ] ) continue;
    if ( is_mu[ i ] && ke >= ke_mu_min ) { ++lt_n_mu; lt_mu_idx = static_cast< int >( i ); }
    if ( is_p[ i ] ) {
      if ( ke > ke_p_hi ) ++lt_n_p_hi;
      else if ( ke >= ke_p_lo ) {
        ++lt_n_p;
        lt_protons.emplace_back( p_of_ke( ke / 1e3, PROTON_MASS ), i );
      }
    }
  }
  for ( size_t i = 0; i < s_pid->size(); ++i ) {
    int pid = std::abs( elem( s_pid, i ) );
    if ( ( pid == ELECTRON_PDG || pid == PHOTON_PDG ) && elem( s_E, i ) > L_SHOWER_E_MAX
      && elem( s_sec, i ) == 0 ) ++lt_n_shw;
  }
  bool vote_L = lt_in_fv && lt_n_mu == 1 && lt_n_p == 2 && lt_n_p_hi == 0
    && lt_n_pi == 0 && lt_n_shw == 0 && lt_vscore >= 0. && lt_cosfrac <= 1.;

  // LANTERN kinematics: the counted muon and the two hardest counted protons
  std::stable_sort( lt_protons.begin(), lt_protons.end(),
    []( const auto& a, const auto& b ) { return a.first > b.first; } );
  bool L_ok = lt_mu_idx >= 0 && lt_protons.size() >= 2u;
  TVector3 L_mu, L_p[ 2 ];
  if ( L_ok ) {
    size_t m = static_cast< size_t >( lt_mu_idx );
    L_mu = kin_dir[ m ] * p_of_ke( elem( t_E, m ) / 1e3, MUON_MASS );
    for ( size_t r = 0; r < 2; ++r ) {
      L_p[ r ] = kin_dir[ lt_protons.at( r ).second ] * lt_protons.at( r ).first;
    }
  }

  // ------------------------------------------------------ kinematic source
  // 1: Pandora (whenever it voted), 2: LANTERN topology with the muon from
  // the aligned Pandora muon candidate (range if contained, else MCS), else
  // the aligned Wire-Cell muon (MCS), else LANTERN; the protons from
  // LANTERN, with the direction of an aligned Pandora proton candidate.
  int kin_source = 0, mu_source = -1;
  TVector3 k_mu, k_lead, k_sub;
  if ( ( vote_P || pandora_default ) && P_ok ) {
    kin_source = 1;
    k_mu = pd_.p3mu; k_lead = pd_.p3_lead_p; k_sub = pd_.p3_sub_p;
  }
  else if ( L_ok ) {
    kin_source = 2;
    k_mu = L_mu;
    mu_source = 0;
    bool wc_ok = wc_has_mu && wc_protons.size() >= 2u
      && ( wc_is_fc ? wc_p3mu.Mag() : wc_p_mcs ) > 0.;
    if ( wc_ok && wc_p_mcs > 0. && aligned( wc_p3mu, L_mu ) ) {
      k_mu = wc_p3mu.Unit() * wc_p_mcs;
      mu_source = 2;
    }
    double pd_mag = pd_.muon_contained ? pd_.p_mu_range : pd_.p_mu_mcs;
    if ( pd_.has_muon_candidate && std::isfinite( pd_mag ) && pd_mag > 0.05
      && pd_mag < 5. && aligned( pd_.p3mu, L_mu ) )
    {
      k_mu = pd_.p3mu.Unit() * pd_mag;
      mu_source = 1;
    }
    // Pandora proton directions are used only when Pandora found the muon
    // candidate as well (otherwise its "proton" candidates are all tracks)
    TVector3 prot[ 2 ] = { L_p[ 0 ], L_p[ 1 ] };
    if ( pd_.has_muon_candidate ) {
      for ( size_t r = 0; r < 2; ++r ) {
        if ( pd_.has_lead_p && aligned( pd_.p3_lead_p, prot[ r ] ) ) {
          prot[ r ] = pd_.p3_lead_p.Unit() * prot[ r ].Mag();
        }
        if ( pd_.has_sub_p && aligned( pd_.p3_sub_p, prot[ r ] ) ) {
          prot[ r ] = pd_.p3_sub_p.Unit() * prot[ r ].Mag();
        }
      }
    }
    bool swap = prot[ 1 ].Mag() > prot[ 0 ].Mag();
    k_lead = swap ? prot[ 1 ] : prot[ 0 ];
    k_sub = swap ? prot[ 0 ] : prot[ 1 ];
  }
  bool in_phase_space = kin_source > 0 && k_mu.Mag() >= MUON_P_MIN
    && k_sub.Mag() >= PROTON_P_MIN && k_lead.Mag() <= PROTON_P_MAX;

  // ------------------------------------------------------------ decision
  int n_votes = int( vote_P ) + int( vote_W ) + int( vote_L );
  bool two_of_three = n_votes >= 2;
  // Events found only by Wire-Cell and LANTERN need a Pandora
  // (SCE-corrected) vertex in the fiducial volume
  bool pandora_vertex_veto = two_of_three && !vote_P
    && !( pd_.has_vertex && pd_.reco_vertex_in_FV );
  bool collinear_rescue = pandora_default && pd_.cos_gamma_pp > RESCUE_COS_GAMMA_MIN;
  bool selected = ( two_of_three && in_phase_space && !pandora_vertex_veto )
    || collinear_rescue;

  // Observables from the assigned kinematics (Pandora's are already stored)
  if ( kin_source == 2 ) {
    fill_observables( out, true, k_mu, k_lead, k_sub );
    out[ "reco_p3_p_vec" ] = std::vector< TVector3 >{ k_lead, k_sub };
  }

  // Sidebands for background constraints, disjoint from the selection and
  // from each other: the Pandora topologies of CC1mu2p0piv2 (CC0pi >= 3p,
  // CC Npi, and CC0pi 1p). The 1p sideband also vetoes events in which
  // Wire-Cell or LANTERN finds a second proton in or above the proton
  // window (as counted in their votes), i.e., signal with a proton that
  // Pandora missed. Sideband codes: 0 none, 1 3p, 2 Npi, 3 1p.
  bool sb_veto_1p = wc_pmag( 1 ) >= PROTON_P_MIN || lt_n_p + lt_n_p_hi >= 2;
  bool sb_3p = pd_.sb_3p && !selected;
  bool sb_npi = pd_.sb_npi && !selected;
  bool sb_1p = pd_.sb_1p && !sb_veto_1p && !selected;
  int sideband = sb_3p ? 1 : sb_npi ? 2 : sb_1p ? 3 : 0;
  in_sideband_ = sideband > 0;
  out[ "sideband_3p" ] = sb_3p;
  out[ "sideband_npi" ] = sb_npi;
  out[ "sideband_1p" ] = sb_1p;
  out[ "sideband" ] = sideband;
  out[ "sb_veto_1p" ] = sb_veto_1p;

  out[ "vote_pandora" ] = vote_P;
  out[ "vote_wirecell" ] = vote_W;
  out[ "vote_lantern" ] = vote_L;
  out[ "n_votes" ] = n_votes;
  out[ "pandora_default" ] = pandora_default;
  out[ "kin_source" ] = kin_source;
  out[ "mu_source" ] = mu_source;
  out[ "in_reco_phase_space" ] = in_phase_space;
  out[ "pandora_vertex_veto" ] = pandora_vertex_veto;
  out[ "collinear_rescue" ] = collinear_rescue;
  out[ "wc_numu_score" ] = static_cast< double >( wc_score );
  out[ "wc_p1" ] = wc_pmag( 0 );
  out[ "wc_p2" ] = wc_pmag( 1 );
  out[ "wc_p3" ] = wc_pmag( 2 );
  out[ "wc_n_pi" ] = wc_n_pi;
  out[ "wc_has_mu" ] = wc_has_mu;
  out[ "wc_vertex_in_FV" ] = wc_in_fv;
  out[ "lt_vertex_in_FV" ] = lt_in_fv;
  out[ "lt_n_mu" ] = lt_n_mu;
  out[ "lt_n_p" ] = lt_n_p;
  out[ "lt_n_p_hi" ] = lt_n_p_hi;
  out[ "lt_n_pi" ] = lt_n_pi;
  out[ "lt_n_shower" ] = lt_n_shw;

  return selected;
}
