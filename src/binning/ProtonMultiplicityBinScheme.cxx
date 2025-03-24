// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/ProtonMultiplicityBinScheme.hh"

ProtonMultiplicityBinScheme::ProtonMultiplicityBinScheme() : BinSchemeBase( "ProtonMultiplicityBinScheme" ) {}

void ProtonMultiplicityBinScheme::DefineBlocks() {

  /////// Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  runs_to_use_ = { 1 };

  // Prefix for the output bin and slice configuration text files
  out_config_prefix_ = "proton_multiplicity_";

  // Selection to use with this binning scheme
  selection_name_ = "CC1muNp0pi";

  // TDirectory file name to use when producing the univmake output histograms
  out_tdir_name_ = "proton_multiplicity";

  /////// Define the blocks of bins in both true and reco space

  // 2D block: proton multiplicity in slices of leading proton momentum
  std::map< double, std::vector<double> > PROTON_2D_BIN_EDGES = {

    { 0.250, { 1, 2, 3, 6} },
    { 0.325, { 1, 2, 3, 4, 6} },
    { 0.4,   { 1, 2, 3, 4, 6} },
    { 0.45,  { 1, 2, 3, 4, 6} },
    { 0.5,   { 1, 2, 3, 4, 6} },
    { 0.550, { 1, 2, 3, 4, 6} },
    { 0.6,   { 1, 2, 3, 4, 6} },
    { 0.65,  { 1, 2, 3, 4, 6} },
    { 0.7,   { 1, 2, 3, 4, 6} },
    { 0.75,  { 1, 2, 3, 4, 6} },
    { 0.8,   { 1, 2, 3, 6} },
    { 0.85,  { 1, 2, 3, 6} },
    { 0.9,   { 1, 2, 3, 6} },
  
  }; //above is the bins for the analysis up to the July 2022 uBooNE CM

  std::string true_branchexpr = "mc_p3_lead_p.Mag(); GeV/c; Sum$(mc_p3_p_vec.Mag() > 0.25); ";
  std::string reco_branchexpr = "p3_lead_p.Mag(); GeV/c; Sum$(p3_p_vec.Mag() > 0.25); ";
  std::string title = "p_{p}; [GeV/c]; p_{mult}; ";
  std::string textitle = "p_{p}; [GeV/c]; p_{mult}; ";
  std::string signal = "mc_is_signal"
  std::string selection = "sel_CCNp0pi";

  Block2D* b1t = new Block2D(true_branchexpr, title, textitle, PROTON_2D_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b1r = new Block2D(reco_branchexpr, title, textitle, PROTON_2D_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b1t, b1r );

  // Dirt sideband
  const std::string DIRT_SIDEBAND_SELECTION =
    "!sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
    " && sel_has_muon_candidate && sel_topo_cut_passed"
    " && sel_no_reco_showers && sel_muon_passed_mom_cuts"
    " && sel_muon_contained && sel_muon_quality_ok"
    " && sel_has_p_candidate && sel_passed_proton_pid_cut"
    " && sel_protons_contained && sel_lead_p_passed_mom_cuts";

  // NC sideband
  const std::string NC_SIDEBAND_SELECTION =
    "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
    " && !sel_has_muon_candidate && sel_topo_cut_passed"
    " && sel_no_reco_showers && sel_has_p_candidate"
    " && sel_passed_proton_pid_cut && sel_protons_contained"
    " && sel_lead_p_passed_mom_cuts";

  // CCNpi sideband
  const std::string CCNPI_SIDEBAND_SELECTION =
    "sel_reco_vertex_in_FV && sel_pfp_starts_in_PCV"
    " && sel_has_muon_candidate && sel_topo_cut_passed"
    " && sel_no_reco_showers && sel_muon_passed_mom_cuts"
    " && sel_muon_contained && sel_muon_quality_ok"
    " && sel_has_p_candidate && sel_protons_contained"
    " && sel_lead_p_passed_mom_cuts"
    " && trk_llr_pid_score_v[ lead_p_candidate_idx ] > 0.2 ";

}
