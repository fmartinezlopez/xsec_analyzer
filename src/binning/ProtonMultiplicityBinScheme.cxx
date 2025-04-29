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
    { 0.400, { 1, 2, 3, 4, 6} },
    { 0.450, { 1, 2, 3, 4, 6} },
    { 0.500, { 1, 2, 3, 4, 6} },
    { 0.550, { 1, 2, 3, 4, 6} },
    { 0.600, { 1, 2, 3, 4, 6} },
    { 0.650, { 1, 2, 3, 4, 6} },
    { 0.700, { 1, 2, 3, 4, 6} },
    { 0.750, { 1, 2, 3, 4, 6} },
    { 0.800, { 1, 2, 3, 6} },
    { 0.850, { 1, 2, 3, 6} },
    { 0.900, { 1, 2, 3, 6} },
  
  }; // above is the bins for the analysis up to the July 2022 uBooNE CM

  //std::string true_branchexpr = "mc_p3_lead_p.Mag(); GeV/c; Sum$(mc_p3_p_vec.Mag() > 0.25); ";
  std::string true_branchexpr = "true_p3_lead_p.Mag();GeV/c;Sum$(true_p3_p_vec.Mag() > 0.25);";
  //std::string reco_branchexpr = "p3_lead_p.Mag(); GeV/c; Sum$(p3_p_vec.Mag() > 0.25); ";
  std::string reco_branchexpr = "reco_p3_lead_p.Mag();GeV/c;Sum$(reco_p3_p_vec.Mag() > 0.25);";
  std::string title = "p_{p};GeV/c;p_{mult};";
  std::string textitle = "p_{p};GeV/c;p_{mult};";
  //std::string signal = "mc_is_signal";
  std::string signal = "MC_Signal";
  //std::string selection = "sel_CCNp0pi";
  std::string selection = "Selected";

  Block2D* b1t = new Block2D(true_branchexpr, title, textitle, PROTON_2D_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b1r = new Block2D(reco_branchexpr, title, textitle, PROTON_2D_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b1t, b1r );

  /* ------------------------------- Side bands ------------------------------- */

  std::vector< double > PROTON_1D_BIN_EDGES = { 0.250, 0.325, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900 };

  std::string side_branchexpr = "reco_p3_lead_p.Mag();GeV/c";
  std::string side_title = "p_{p};GeV/c";
  std::string side_textitle = "p_{p};GeV/c";

  // Dirt sideband
  const std::string DIRT_SIDEBAND_SELECTION =
    "!reco_vertex_in_FV && pfp_starts_in_PCV"
    " && has_muon_candidate && topo_cut_passed"
    " && no_reco_showers && muon_passed_mom_cuts"
    " && muon_contained && muon_quality_ok"
    " && has_p_candidate && passed_proton_pid_cut"
    " && protons_contained && lead_p_passed_mom_cuts";

  Block1D* b2s = new Block1D(side_branchexpr, side_title, side_textitle, PROTON_1D_BIN_EDGES, DIRT_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b2s );

  // NC sideband
  const std::string NC_SIDEBAND_SELECTION =
    "reco_vertex_in_FV && pfp_starts_in_PCV"
    " && !has_muon_candidate && topo_cut_passed"
    " && no_reco_showers && has_p_candidate"
    " && passed_proton_pid_cut && protons_contained"
    " && lead_p_passed_mom_cuts";

  Block1D* b3s = new Block1D(side_branchexpr, side_title, side_textitle, PROTON_1D_BIN_EDGES, NC_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b3s );

  // CCNpi sideband
  const std::string CCNPI_SIDEBAND_SELECTION =
    "reco_vertex_in_FV && pfp_starts_in_PCV"
    " && has_muon_candidate && topo_cut_passed"
    " && no_reco_showers && muon_passed_mom_cuts"
    " && muon_contained && muon_quality_ok"
    " && has_p_candidate && protons_contained"
    " && lead_p_passed_mom_cuts"
    " && trk_llr_pid_score_v[ lead_p_candidate_idx ] > 0.2 ";

  Block1D* b4s = new Block1D(side_branchexpr, side_title, side_textitle, PROTON_1D_BIN_EDGES, CCNPI_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b4s );

  /* -------------------------- Background categories ------------------------- */

  // Add relevant background categories
  //CATEGORY = "category";
  CATEGORY = "EventCategory";
  //background_index = {5, 6, 7, 8, 9, 10, 11};
  background_index = {1, 2, 3, 4, 17, 18, 19, 20, 21, 22};

}
