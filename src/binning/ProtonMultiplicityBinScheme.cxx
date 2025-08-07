// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/ProtonMultiplicityBinScheme.hh"

ProtonMultiplicityBinScheme::ProtonMultiplicityBinScheme() : BinSchemeBase( "ProtonMultiplicityBinScheme" ) {}

void ProtonMultiplicityBinScheme::DefineBlocks() {

  // Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  runs_to_use_ = { 1 };

  // Prefix for the output bin and slice configuration text files
  out_config_prefix_ = "proton_multiplicity_";

  // Selection to use with this binning scheme
  selection_name_ = "CC1muNp0piEnriched";

  // TDirectory file name to use when producing the univmake output histograms
  out_tdir_name_ = "proton_multiplicity";

  // Signal definition and selection branch names
  std::string signal = selection_name_ + "_" + "MC_Signal";
  std::string selection = selection_name_ + "_" + "Selected";
  //std::string signal = "MC_Signal";
  //std::string selection = "Selected";

  /* -------------------------------------------------------------------------- */
  /*                          Analysis bin definitions                          */
  /* -------------------------------------------------------------------------- */
  
  // Define the blocks of bins in both true and reco space

  // Block 0: leading proton momentum (1D block)
  std::vector< double > BLOCK_0_BIN_EDGES = { 0.250, 0.325, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 1.000 };

  std::string block0_true_branchexpr = "CC1muNp0piEnriched_true_p3_lead_p.Mag();GeV/c";
  std::string block0_reco_branchexpr = "CC1muNp0piEnriched_reco_p3_lead_p.Mag();GeV/c";
  std::string block0_title = "p_{l_{1}};GeV/c";
  std::string block0_textitle = "p_{l_{1}};GeV/c";

  Block1D* b0t = new Block1D(block0_true_branchexpr, block0_title, block0_textitle, BLOCK_0_BIN_EDGES, signal, kSignalTrueBin);
  Block1D* b0r = new Block1D(block0_reco_branchexpr, block0_title, block0_textitle, BLOCK_0_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b0t, b0r );

  // Block 1: proton multiplicity in slices of leading proton momentum (2D block)
  std::map< double, std::vector<double> > BLOCK_1_BIN_EDGES = {

    { 0.250, { 1, 2, 3, 10} },
    { 0.325, { 1, 2, 3, 4, 10} },
    { 0.400, { 1, 2, 3, 4, 10} },
    { 0.450, { 1, 2, 3, 4, 10} },
    { 0.500, { 1, 2, 3, 4, 10} },
    { 0.550, { 1, 2, 3, 4, 10} },
    { 0.600, { 1, 2, 3, 4, 10} },
    { 0.650, { 1, 2, 3, 4, 10} },
    { 0.700, { 1, 2, 3, 4, 10} },
    { 0.750, { 1, 2, 3, 4, 10} },
    { 0.800, { 1, 2, 3, 10} },
    { 0.850, { 1, 2, 3, 10} },
    { 0.900, { 1, 2, 3, 10} },
    { 1.000, { 1, 2, 3, 10} },
  
  }; // above is the bins for the analysis up to the July 2022 uBooNE CM

  std::string block1_true_branchexpr = "CC1muNp0piEnriched_true_p3_lead_p.Mag();GeV/c;Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);";
  std::string block1_reco_branchexpr = "CC1muNp0piEnriched_reco_p3_lead_p.Mag();GeV/c;Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);";
  std::string block1_title = "p_{l_{1}};GeV/c;p_{mult};";
  std::string block1_textitle = "p_{l_{1}};GeV/c;p_{mult};";

  Block2D* b1t = new Block2D(block1_true_branchexpr, block1_title, block1_textitle, BLOCK_1_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b1r = new Block2D(block1_reco_branchexpr, block1_title, block1_textitle, BLOCK_1_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b1t, b1r );

  // Block 2: leading proton momentum in slices of proton multiplicity (2D block)
  std::map< double, std::vector<double> > BLOCK_2_BIN_EDGES = {

    { 1,  { 0.250, 0.350, 0.450, 0.550, 0.650, 0.750, 0.850, 1.000 } },
    { 2,  { 0.250, 0.350, 0.450, 0.550, 0.650, 0.750, 0.850, 1.000 } },
    { 3,  { 0.250, 0.350, 0.450, 0.550, 0.650, 0.750, 1.000 } },
    { 4,  { 0.250, 0.450, 0.550, 0.650, 0.750, 1.000 } },
    { 10, { 0.250, 0.450, 0.550, 0.650, 0.750, 1.000 } },

  };

  std::string block2_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_p3_p_vec[0].Mag();GeV/c";
  std::string block2_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p3_p_vec[0].Mag();GeV/c";
  std::string block2_title = "p_{mult};;p_{l_{1}};GeV/c";
  std::string block2_textitle = "p_{mult};;p_{l_{1}};GeV/c";

  Block2D* b2t = new Block2D(block2_true_branchexpr, block2_title, block2_textitle, BLOCK_2_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b2r = new Block2D(block2_reco_branchexpr, block2_title, block2_textitle, BLOCK_2_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b2t, b2r );

  // Block 3: subleading proton momentum in slices of proton multiplicity (2D block)
  std::map< double, std::vector<double> > BLOCK_3_BIN_EDGES = {

    { 2,  { 0.250, 0.350, 0.450, 0.550, 1.000 } },
    { 3,  { 0.250, 0.350, 0.450, 0.550, 1.000 } },
    { 4,  { 0.250, 0.350, 0.450, 0.550, 1.000 } },
    { 10, { 0.250, 0.350, 0.450, 0.550, 1.000 } },

  };

  std::string block3_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_p3_p_vec[1].Mag();GeV/c";
  std::string block3_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p3_p_vec[1].Mag();GeV/c";
  std::string block3_title = "p_{mult};;p_{l_{2}};GeV/c";
  std::string block3_textitle = "p_{mult};;p_{l_{2}};GeV/c";

  Block2D* b3t = new Block2D(block3_true_branchexpr, block3_title, block3_textitle, BLOCK_3_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b3r = new Block2D(block3_reco_branchexpr, block3_title, block3_textitle, BLOCK_3_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b3t, b3r );

  // Block 4: third-leading proton momentum in slices of proton multiplicity (2D block)
  std::map< double, std::vector<double> > BLOCK_4_BIN_EDGES = {

    { 3,  { 0.250, 0.350, 0.450, 1.000 } },
    { 4,  { 0.250, 0.350, 0.450, 1.000 } },
    { 10, { 0.250, 0.350, 0.450, 1.000 } },

  };

  std::string block4_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_p3_p_vec[2].Mag();GeV/c";
  std::string block4_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p3_p_vec[2].Mag();GeV/c";
  std::string block4_title = "p_{mult};;p_{l_{3}};GeV/c";
  std::string block4_textitle = "p_{mult};;p_{l_{3}};GeV/c";

  Block2D* b4t = new Block2D(block4_true_branchexpr, block4_title, block4_textitle, BLOCK_4_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b4r = new Block2D(block4_reco_branchexpr, block4_title, block4_textitle, BLOCK_4_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b4t, b4r );

  // Block 5: fourth-leading proton momentum in slices of proton multiplicity (2D block)
  std::map< double, std::vector<double> > BLOCK_5_BIN_EDGES = {

    { 4,  { 0.250, 0.350, 0.450, 1.000 } },
    { 10, { 0.250, 0.350, 0.450, 1.000 } },

  };

  std::string block5_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_p3_p_vec[3].Mag();GeV/c";
  std::string block5_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p3_p_vec[3].Mag();GeV/c";
  std::string block5_title = "p_{mult};;p_{l_{4}};GeV/c";
  std::string block5_textitle = "p_{mult};;p_{l_{4}};GeV/c";

  Block2D* b5t = new Block2D(block5_true_branchexpr, block5_title, block5_textitle, BLOCK_5_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b5r = new Block2D(block5_reco_branchexpr, block5_title, block5_textitle, BLOCK_5_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b5t, b5r );

  //Block 6: available energy (1D block)
  std::vector< double > BLOCK_6_BIN_EDGES = { 0.000, 0.050, 0.100, 0.200, 0.300, 0.400, 0.600 };

  std::string block6_true_branchexpr = "CC1muNp0piEnriched_true_e_had_vis;GeV";
  std::string block6_reco_branchexpr = "CC1muNp0piEnriched_reco_e_had_vis;GeV";
  std::string block6_title = "E_{avail};GeV";
  std::string block6_textitle = "E_{avail};GeV";

  Block1D* b6t = new Block1D(block6_true_branchexpr, block6_title, block6_textitle, BLOCK_6_BIN_EDGES, signal, kSignalTrueBin);
  Block1D* b6r = new Block1D(block6_reco_branchexpr, block6_title, block6_textitle, BLOCK_6_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b6t, b6r );

  // Block 7: proton multiplicity in slices of available energy (2D block)
  std::map< double, std::vector<double> > BLOCK_7_BIN_EDGES = {

    //{ 0.000, { 1, 10} },
    { 0.050, { 1, 2, 10} },
    { 0.100, { 1, 2, 3, 10} },
    { 0.200, { 1, 2, 3, 4, 10} },
    { 0.300, { 1, 2, 3, 4, 10} },
    { 0.400, { 1, 2, 3, 4, 10} },
    { 0.600, { 1, 2, 3, 4, 10} },
  
  };

  std::string block7_true_branchexpr = "CC1muNp0piEnriched_true_e_had_vis;GeV;Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);";
  std::string block7_reco_branchexpr = "CC1muNp0piEnriched_reco_e_had_vis;GeV;Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);";
  std::string block7_title = "E_{avail};GeV;p_{mult};";
  std::string block7_textitle = "E_{avail};GeV;p_{mult};";

  Block2D* b7t = new Block2D(block7_true_branchexpr, block7_title, block7_textitle, BLOCK_7_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b7r = new Block2D(block7_reco_branchexpr, block7_title, block7_textitle, BLOCK_7_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b7t, b7r );

  // Block 8: opening angle between muon and leading proton (1D block)

  std::vector< double > BLOCK_8_BIN_EDGES = { -1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00 };

  std::string block8_true_branchexpr = "CC1muNp0piEnriched_true_cos_theta_mu_p;";
  std::string block8_reco_branchexpr = "CC1muNp0piEnriched_reco_cos_theta_mu_p;";
  std::string block8_title = "cos #theta_{l_{1}#mu};";
  std::string block8_textitle = "\\cos \\theta_{l_{1}\\mu};";

  Block1D* b8t = new Block1D(block8_true_branchexpr, block8_title, block8_textitle, BLOCK_8_BIN_EDGES, signal, kSignalTrueBin);
  Block1D* b8r = new Block1D(block8_reco_branchexpr, block8_title, block8_textitle, BLOCK_8_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b8t, b8r );

  // Block 9: proton multiplicity in slices of opening angle between muon and leading proton (2D block)
  std::map< double, std::vector<double> > BLOCK_9_BIN_EDGES = {

    { -1.00, { 1, 2, 3, 4, 10} },
    { -0.75, { 1, 2, 3, 4, 10} },
    { -0.50, { 1, 2, 3, 4, 10} },
    { -0.25, { 1, 2, 3, 4, 10} },
    {  0.00, { 1, 2, 3, 4, 10} },
    {  0.25, { 1, 2, 3, 4, 10} },
    {  0.50, { 1, 2, 3, 4, 10} },
    {  0.75, { 1, 2, 3, 4, 10} },
    {  1.00, { 1, 2, 3, 4, 10} },
  
  };

  std::string block9_true_branchexpr = "CC1muNp0piEnriched_true_cos_theta_mu_p;;Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);";
  std::string block9_reco_branchexpr = "CC1muNp0piEnriched_reco_cos_theta_mu_p;;Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);";
  std::string block9_title = "cos #theta_{l_{1}#mu};;p_{mult};";
  std::string block9_textitle = "\\cos \\theta_{l_{1}\\mu};;p_{mult};";

  Block2D* b9t = new Block2D(block9_true_branchexpr, block9_title, block9_textitle, BLOCK_9_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b9r = new Block2D(block9_reco_branchexpr, block9_title, block9_textitle, BLOCK_9_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b9t, b9r );

  // Block 10: leading proton opening angle in slices of proton multiplicity (2D block)
  std::map< double, std::vector<double> > BLOCK_10_BIN_EDGES = {

    { 1,  { -1.00, -0.66, -0.33, 0.00, 0.33, 0.66, 1.00 } },
    { 2,  { -1.00, -0.66, -0.33, 0.00, 0.33, 0.66, 1.00 } },
    { 3,  { -1.00, -0.66, -0.33, 0.00, 0.33, 0.66, 1.00 } },
    { 4,  { -1.00, -0.66, -0.33, 0.00, 0.33, 0.66, 1.00 } },
    { 10, { -1.00, -0.66, -0.33, 0.00, 0.33, 0.66, 1.00 } },

  };

  std::string block10_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_cos_theta_mu_p;";
  std::string block10_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_p;";
  std::string block10_title = "p_{mult};;cos #theta_{l_{1}#mu};";
  std::string block10_textitle = "p_{mult};;\\cos \\theta_{l_{1}\\mu};";

  Block2D* b10t = new Block2D(block10_true_branchexpr, block10_title, block10_textitle, BLOCK_10_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b10r = new Block2D(block10_reco_branchexpr, block10_title, block10_textitle, BLOCK_10_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b10t, b10r );

  // Block 11: subleading proton opening angle in slices of proton multiplicity (2D block)
  std::map< double, std::vector<double> > BLOCK_11_BIN_EDGES = {

    { 2,  { -1.00, -0.66, -0.33, 0.00, 0.33, 0.66, 1.00 } },
    { 3,  { -1.00, -0.66, -0.33, 0.00, 0.33, 0.66, 1.00 } },
    { 4,  { -1.00, -0.66, -0.33, 0.00, 0.33, 0.66, 1.00 } },
    { 10, { -1.00, -0.66, -0.33, 0.00, 0.33, 0.66, 1.00 } },

  };

  std::string block11_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_cos_theta_mu_p2;";
  std::string block11_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_p2;";
  std::string block11_title = "p_{mult};;cos #theta_{l_{2}#mu};";
  std::string block11_textitle = "p_{mult};;\\cos \\theta_{l_{2}\\mu};";

  Block2D* b11t = new Block2D(block11_true_branchexpr, block11_title, block11_textitle, BLOCK_11_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b11r = new Block2D(block11_reco_branchexpr, block11_title, block11_textitle, BLOCK_11_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b11t, b11r );

  // Block 12: third-leading proton opening angle in slices of proton multiplicity (2D block)
  std::map< double, std::vector<double> > BLOCK_12_BIN_EDGES = {

    { 3,  { -1.00, -0.5, 0.00, 0.50, 1.00 } },
    { 4,  { -1.00, -0.5, 0.00, 0.50, 1.00 } },
    { 10, { -1.00, -0.5, 0.00, 0.50, 1.00 } },

  };

  std::string block12_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_cos_theta_mu_p3;";
  std::string block12_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_p3;";
  std::string block12_title = "p_{mult};;cos #theta_{l_{3}#mu};";
  std::string block12_textitle = "p_{mult};;\\cos \\theta_{l_{3}\\mu};";

  Block2D* b12t = new Block2D(block12_true_branchexpr, block12_title, block12_textitle, BLOCK_12_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b12r = new Block2D(block12_reco_branchexpr, block12_title, block12_textitle, BLOCK_12_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b12t, b12r );

  // Block 13: fourth-leading proton opening angle in slices of proton multiplicity (2D block)
  std::map< double, std::vector<double> > BLOCK_13_BIN_EDGES = {

    { 4,  { -1.00, -0.5, 0.00, 0.50, 1.00 } },
    { 10, { -1.00, -0.5, 0.00, 0.50, 1.00 } },

  };

  std::string block13_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_cos_theta_mu_p4;";
  std::string block13_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_p4;";
  std::string block13_title = "p_{mult};;cos #theta_{l_{4}#mu};";
  std::string block13_textitle = "p_{mult};;\\cos \\theta_{l_{4}\\mu};";

  Block2D* b13t = new Block2D(block13_true_branchexpr, block13_title, block13_textitle, BLOCK_13_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b13r = new Block2D(block13_reco_branchexpr, block13_title, block13_textitle, BLOCK_13_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b13t, b13r );

  /* -------------------------------------------------------------------------- */
  /*                           Side bands definitions                           */
  /* -------------------------------------------------------------------------- */

  // Don't add side bands for now...

  /* std::vector< double > PROTON_1D_BIN_EDGES = { 0.250, 0.325, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900 };

  std::string side_branchexpr = "CC1muNp0piEnriched_reco_p3_p_vec.Mag();GeV/c";
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

  vect_sideband.emplace_back( b4s ); */

  /* -------------------------------------------------------------------------- */
  /*                            Background categories                           */
  /* -------------------------------------------------------------------------- */

  // Add relevant background categories
  CATEGORY = selection_name_ + "_" + "EventCategory";
  //CATEGORY = "EventCategory";
  background_index = {1, 2, 3, 4, 17, 18, 19, 20, 21, 22};

}
