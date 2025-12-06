// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/ProtonMultiplicityBinScheme.hh"

ProtonMultiplicityBinScheme::ProtonMultiplicityBinScheme() : BinSchemeBase( "ProtonMultiplicityBinScheme" ) {}

void ProtonMultiplicityBinScheme::DefineBlocks() {

  // Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  //runs_to_use_ = { 1, 2, 3, 40, 41, 42, 43, 5 };
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

  /* // Block 1: proton multiplicity in slices of leading proton momentum (2D block)

  std::map< double, std::vector<double> > BLOCK_1_BIN_EDGES = {

    { 0.250, { 1, 2, 10} },
    { 0.375, { 1, 2, 3, 10} },
    { 0.475, { 1, 2, 3, 4, 10} },
    { 0.575, { 1, 2, 3, 4, 10} },
    { 0.650, { 1, 2, 3, 4, 10} },
    { 0.725, { 1, 2, 3, 4, 10} },
    { 0.800, { 1, 2, 3, 10} },
    { 0.850, { 1, 2, 3, 10} },
    { 1.000, { 1, 2, 3, 10} },
  
  }; // version in technote (uncertainty opt.)

  std::string block1_true_branchexpr = "CC1muNp0piEnriched_true_p_li_vec[0];GeV/c;Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);";
  std::string block1_reco_branchexpr = "CC1muNp0piEnriched_reco_p_li_vec[0];GeV/c;Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);";
  std::string block1_title = "p_{l_{1}};GeV/c;p_{mult};";
  std::string block1_textitle = "p_{l_{1}};GeV/c;p_{mult};";

  Block2D* b1t = new Block2D(block1_true_branchexpr, block1_title, block1_textitle, BLOCK_1_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b1r = new Block2D(block1_reco_branchexpr, block1_title, block1_textitle, BLOCK_1_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b1t, b1r );

  // Block 2: proton multiplicity in slices of opening angle between muon and leading proton (2D block)

  std::map< double, std::vector<double> > BLOCK_2_BIN_EDGES = {

    { -1.00, { 1, 2, 3, 10} },
    { -0.85, { 1, 2, 3, 10} },
    { -0.65, { 1, 2, 3, 4, 10} },
    {  0.00, { 1, 2, 3, 4, 10} },
    {  0.65, { 1, 2, 3, 10} },
    {  0.85, { 1, 2, 3, 10} },
    {  1.00, { 1, 2, 3, 10} },
  
  }; // version in technote

  std::string block2_true_branchexpr = "CC1muNp0piEnriched_true_cos_theta_mu_li_vec[0];;Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);";
  std::string block2_reco_branchexpr = "CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[0];;Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);";
  std::string block2_title = "cos #theta_{l_{1}#mu};;p_{mult};";
  std::string block2_textitle = "\\cos \\theta_{l_{1}\\mu};;p_{mult};";

  Block2D* b2t = new Block2D(block2_true_branchexpr, block2_title, block2_textitle, BLOCK_2_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b2r = new Block2D(block2_reco_branchexpr, block2_title, block2_textitle, BLOCK_2_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b2t, b2r );

  // Block 3: proton multiplicity in slices of available energy (2D block)

  std::map< double, std::vector<double> > BLOCK_3_BIN_EDGES = {

    { 0.050, { 1, 2, 10} },
    { 0.125, { 1, 2, 3, 10} },
    { 0.285, { 1, 2, 3, 10} },
    { 0.400, { 1, 2, 3, 4, 10} },
    { 1.000, { 1, 2, 3, 4, 10} },
  
  }; // version in technote (uncertainty opt.)

  std::string block3_true_branchexpr = "CC1muNp0piEnriched_true_e_had_vis;GeV;Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);";
  std::string block3_reco_branchexpr = "CC1muNp0piEnriched_reco_e_had_vis;GeV;Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);";
  std::string block3_title = "E_{avail};GeV;p_{mult};";
  std::string block3_textitle = "E_{avail};GeV;p_{mult};";

  Block2D* b3t = new Block2D(block3_true_branchexpr, block3_title, block3_textitle, BLOCK_3_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b3r = new Block2D(block3_reco_branchexpr, block3_title, block3_textitle, BLOCK_3_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b3t, b3r );

  // Block 4: leading proton momentum in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_4_BIN_EDGES = {

    { 1,  { 0.250, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.900, 1.000 } },
    { 2,  { 0.250, 0.400, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 1.000 } },
    { 3,  { 0.250, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 1.000 } },
    { 4,  { 0.250, 0.600, 0.700, 0.800, 1.000 } },
    { 10, { 0.250, 0.600, 0.700, 0.800, 1.000 } },

  }; // version in technote (uncertainty opt.)

  std::string block4_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_p_li_vec[0];GeV/c";
  std::string block4_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p_li_vec[0];GeV/c";
  std::string block4_title = "p_{mult};;p_{l_{1}};GeV/c";
  std::string block4_textitle = "p_{mult};;p_{l_{1}};GeV/c";

  Block2D* b4t = new Block2D(block4_true_branchexpr, block4_title, block4_textitle, BLOCK_4_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b4r = new Block2D(block4_reco_branchexpr, block4_title, block4_textitle, BLOCK_4_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b4t, b4r );

  // Block 5: subleading proton momentum in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_5_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 2,  { 0.250, 0.350, 0.400, 0.450, 0.500, 0.550, 1.000 } },
    { 3,  { 0.250, 0.450, 0.500, 0.550, 0.600, 0.650, 1.000 } },
    { 4,  { 0.250, 0.550, 0.600, 1.000 } },
    { 10, { 0.250, 0.550, 0.600, 1.000 } },
  
  }; // v3

  std::string block5_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_p_li_vec[1];GeV/c";
  std::string block5_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p_li_vec[1];GeV/c";
  std::string block5_title = "p_{mult};;p_{l_{2}};GeV/c";
  std::string block5_textitle = "p_{mult};;p_{l_{2}};GeV/c";

  Block2D* b5t = new Block2D(block5_true_branchexpr, block5_title, block5_textitle, BLOCK_5_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b5r = new Block2D(block5_reco_branchexpr, block5_title, block5_textitle, BLOCK_5_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b5t, b5r );

  // Block 6: third-leading proton momentum in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_6_BIN_EDGES = {

    { 1,   { -10000.0, 10000.0 } },
    { 3,   { 0.250, 0.350, 0.400, 0.450, 0.500, 1.000 } },
    { 4,   { 0.250, 0.450, 0.500, 1.000 } },
    { 10,  { 0.250, 0.450, 0.500, 1.000 } },

  }; // v3

  std::string block6_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_p_li_vec[2];GeV/c";
  std::string block6_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p_li_vec[2];GeV/c";
  std::string block6_title = "p_{mult};;p_{l_{3}};GeV/c";
  std::string block6_textitle = "p_{mult};;p_{l_{3}};GeV/c";

  Block2D* b6t = new Block2D(block6_true_branchexpr, block6_title, block6_textitle, BLOCK_6_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b6r = new Block2D(block6_reco_branchexpr, block6_title, block6_textitle, BLOCK_6_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b6t, b6r );

  // Block 7: fourth-leading proton momentum in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_7_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 4,  { 0.250, 0.400, 1.000 } },
    { 10, { 0.250, 0.400, 1.000 } },

  }; // v3
  
  std::string block7_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_p_li_vec[3];GeV/c";
  std::string block7_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p_li_vec[3];GeV/c";
  std::string block7_title = "p_{mult};;p_{l_{4}};GeV/c";
  std::string block7_textitle = "p_{mult};;p_{l_{4}};GeV/c";

  Block2D* b7t = new Block2D(block7_true_branchexpr, block7_title, block7_textitle, BLOCK_7_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b7r = new Block2D(block7_reco_branchexpr, block7_title, block7_textitle, BLOCK_7_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b7t, b7r );

  // Block 8: leading proton opening angle in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_8_BIN_EDGES = {

    { 1,  { -1.000, -0.800, -0.600, -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000 } },
    { 2,  { -1.000, -0.800, -0.600, -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000 } },
    { 3,  { -1.000, -0.750, -0.500, -0.250, 0.000, 0.250, 0.500, 0.750, 1.000 } },
    { 4,  { -1.000, -0.500, 0.000, 0.500, 1.000 } },
    { 10, { -1.000, -0.500, 0.000, 0.500, 1.000 } },

  }; // version in technote

  std::string block8_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_cos_theta_mu_li_vec[0];";
  std::string block8_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[0];";
  std::string block8_title = "p_{mult};;cos #theta_{l_{1}#mu};";
  std::string block8_textitle = "p_{mult};;\\cos \\theta_{l_{1}\\mu};";

  Block2D* b8t = new Block2D(block8_true_branchexpr, block8_title, block8_textitle, BLOCK_8_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b8r = new Block2D(block8_reco_branchexpr, block8_title, block8_textitle, BLOCK_8_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b8t, b8r );

  // Block 9: subleading proton opening angle in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_9_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 2,  { -1.000, -0.800, -0.600, -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000 } },
    { 3,  { -1.000, -0.666, -0.333, 0.000, 0.333, 0.666, 1.000 } },
    { 4,  { -1.000, -0.500, 0.000, 0.500, 1.000 } },
    { 10, { -1.000, -0.500, 0.000, 0.500, 1.000 } },

  }; // v3

  std::string block9_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_cos_theta_mu_li_vec[1];";
  std::string block9_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[1];";
  std::string block9_title = "p_{mult};;cos #theta_{l_{2}#mu};";
  std::string block9_textitle = "p_{mult};;\\cos \\theta_{l_{2}\\mu};";

  Block2D* b9t = new Block2D(block9_true_branchexpr, block9_title, block9_textitle, BLOCK_9_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b9r = new Block2D(block9_reco_branchexpr, block9_title, block9_textitle, BLOCK_9_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b9t, b9r );

  // Block 10: third-leading proton opening angle in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_10_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 3,  { -1.000, -0.500, 0.000, 0.500, 1.000 } },
    { 4,  { -1.000, 0.000, 1.000 } },
    { 10, { -1.000, 0.000, 1.000 } },

  }; // v3

  std::string block10_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_cos_theta_mu_li_vec[2];";
  std::string block10_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[2];";
  std::string block10_title = "p_{mult};;cos #theta_{l_{3}#mu};";
  std::string block10_textitle = "p_{mult};;\\cos \\theta_{l_{3}\\mu};";

  Block2D* b10t = new Block2D(block10_true_branchexpr, block10_title, block10_textitle, BLOCK_10_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b10r = new Block2D(block10_reco_branchexpr, block10_title, block10_textitle, BLOCK_10_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b10t, b10r );

  // Block 11: fourth-leading proton opening angle in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_11_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 4,  { -1.000, 0.000, 1.000 } },
    { 10, { -1.000, 0.000, 1.000 } },

  }; // v3

  std::string block11_true_branchexpr = "Sum$(CC1muNp0piEnriched_true_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_true_cos_theta_mu_li_vec[3];";
  std::string block11_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[3];";
  std::string block11_title = "p_{mult};;cos #theta_{l_{4}#mu};";
  std::string block11_textitle = "p_{mult};;\\cos \\theta_{l_{4}\\mu};";

  Block2D* b11t = new Block2D(block11_true_branchexpr, block11_title, block11_textitle, BLOCK_11_BIN_EDGES, signal, kSignalTrueBin);
  Block2D* b11r = new Block2D(block11_reco_branchexpr, block11_title, block11_textitle, BLOCK_11_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b11t, b11r ); */

  /* -------------------------------------------------------------------------- */
  /*                                1D bin blocks                               */
  /* -------------------------------------------------------------------------- */

  // Don't add 1D blocks for now...

  /* // Block 1b: leading proton momentum (1D block)
  std::vector< double > BLOCK_1B_BIN_EDGES = { 0.250, 0.325, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 1.000 };

  std::string block1b_true_branchexpr = "CC1muNp0piEnriched_true_p3_lead_p.Mag();GeV/c";
  std::string block1b_reco_branchexpr = "CC1muNp0piEnriched_reco_p3_lead_p.Mag();GeV/c";
  std::string block1b_title = "p_{l_{1}};GeV/c";
  std::string block1b_textitle = "p_{l_{1}};GeV/c";

  Block1D* b1bt = new Block1D(block1b_true_branchexpr, block1b_title, block1b_textitle, BLOCK_1B_BIN_EDGES, signal, kSignalTrueBin);
  Block1D* b1br = new Block1D(block1b_reco_branchexpr, block1b_title, block1b_textitle, BLOCK_1B_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b1bt, b1br );

  // Block 2b: opening angle between muon and leading proton (1D block)
  std::vector< double > BLOCK_2B_BIN_EDGES = { -1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00 };

  std::string block2b_true_branchexpr = "CC1muNp0piEnriched_true_cos_theta_mu_li_vec[0];";
  std::string block2b_reco_branchexpr = "CC1muNp0piEnriched_true_cos_theta_mu_li_vec[0];";
  std::string block2b_title = "cos #theta_{l_{1}#mu};";
  std::string block2b_textitle = "\\cos \\theta_{l_{1}\\mu};";

  Block1D* b2bt = new Block1D(block2b_true_branchexpr, block2b_title, block2b_textitle, BLOCK_2B_BIN_EDGES, signal, kSignalTrueBin);
  Block1D* b2br = new Block1D(block2b_reco_branchexpr, block2b_title, block2b_textitle, BLOCK_2B_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b2bt, b2br );

  //Block 3b: available energy (1D block)
  std::vector< double > BLOCK_3B_BIN_EDGES = { 0.000, 0.050, 0.100, 0.200, 0.300, 0.400, 0.600 };

  std::string block3b_true_branchexpr = "CC1muNp0piEnriched_true_e_had_vis;GeV";
  std::string block3b_reco_branchexpr = "CC1muNp0piEnriched_reco_e_had_vis;GeV";
  std::string block3b_title = "E_{avail};GeV";
  std::string block3b_textitle = "E_{avail};GeV";

  Block1D* b3bt = new Block1D(block3b_true_branchexpr, block3b_title, block3b_textitle, BLOCK_3B_BIN_EDGES, signal, kSignalTrueBin);
  Block1D* b3br = new Block1D(block3b_reco_branchexpr, block3b_title, block3b_textitle, BLOCK_3B_BIN_EDGES, selection, kOrdinaryRecoBin);

  vect_block.emplace_back( b3bt, b3br ); */

  /* -------------------------------------------------------------------------- */
  /*                         Side bands definitions (1D)                        */
  /* -------------------------------------------------------------------------- */

  /* std::vector< double > PROTON_1D_BIN_EDGES = { 0.250, 0.325, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 1.000 };

  std::string side_branchexpr = "CC1muNp0piEnriched_reco_p3_p_vec.Mag();GeV/c";
  std::string side_title = "p_{l_{1}};GeV/c";
  std::string side_textitle = "p_{l_{1}};GeV/c";

  // Dirt sideband
  const std::string DIRT_SIDEBAND_SELECTION =
    "!" + selection_name_ + "_reco_vertex_in_FV && " + selection_name_ + "_pfp_starts_in_PCV"
    " && " + selection_name_ + "_has_muon_candidate && " + selection_name_ + "_topo_cut_passed"
    " && " + selection_name_ + "_no_reco_showers && " + selection_name_ + "_muon_passed_mom_cuts"
    " && " + selection_name_ + "_has_p_candidate && " + selection_name_ + "_passed_proton_pid_cut"
    " && " + selection_name_ + "_protons_contained && " + selection_name_ + "_lead_p_passed_mom_cuts";

  Block1D* b1s = new Block1D(side_branchexpr, side_title, side_textitle, PROTON_1D_BIN_EDGES, DIRT_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b1s );

  // NC sideband
  const std::string NC_SIDEBAND_SELECTION =
    selection_name_ + "_reco_vertex_in_FV && " + selection_name_ + "_pfp_starts_in_PCV"
    " && !" + selection_name_ + "_has_muon_candidate && " + selection_name_ + "_topo_cut_passed"
    " && " + selection_name_ + "_no_reco_showers && " + selection_name_ + "_has_p_candidate"
    " && " + selection_name_ + "_passed_proton_pid_cut && " + selection_name_ + "_protons_contained"
    " && " + selection_name_ + "_lead_p_passed_mom_cuts";

  Block1D* b2s = new Block1D(side_branchexpr, side_title, side_textitle, PROTON_1D_BIN_EDGES, NC_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b2s );

  // CCNpi sideband
  const std::string CCNPI_SIDEBAND_SELECTION =
    selection_name_ + "_reco_vertex_in_FV && " + selection_name_ + "_pfp_starts_in_PCV"
    " && " + selection_name_ + "_has_muon_candidate && " + selection_name_ + "_topo_cut_passed"
    " && " + selection_name_ + "_no_reco_showers && " + selection_name_ + "_muon_passed_mom_cuts"
    " && " + selection_name_ + "_has_p_candidate && " + selection_name_ + "_protons_contained"
    " && " + selection_name_ + "_lead_p_passed_mom_cuts"
    " && !" + selection_name_ + "_passed_proton_pid_cut";
    //" && trk_llr_pid_score_v[ " + selection_name_ + "_lead_p_candidate_idx ] > 0.2";

  Block1D* b3s = new Block1D(side_branchexpr, side_title, side_textitle, PROTON_1D_BIN_EDGES, CCNPI_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b3s ); */

  /* -------------------------------------------------------------------------- */
  /*                         Side bands definitions (2D)                        */
  /* -------------------------------------------------------------------------- */

  /* const std::string DIRT_SIDEBAND_SELECTION =
    "!" + selection_name_ + "_reco_vertex_in_FV && " + selection_name_ + "_pfp_starts_in_PCV"
    " && " + selection_name_ + "_has_muon_candidate && " + selection_name_ + "_topo_cut_passed"
    " && " + selection_name_ + "_no_reco_showers && " + selection_name_ + "_muon_passed_mom_cuts"
    " && " + selection_name_ + "_has_p_candidate && " + selection_name_ + "_passed_proton_pid_cut"
    " && " + selection_name_ + "_protons_contained && " + selection_name_ + "_lead_p_passed_mom_cuts";

  const std::string NC_SIDEBAND_SELECTION =
    selection_name_ + "_reco_vertex_in_FV && " + selection_name_ + "_pfp_starts_in_PCV"
    " && !" + selection_name_ + "_has_muon_candidate && " + selection_name_ + "_topo_cut_passed"
    " && " + selection_name_ + "_no_reco_showers && " + selection_name_ + "_has_p_candidate"
    " && " + selection_name_ + "_passed_proton_pid_cut && " + selection_name_ + "_protons_contained"
    " && " + selection_name_ + "_lead_p_passed_mom_cuts";

  const std::string CCNPI_SIDEBAND_SELECTION =
    selection_name_ + "_reco_vertex_in_FV && " + selection_name_ + "_pfp_starts_in_PCV"
    " && " + selection_name_ + "_has_muon_candidate && " + selection_name_ + "_topo_cut_passed"
    " && " + selection_name_ + "_no_reco_showers && " + selection_name_ + "_muon_passed_mom_cuts"
    " && " + selection_name_ + "_has_p_candidate && " + selection_name_ + "_protons_contained"
    " && " + selection_name_ + "_lead_p_passed_mom_cuts"
    " && !" + selection_name_ + "_passed_proton_pid_cut"; */

  const std::string COMBINED_SIDEBAND_SELECTION = 
    "CC1muNp0piEnriched_pfp_starts_in_PCV && "
    "CC1muNp0piEnriched_topo_cut_passed && "
    "CC1muNp0piEnriched_no_reco_showers && "
    "CC1muNp0piEnriched_has_p_candidate && "
    "CC1muNp0piEnriched_protons_contained && "
    "CC1muNp0piEnriched_lead_p_passed_mom_cuts && "
    "((!CC1muNp0piEnriched_reco_vertex_in_FV && CC1muNp0piEnriched_has_muon_candidate && CC1muNp0piEnriched_muon_passed_mom_cuts && CC1muNp0piEnriched_passed_proton_pid_cut) || "
    "(CC1muNp0piEnriched_reco_vertex_in_FV && !CC1muNp0piEnriched_has_muon_candidate && CC1muNp0piEnriched_passed_proton_pid_cut) || "
    "(CC1muNp0piEnriched_reco_vertex_in_FV && CC1muNp0piEnriched_has_muon_candidate && CC1muNp0piEnriched_muon_passed_mom_cuts && !CC1muNp0piEnriched_passed_proton_pid_cut))";

  // Block 1: proton multiplicity in slices of leading proton momentum (2D block)

  std::map< double, std::vector<double> > BLOCK_1_BIN_EDGES = {

    { 0.250, { 1, 2, 10} },
    { 0.375, { 1, 2, 3, 10} },
    { 0.475, { 1, 2, 3, 4, 10} },
    { 0.575, { 1, 2, 3, 4, 10} },
    { 0.650, { 1, 2, 3, 4, 10} },
    { 0.725, { 1, 2, 3, 4, 10} },
    { 0.800, { 1, 2, 3, 10} },
    { 0.850, { 1, 2, 3, 10} },
    { 1.000, { 1, 2, 3, 10} },
  
  };

  std::string block1_reco_branchexpr = "CC1muNp0piEnriched_reco_p_li_vec[0];GeV/c;Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);";
  std::string block1_title = "p_{l_{1}};GeV/c;p_{mult};";
  std::string block1_textitle = "p_{l_{1}};GeV/c;p_{mult};";

  Block2D* b1s = new Block2D(block1_reco_branchexpr, block1_title, block1_textitle, BLOCK_1_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b1s );

  // Block 2: proton multiplicity in slices of opening angle between muon and leading proton (2D block)

  std::map< double, std::vector<double> > BLOCK_2_BIN_EDGES = {

    { -1.00, { 1, 2, 3, 10} },
    { -0.85, { 1, 2, 3, 10} },
    { -0.65, { 1, 2, 3, 4, 10} },
    {  0.00, { 1, 2, 3, 4, 10} },
    {  0.65, { 1, 2, 3, 10} },
    {  0.85, { 1, 2, 3, 10} },
    {  1.00, { 1, 2, 3, 10} },
  
  };

  std::string block2_reco_branchexpr = "CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[0];;Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);";
  std::string block2_title = "cos #theta_{l_{1}#mu};;p_{mult};";
  std::string block2_textitle = "\\cos \\theta_{l_{1}\\mu};;p_{mult};";

  Block2D* b2s = new Block2D(block2_reco_branchexpr, block2_title, block2_textitle, BLOCK_2_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b2s );

  // Block 3: proton multiplicity in slices of available energy (2D block)

  std::map< double, std::vector<double> > BLOCK_3_BIN_EDGES = {

    { 0.050, { 1, 2, 10} },
    { 0.125, { 1, 2, 3, 10} },
    { 0.285, { 1, 2, 3, 10} },
    { 0.400, { 1, 2, 3, 4, 10} },
    { 1.000, { 1, 2, 3, 4, 10} },
  
  };

  std::string block3_reco_branchexpr = "CC1muNp0piEnriched_reco_e_had_vis;GeV;Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);";
  std::string block3_title = "E_{avail};GeV;p_{mult};";
  std::string block3_textitle = "E_{avail};GeV;p_{mult};";

  Block2D* b3s = new Block2D(block3_reco_branchexpr, block3_title, block3_textitle, BLOCK_3_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b3s );

  // Block 4: leading proton momentum in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_4_BIN_EDGES = {

    { 1,  { 0.250, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.900, 1.000 } },
    { 2,  { 0.250, 0.400, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 1.000 } },
    { 3,  { 0.250, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 1.000 } },
    { 4,  { 0.250, 0.600, 0.700, 0.800, 1.000 } },
    { 10, { 0.250, 0.600, 0.700, 0.800, 1.000 } },

  };

  std::string block4_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p_li_vec[0];GeV/c";
  std::string block4_title = "p_{mult};;p_{l_{1}};GeV/c";
  std::string block4_textitle = "p_{mult};;p_{l_{1}};GeV/c";

  Block2D* b4s = new Block2D(block4_reco_branchexpr, block4_title, block4_textitle, BLOCK_4_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b4s );

  // Block 5: subleading proton momentum in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_5_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 2,  { 0.250, 0.350, 0.400, 0.450, 0.500, 0.550, 1.000 } },
    { 3,  { 0.250, 0.450, 0.500, 0.550, 0.600, 0.650, 1.000 } },
    { 4,  { 0.250, 0.550, 0.600, 1.000 } },
    { 10, { 0.250, 0.550, 0.600, 1.000 } },
  
  };

  std::string block5_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p_li_vec[1];GeV/c";
  std::string block5_title = "p_{mult};;p_{l_{2}};GeV/c";
  std::string block5_textitle = "p_{mult};;p_{l_{2}};GeV/c";

  Block2D* b5s = new Block2D(block5_reco_branchexpr, block5_title, block5_textitle, BLOCK_5_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b5s );

  // Block 6: third-leading proton momentum in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_6_BIN_EDGES = {

    { 1,   { -10000.0, 10000.0 } },
    { 3,   { 0.250, 0.350, 0.400, 0.450, 0.500, 1.000 } },
    { 4,   { 0.250, 0.450, 0.500, 1.000 } },
    { 10,  { 0.250, 0.450, 0.500, 1.000 } },

  };

  std::string block6_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p_li_vec[2];GeV/c";
  std::string block6_title = "p_{mult};;p_{l_{3}};GeV/c";
  std::string block6_textitle = "p_{mult};;p_{l_{3}};GeV/c";

  Block2D* b6s = new Block2D(block6_reco_branchexpr, block6_title, block6_textitle, BLOCK_6_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b6s );

  // Block 7: fourth-leading proton momentum in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_7_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 4,  { 0.250, 0.400, 1.000 } },
    { 10, { 0.250, 0.400, 1.000 } },

  };
  
  std::string block7_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_p_li_vec[3];GeV/c";
  std::string block7_title = "p_{mult};;p_{l_{4}};GeV/c";
  std::string block7_textitle = "p_{mult};;p_{l_{4}};GeV/c";

  Block2D* b7s = new Block2D(block7_reco_branchexpr, block7_title, block7_textitle, BLOCK_7_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b7s );

  // Block 8: leading proton opening angle in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_8_BIN_EDGES = {

    { 1,  { -1.000, -0.800, -0.600, -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000 } },
    { 2,  { -1.000, -0.800, -0.600, -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000 } },
    { 3,  { -1.000, -0.750, -0.500, -0.250, 0.000, 0.250, 0.500, 0.750, 1.000 } },
    { 4,  { -1.000, -0.500, 0.000, 0.500, 1.000 } },
    { 10, { -1.000, -0.500, 0.000, 0.500, 1.000 } },

  };

  std::string block8_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[0];";
  std::string block8_title = "p_{mult};;cos #theta_{l_{1}#mu};";
  std::string block8_textitle = "p_{mult};;\\cos \\theta_{l_{1}\\mu};";

  Block2D* b8s = new Block2D(block8_reco_branchexpr, block8_title, block8_textitle, BLOCK_8_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b8s );

  // Block 9: subleading proton opening angle in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_9_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 2,  { -1.000, -0.800, -0.600, -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000 } },
    { 3,  { -1.000, -0.666, -0.333, 0.000, 0.333, 0.666, 1.000 } },
    { 4,  { -1.000, -0.500, 0.000, 0.500, 1.000 } },
    { 10, { -1.000, -0.500, 0.000, 0.500, 1.000 } },

  };

  std::string block9_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[1];";
  std::string block9_title = "p_{mult};;cos #theta_{l_{2}#mu};";
  std::string block9_textitle = "p_{mult};;\\cos \\theta_{l_{2}\\mu};";

  Block2D* b9s = new Block2D(block9_reco_branchexpr, block9_title, block9_textitle, BLOCK_9_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b9s );

  // Block 10: third-leading proton opening angle in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_10_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 3,  { -1.000, -0.500, 0.000, 0.500, 1.000 } },
    { 4,  { -1.000, 0.000, 1.000 } },
    { 10, { -1.000, 0.000, 1.000 } },

  };

  std::string block10_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[2];";
  std::string block10_title = "p_{mult};;cos #theta_{l_{3}#mu};";
  std::string block10_textitle = "p_{mult};;\\cos \\theta_{l_{3}\\mu};";

  Block2D* b10s = new Block2D(block10_reco_branchexpr, block10_title, block10_textitle, BLOCK_10_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b10s );

  // Block 11: fourth-leading proton opening angle in slices of proton multiplicity (2D block)

  std::map< double, std::vector<double> > BLOCK_11_BIN_EDGES = {

    { 1,  { -10000.0, 10000.0 } },
    { 4,  { -1.000, 0.000, 1.000 } },
    { 10, { -1.000, 0.000, 1.000 } },

  };

  std::string block11_reco_branchexpr = "Sum$(CC1muNp0piEnriched_reco_p3_p_vec.Mag() > 0.25);;CC1muNp0piEnriched_reco_cos_theta_mu_li_vec[3];";
  std::string block11_title = "p_{mult};;cos #theta_{l_{4}#mu};";
  std::string block11_textitle = "p_{mult};;\\cos \\theta_{l_{4}\\mu};";

  Block2D* b11s = new Block2D(block11_reco_branchexpr, block11_title, block11_textitle, BLOCK_11_BIN_EDGES, COMBINED_SIDEBAND_SELECTION, kSidebandRecoBin);

  vect_sideband.emplace_back( b11s );

  /* -------------------------------------------------------------------------- */
  /*                            Background categories                           */
  /* -------------------------------------------------------------------------- */

  // Add relevant background categories
  CATEGORY = selection_name_ + "_" + "EventCategory";
  //CATEGORY = "EventCategory";
  background_index = {1, 2, 3, 4, 17, 18, 19, 20, 21, 22};

}
