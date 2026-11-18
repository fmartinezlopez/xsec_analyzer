#pragma once

// Helpers describing the layout of the post-processed ntuple files written by
// ProcessNTuples, so that downstream code (UniverseMaker, univmake,
// SystematicsCalculator, MakeConfig) can work with both the current ("slim")
// format and the legacy "stv_tree" format.
//
// Slim format (written by ProcessNTuples):
//   XSecAnalyzer/summed_pot   TParameter<float>, simulated POT (MC only)
//   XSecAnalyzer/is_mc        TParameter<bool>
//   XSecAnalyzer/<selection>  TTree with the selection results (one per
//                             selection: <sel>_Selected, <sel>_MC_Signal,
//                             <sel>_Category, reco_ / true_ quantities, ...)
//   XSecAnalyzer/weights      TTree aligned entry by entry with the selection
//                             TTrees: is_mc, run, subrun, event, entry (in the
//                             input ntuple), weightSplineTimesTune and one
//                             std::vector<double> branch weight_<name> per
//                             entry of the PeLEE "weights" map
//
// Legacy format: a single TTree "stv_tree" at the top level holding all of the
// above, with the POT stored as "summed_pot" at the top level.

// Standard library includes
#include <string>

// ROOT includes
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"

const std::string XSEC_DIR_NAME = "XSecAnalyzer";
const std::string WEIGHTS_TREE_NAME = XSEC_DIR_NAME + "/weights";
const std::string LEGACY_TREE_NAME = "stv_tree";

// Returns the simulated POT stored in a post-processed ntuple file (either
// format), or a negative value if it cannot be found
inline double get_summed_pot( TFile& file ) {
  TParameter<float>* pot = nullptr;
  file.GetObject( ( XSEC_DIR_NAME + "/summed_pot" ).c_str(), pot );
  if ( !pot ) file.GetObject( "summed_pot", pot );
  if ( !pot ) return -1.;
  return pot->GetVal();
}

// Returns the TTree holding the event weights and the is_mc flag: the
// XSecAnalyzer/weights TTree for slim files, or stv_tree for legacy files.
// Returns nullptr if neither is present.
inline TTree* get_weights_tree( TFile& file ) {
  TTree* tree = nullptr;
  file.GetObject( WEIGHTS_TREE_NAME.c_str(), tree );
  if ( !tree ) file.GetObject( LEGACY_TREE_NAME.c_str(), tree );
  return tree;
}
