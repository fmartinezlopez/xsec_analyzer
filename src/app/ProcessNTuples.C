// Post-processing program for the MicroBooNE xsec_analyzer framework. This is
// currently designed for use with the PeLEE group's "searchingfornues" ntuples
//
// Updated 24 September 2024
// Steven Gardiner <gardiner@fnal.gov>
// Daniel Barrow <daniel.barrow@physics.ox.ac.uk>

// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/AnalysisEvent.hh"
#include "XSecAnalyzer/Branches.hh"
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/FiducialVolume.hh"

#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"

void analyze( const std::string& input_filename,
  const std::string& file_type,
  const std::vector< std::string >& selection_names,
  const std::string& output_filename )
{
  std::cout << "\nRunning ProcessNTuples with options:\n";
  std::cout << "\tinput_filename: " << input_filename << '\n';
  std::cout << "\tinput_file_type: " << file_type << '\n';
  std::cout << "\toutput_filename: " << output_filename << '\n';
  std::cout << "\n\nselection names:\n";
  for ( const auto& sel_name : selection_names ) {
    std::cout << "\t\t- " << sel_name << '\n';
  }

  // Get the TTrees containing the event ntuples and subrun POT information
  // Use TChain objects for simplicity in manipulating multiple files
  TChain events_ch( "nuselection/NeutrinoSelectionFilter" );
  TChain subruns_ch( "nuselection/SubRun" );
  events_ch.Add( input_filename.c_str() );
  subruns_ch.Add( input_filename.c_str() );

  // OUTPUT TTREE
  // Make an output TTree for plotting (one entry per event)
  TFile* out_file = new TFile( output_filename.c_str(), "recreate" );
  out_file->cd();
  TTree* out_tree = new TTree( "stv_tree", "STV analysis tree" );

  // Get the total POT from the subruns TTree. Save it in the output
  // TFile as a TParameter<float>. Real data doesn't have this TTree,
  // so check that it exists first.
  float pot;
  float summed_pot = 0.;
  bool has_pot_branch = ( subruns_ch.GetBranch("pot") != nullptr );
  if ( has_pot_branch ) {
    subruns_ch.SetBranchAddress( "pot", &pot );
    for ( int se = 0; se < subruns_ch.GetEntries(); ++se ) {
      subruns_ch.GetEntry( se );
      summed_pot += pot;
    }
  }

  TParameter<float>* summed_pot_param = new TParameter<float>( "summed_pot",
    summed_pot );

  summed_pot_param->Write();

  std::vector< std::unique_ptr<SelectionBase> > selections;

  SelectionFactory sf;
  for ( const auto& sel_name : selection_names ) {
    selections.emplace_back().reset( sf.CreateSelection(sel_name) );
  }

  out_file->cd();
  for ( auto& sel : selections ) {
    sel->setup( out_tree );
  }

  // Active volume definition
  // required for correctly incorporating signal enhanced samples 
  // generated only in active volume rather than full cryostat volume
  FiducialVolume AV = { 0.0, 256.0, -120.0, 120.0, 0.0, 1076.0 };

  // EVENT LOOP
  // TChains can potentially be really big (and spread out over multiple
  // files). When that's the case, calling TChain::GetEntries() can be very
  // slow. I get around this by using a while loop instead of a for loop.
  bool created_output_branches = false;
  long events_entry = 0;

  while ( true ) {

    //if ( events_entry > 1000) break;

    if ( events_entry % 1000 == 0 ) {
      std::cout << "Processing event #" << events_entry << '\n';
    }

    // Create a new AnalysisEvent object. This will reset all analysis
    // variables for the current event.
    AnalysisEvent cur_event;

    // Set branch addresses for the member variables that will be read
    // directly from the Event TTree.
    set_event_branch_addresses( events_ch, cur_event );

    // TChain::LoadTree() returns the entry number that should be used with
    // the current TTree object, which (together with the TBranch objects
    // that it owns) doesn't know about the other TTrees in the TChain.
    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of the TChain.
    int local_entry = events_ch.LoadTree( events_entry );

    // If we've reached the end of the TChain (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load all of the branches for which we've called
    // TChain::SetBranchAddress() above
    events_ch.GetEntry( events_entry );

    // Handle integrating signal enhanced samples
    // note that these are typically generated only in the active volume
    // compared with full overlay that is generated for the whole cryostat
    // and may only be generated for CC events, excluding NC
    
    // *** Intrinsic Nue ***
    if (file_type == "nueMC" || file_type == "nueDV") {
      // inverse cut, to avoid any accidental double-counting
      if ( !(std::abs(cur_event.mc_nu_pdg_) == 12 && cur_event.mc_nu_ccnc_ == 0 && point_inside_FV(AV, cur_event.mc_nu_vx_, cur_event.mc_nu_vy_, cur_event.mc_nu_vz_)) ) {
        ++events_entry;
        continue;
      }
    }
    if (file_type == "numuMC") {
      if ( (std::abs(cur_event.mc_nu_pdg_) == 12 && cur_event.mc_nu_ccnc_ == 0 && point_inside_FV(AV, cur_event.mc_nu_vx_, cur_event.mc_nu_vy_, cur_event.mc_nu_vz_)) ) {
        ++events_entry;
        continue;
      }
    }

    // *** Add any other signal enhanced samples here ***

    // NuMI specific: configure normalisation weight
    // dirt scaling
    if (useNuMI) {
      if (file_type == "dirtMC") cur_event.normalisation_weight_ = 0.65;
      else cur_event.normalisation_weight_ = 1.0;
    }

    // Set the output TTree branch addresses, creating the branches if needed
    // (during the first event loop iteration)
    bool create_them = false;
    if ( !created_output_branches ) {
      create_them = true;
      created_output_branches = true;
    }
    set_event_output_branch_addresses(*out_tree, cur_event, create_them );

    for ( auto& sel : selections ) {
      sel->apply_selection( &cur_event );
    }

    // We're done. Save the results and move on to the next event.
    out_tree->Fill();
    ++events_entry;
  }

  for ( auto& sel : selections ) {
    sel->summary();
  }
  std::cout << "Wrote output to:" << output_filename << std::endl;

  for ( auto& sel : selections ) {
    sel->final_tasks();
  }

  out_tree->Write();
  out_file->Close();
  delete out_file;
}

int main( int argc, char* argv[] ) {

  if ( argc != 5 ) {
    std::cout << "Usage: " << argv[0]
      << " INPUT_PELEE_NTUPLE_FILE FILE_TYPE SELECTION_NAMES OUTPUT_FILE\n";
    return 1;
  }

  std::string input_file_name( argv[1] );
  std::string output_file_name( argv[4] );

  std::vector< std::string > selection_names;

  std::stringstream sel_ss( argv[3] );
  std::string sel_name;
  while ( std::getline(sel_ss, sel_name, ',') ) {
    selection_names.push_back( sel_name );
  }

  std::string file_type( argv[2] );

  analyze( input_file_name, file_type, selection_names, output_file_name );

  return 0;
}
