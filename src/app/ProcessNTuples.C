// Post-processing program for the MicroBooNE xsec_analyzer framework
//
// Updated 28 April 2025
// Steven Gardiner <gardiner@fnal.gov>
// Daniel Barrow <daniel.barrow@physics.ox.ac.uk>
// Patrick Green <patrick.green@physics.ox.ac.uk>
//
// Updated November 2026: the output is now a slim file that no longer copies
// the full input ntuple. See XSecAnalyzer/NtupleFormat.hh for the layout:
// one TTree per selection, plus an entry-aligned XSecAnalyzer/weights TTree
// holding is_mc, event IDs, the CV weight and the flattened systematic weight
// map. By default only events that are selected, or that are true signal (MC),
// for at least one requested selection are written; use --keep-all to write
// every processed event. --no-syst-weights omits the systematic weight
// branches (weight_<name>), keeping only the CV weight; combined with
// --keep-all this gives a compact all-events file for selection studies.

// Standard library includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TTree.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Constants.hh"
#include "XSecAnalyzer/Functions.hh"
#include "XSecAnalyzer/FiducialVolume.hh"
#include "XSecAnalyzer/NtupleFormat.hh"

#include "XSecAnalyzer/Selections/SelectionBase.hh"
#include "XSecAnalyzer/Selections/SelectionFactory.hh"
#include "XSecAnalyzer/TreeHandler.hh"

namespace {

  // Manages the XSecAnalyzer/weights output TTree
  class WeightsTreeWriter {

    public:

      WeightsTreeWriter( TDirectory* dir, TTree* pelee_tree, bool is_mc,
        bool syst_weights = true ) : is_mc_( is_mc )
      {
        dir->cd();
        tree_ = new TTree( "weights", "Event weights and IDs (aligned with"
          " the selection TTrees)" );
        tree_->Branch( "is_mc", &is_mc_ );
        tree_->Branch( "run", &run_ );
        tree_->Branch( "subrun", &subrun_ );
        tree_->Branch( "event", &event_ );
        tree_->Branch( "entry", &entry_ );

        if ( !pelee_tree ) return;
        pelee_ = pelee_tree;

        // Event IDs. These scalars are read through TLeaf::GetValue() so
        // that they work regardless of any branch addresses set elsewhere
        // (e.g., by a selection that also reads them).
        leaf_run_ = pelee_->GetLeaf( "run" );
        leaf_subrun_ = pelee_->GetLeaf( "sub" );
        leaf_event_ = pelee_->GetLeaf( "evt" );

        if ( !is_mc_ ) return;

        // CV weight (scalar) for quick plotting
        leaf_cv_weight_ = pelee_->GetLeaf( "weightSplineTimesTune" );
        if ( leaf_cv_weight_ ) {
          tree_->Branch( "weightSplineTimesTune", &cv_weight_ );
        }

        // Systematic weight map. The output branches are created on the
        // first event, once the keys are known.
        if ( syst_weights && pelee_->GetBranch( "weights" ) ) {
          pelee_->SetBranchAddress( "weights", &in_map_ );
          has_map_ = true;
        }
      }

      void fill( long long entry ) {
        entry_ = entry;
        run_ = leaf_run_ ? static_cast< int >( leaf_run_->GetValue() ) : 0;
        subrun_ = leaf_subrun_ ? static_cast< int >( leaf_subrun_->GetValue() )
          : 0;
        event_ = leaf_event_ ? static_cast< int >( leaf_event_->GetValue() )
          : 0;
        cv_weight_ = leaf_cv_weight_ ? leaf_cv_weight_->GetValue() : 1.f;

        if ( has_map_ && in_map_ ) {
          if ( !branches_made_ ) this->make_branches();
          for ( auto& pair : out_vecs_ ) {
            const std::string& key = pair.first;
            auto& vec = *pair.second;
            auto iter = in_map_->find( key );
            if ( iter != in_map_->end() ) vec = iter->second;
            else vec.assign( sizes_.at( key ), 1. );
          }
        }

        tree_->Fill();
      }

      void write() { tree_->Write( "", TObject::kOverwrite ); }

    protected:

      void make_branches() {
        for ( const auto& pair : *in_map_ ) {
          const std::string& key = pair.first;
          auto& vec_ptr = out_vecs_[ key ];
          vec_ptr = std::make_unique< std::vector< double > >();
          sizes_[ key ] = pair.second.size();
          std::string br_name = "weight_" + key;
          tree_->Branch( br_name.c_str(), vec_ptr.get() );
        }
        branches_made_ = true;
        std::cout << "Writing " << out_vecs_.size() << " systematic weight"
          << " branches to " << WEIGHTS_TREE_NAME << '\n';
      }

      TTree* tree_ = nullptr;
      TTree* pelee_ = nullptr;

      bool is_mc_ = false;
      int run_ = 0, subrun_ = 0, event_ = 0;
      long long entry_ = 0;
      float cv_weight_ = 1.f;

      TLeaf* leaf_run_ = nullptr;
      TLeaf* leaf_subrun_ = nullptr;
      TLeaf* leaf_event_ = nullptr;
      TLeaf* leaf_cv_weight_ = nullptr;

      std::map< std::string, std::vector< double > >* in_map_ = nullptr;
      bool has_map_ = false;
      bool branches_made_ = false;
      std::map< std::string, std::unique_ptr< std::vector< double > > >
        out_vecs_;
      std::map< std::string, size_t > sizes_;
  };

}

int main( int argc, char* argv[] ) {

  // Parse the optional flags (may appear anywhere)
  bool keep_all = false;
  bool syst_weights = true;
  std::vector< std::string > args;
  for ( int a = 1; a < argc; ++a ) {
    std::string arg( argv[ a ] );
    if ( arg == "--keep-all" ) keep_all = true;
    else if ( arg == "--no-syst-weights" ) syst_weights = false;
    else args.push_back( arg );
  }

  if ( args.size() != 4u ) {
    std::cout << "Usage: " << argv[0]
      << " INPUT_NTUPLE_FILE FILE_TYPE SELECTION_NAMES OUTPUT_FILE"
      << " [--keep-all] [--no-syst-weights]\n"
      << "  By default only events that are selected, or true signal (MC),"
      << " for at least one\n  selection are written. --keep-all writes every"
      << " processed event.\n  --no-syst-weights keeps only the CV weight"
      << " (no weight_<name> branches).\n";
    return 1;
  }

  std::string input_file_name( args[ 0 ] );
  std::string file_type( args[ 1 ] );
  std::stringstream sel_ss( args[ 2 ] );
  std::string output_file_name( args[ 3 ] );

  std::vector< std::string > selection_names;

  std::string sel_name;
  while ( std::getline( sel_ss, sel_name, ',' ) ) {
    selection_names.push_back( sel_name );
  }

  std::cout << "\nRunning ProcessNTuples with options:\n";
  std::cout << "\tinput_file_name: " << input_file_name << '\n';
  std::cout << "\tinput_file_type: " << file_type << '\n';
  std::cout << "\toutput_file_name: " << output_file_name << '\n';
  std::cout << "\tkeep: " << ( keep_all ? "all events"
    : "selected or true-signal events only" ) << '\n';
  std::cout << "\tsystematic weights: " << ( syst_weights ? "yes" : "no (CV only)" )
    << '\n';
  std::cout << "\n\nselection names:\n";
  for ( const auto& sel_name : selection_names ) {
    std::cout << "\t\t- " << sel_name << '\n';
  }

  // Open the input ntuple (read only) and create a new, slim output file
  TFile in_file( input_file_name.c_str(), "read" );
  if ( in_file.IsZombie() ) {
    throw std::runtime_error( "Could not open input file "
      + input_file_name );
  }

  TFile out_file( output_file_name.c_str(), "recreate" );
  if ( out_file.IsZombie() ) {
    throw std::runtime_error( "Could not create output file "
      + output_file_name );
  }
  TDirectory* sub_dir = out_file.mkdir( XSEC_DIR_NAME.c_str() );
  if ( !sub_dir ) {
    throw std::runtime_error( "Failed to create XSecAnalyzer subdirectory" );
  }
  sub_dir->cd();

  // Get the PeLEE TTree containing subrun POT information. Note that real data
  // does not have this TTree, so just skip summing the POT if it is missing.
  TTree* subruns_tree = nullptr;
  in_file.GetObject( "nuselection/SubRun", subruns_tree );

  // Get the total POT from the subruns TTree
  float summed_pot = 0.;
  if ( subruns_tree ) {
    TBranch* pot_br = subruns_tree->GetBranch( "pot" );
    if ( pot_br ) {
      float pot;
      pot_br->SetAddress( &pot );
      for ( int se = 0; se < pot_br->GetEntries(); ++se ) {
        pot_br->GetEntry( se );
        summed_pot += pot;
      }
    }
  }

  // Save the total POT in the output subdirectory as a TParameter< float >
  auto summed_pot_param = std::make_unique< TParameter< float > >(
    "summed_pot", summed_pot );

  sub_dir->cd();
  summed_pot_param->Write();

  // Determine whether we are working with data or MC using the presence
  // or absence of the SubRun TTree
  bool is_mc = ( subruns_tree != nullptr );

  // Also save the result as a TParameter< bool >
  auto is_mc_param = std::make_unique< TParameter< bool > >( "is_mc", is_mc );
  is_mc_param->Write();

  // Instantiate the requested selections to use for post-processing
  std::vector< std::unique_ptr< SelectionBase > > selections;

  SelectionFactory sf;
  for ( const auto& sel_name : selection_names ) {
    selections.emplace_back().reset( sf.CreateSelection( sel_name ) );
  }

  // Create a TreeHandler object to manage TTree input/output. Apply
  // setup procedures to it from each of the configured selections
  TreeHandler th;

  for ( auto& sel : selections ) {

    // Configure an output TTree to store the selection results, and give
    // it the same name as the selection object itself
    const std::string& sel_name = sel->name();
    th.add_output_tree( sub_dir, sel_name, sel_name + " selection results" );

    // Load all input TTrees required by the current selection object
    const auto& in_tree_names = sel->input_tree_names();
    for ( const auto& tn : in_tree_names ) {
      TTree* temp_tree = nullptr;
      in_file.GetObject( tn.c_str(), temp_tree );

      if ( !temp_tree ) throw std::runtime_error( "Could not load TTree \""
        + tn + "\" from the file " + input_file_name + " (this TTree"
        + " is required to apply the \"" + sel_name + "\" selection)" );

      th.add_input_tree( temp_tree, tn );
    }
  }

  // Event weights and IDs are taken from the PeLEE TTree if it is present
  const std::string PELEE_TREE_NAME = "nuselection/NeutrinoSelectionFilter";
  TTree* pelee_tree = nullptr;
  in_file.GetObject( PELEE_TREE_NAME.c_str(), pelee_tree );
  WeightsTreeWriter weights_writer( sub_dir, pelee_tree, is_mc,
    syst_weights );

  // EVENT LOOP
  long long events_entry = 0;
  long long num_written = 0;
  while ( true ) {

    if ( events_entry % 1000 == 0 ) {
      std::cout << "Processing event #" << events_entry << '\n';
    }

    // If the return value is negative, there was an I/O error, or we've
    // attempted to read past the end of a TTree
    int local_entry = th.load_tree( events_entry );

    // If we've reached the end of the TTrees (or encountered an I/O error),
    // then terminate the event loop
    if ( local_entry < 0 ) break;

    // Load the current event using branches from the input TTree(s)
    th.get_entry( events_entry );

    // Implement Patrick's recipe for event filtering to avoid double-counting
    // when integrating signal-enhanced samples
    // NOTE: The current treatment is entirely PeLEE-ntuple-dependent
    // TODO: Revisit how this is handled to make it more general/configurable
    // Get a pointer to the PeLEE ntuple TTree if it is available
    auto p_tree = th.in_tree( PELEE_TREE_NAME );

    // Only do the filtering at all if we can access the PeLEE ntuple
    if ( p_tree ) {

      // Active volume definition needed to correctly incorporate
      // signal-enhanced samples generated only in the active volume rather
      // than the full cryostat volume
      FiducialVolume av = { 0.0, 256.0, -120.0, 120.0, 0.0, 1076.0 };

      // Get access to some truth information for the current event that
      // we need for filtering
      auto p_map = th.in_map( PELEE_TREE_NAME );

      int nu_pdg, ccnc;
      float true_nu_x, true_nu_y, true_nu_z;
      p_map.at( "nu_pdg" ) >> nu_pdg;
      p_map.at( "ccnc" ) >> ccnc;
      p_map.at( "true_nu_vtx_x" ) >> true_nu_x;
      p_map.at( "true_nu_vtx_y" ) >> true_nu_y;
      p_map.at( "true_nu_vtx_z" ) >> true_nu_z;

      bool is_nue_or_nuebar = ( std::abs( nu_pdg ) == 12 );
      bool is_cc = ( ccnc == CHARGED_CURRENT );
      bool is_in_active_vol = av.is_inside( true_nu_x, true_nu_y, true_nu_z );
      bool is_nue_or_nuebar_cc_av = is_nue_or_nuebar
        && is_cc && is_in_active_vol;

      // Note that the signal-enhanced samples are typically generated only in
      // the active volume compared with full overlay that is generated for the
      // whole cryostat, and they may only be generated for CC events,
      // excluding NC

      // *** Intrinsic Nue ***
      if ( file_type == "nueMC" || file_type == "nueDV" ) {
        // Inverse cut, to avoid any accidental double-counting
        if ( !is_nue_or_nuebar_cc_av ) {
          ++events_entry;
          continue;
        }
      }
      else if ( file_type == "numuMC" && is_nue_or_nuebar_cc_av ) {
        ++events_entry;
        continue;
      }

      // *** Add any other signal enhanced samples here ***

      // NuMI specific: configure normalisation weight for dirt scaling
      if ( useNuMI ) {
        // SG: This is a bit of a hack, but it should work fine. The
        // normalization weight isn't an actual branch of the input
        // ntuple, but we can add its value to the input TreeMap anyway.
        // Since it's an input TreeMap and the branch hasn't been marked
        // as a variable-size array, the branch address for our new
        // entry will never be set, avoiding any conflict with the real
        // TTree contents. Any selection that wants to use it, however,
        // will need to make an actual branch to store the normalization
        // weight value in the output.
        // NOTE: I adopt the British spelling below to match the convention
        // in Patrick's code.
        float norm_wgt = 1.0f;
        if ( file_type == "dirtMC" ) norm_wgt = 0.65f;
        p_map[ "normalisation_weight" ] = norm_wgt;
      }

    } // Filtering done only if the PeLEE ntuple TTree is present

    bool keep = keep_all;
    for ( auto& sel : selections ) {
      // Clear any internal state as needed before processing the new event
      sel->reset();

      AnalysisEvent ev = th.get_event( sel->input_tree_names(), sel->name() );
      sel->apply_selection( is_mc, ev );

      keep = keep || sel->last_selected() || ( is_mc && sel->last_mc_signal() );
    }

    // We're done. Save the results (if requested) and move on to the next
    // event. All output TTrees are filled together, so they stay aligned.
    if ( keep ) {
      th.fill();
      weights_writer.fill( events_entry );
      ++num_written;
    }
    ++events_entry;
  }

  for ( auto& sel : selections ) {
    std::cout << sel->name() << " has " << sel->passed_events()
      << " events which passed\n";
  }
  std::cout << "Wrote " << num_written << " of " << events_entry
    << " processed events\n";
  std::cout << "Wrote output to: " << output_file_name << '\n';

  sub_dir->cd();
  th.write();
  weights_writer.write();
  out_file.Close();
}
