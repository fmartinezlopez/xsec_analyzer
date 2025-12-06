const std::string selection_name = "CC1muNp0piEnriched";

/* -------------------------------------------------------------------------- */
/*                      Utilities to handle input ntuples                     */
/* -------------------------------------------------------------------------- */

struct RunFile {
    const char* filename;      // full path
    int run;                   // run number
    const char* type;          // "onBNB", "extBNB", or "MC"
    double pot;                // POT for onBNB, 0 otherwise
    int triggers;              // triggers for on/off data
    bool is_reweightable_mc;   //
    bool is_detvar;            //
};

// Replace forward slashes by plus signs in path to use in TFile
TString folder_from_file_name(const std::string& file_name) {
    constexpr char FORWARD_SLASH = '/';
    constexpr char PLUS = '+';
    TString result = file_name.c_str();
    result.ReplaceAll(FORWARD_SLASH, PLUS);
    return result;
}

// Checks whether a given file represents an ntuple from a reweightable MC sample
Bool_t is_reweightable_mc_ntuple(const char* input_file_name, const char* input_tree_name="stv_tree",
                                 const char* branch_name="weight_TunedCentralValue_UBGenie") {
    // Open the input file
    TFile *temp_file = TFile::Open(input_file_name, "READ");
    if (!temp_file || temp_file->IsZombie()) {
        std::cerr << "Error opening file " << input_file_name << std::endl;
        return false;
    }

    // Get the input tree
    TTree *temp_tree = (TTree*)temp_file->Get(input_tree_name);
    if (!temp_tree) {
        std::cerr << "Error getting tree " << input_tree_name << std::endl;
        temp_file->Close();
        return false;
    }

    // Check if branch exists
    TBranch* cv_weight_br = temp_tree->GetBranch(branch_name);
    Bool_t has_cv_weights = (cv_weight_br != nullptr);
    return has_cv_weights;
}

// Parse a single line
// Returns true if the line is valid, false if ignored
bool parseLine(const std::string& line, RunFile &outFile) {
    std::istringstream iss(line);
    std::string filenameStr, typeStr;
    int runNum = 0;
    int triggers = 0;
    double pot = 0;
    bool is_detvar = false;

    iss >> filenameStr >> runNum >> typeStr;

    const char* finalType = nullptr;

    if (typeStr == "onBNB") {
        iss >> triggers >> pot;
        finalType = strdup("onBNB");
    } else if (typeStr == "extBNB") {
        iss >> triggers;
        finalType = strdup("extBNB");
    } else if (strstr(typeStr.c_str(), "MC") != nullptr) {
        finalType = strdup("MC");
    } else if (strstr(typeStr.c_str(), "detVar") != nullptr) {
        finalType = strdup(typeStr.c_str());
        is_detvar = true;
    } else {
        // unknown type -- ignore
        return false;
    }

    char* fname_c = strdup(filenameStr.c_str());
    bool is_reweightable_mc = is_reweightable_mc_ntuple(fname_c);

    outFile = {fname_c, runNum, finalType, pot, triggers, is_reweightable_mc, is_detvar};
    return true;
}

// Parse the full file list
std::vector<RunFile> parseFileList(const char* listFile) {
    std::ifstream infile(listFile);
    if (!infile.is_open()) {
        throw std::runtime_error("Cannot open file list");
    }

    std::vector<RunFile> files;
    std::string line;
    while (std::getline(infile, line)) {
        // skip empty lines or comments
        if (line.empty() || line[0]=='#') continue;

        RunFile f;
        if (parseLine(line, f)) {
            files.push_back(f);
        }
        // else ignored
    }

    return files;
}

/* -------------------------------------------------------------------------- */
/*                   Container for histograms and universes                   */
/* -------------------------------------------------------------------------- */

struct Universe {
    std::vector<Float_t> weights_vector;
    Float_t cv_correction;
    Bool_t average;
};

const Float_t MIN_WEIGHT = 0.;
const Float_t MAX_WEIGHT = 30.;

Float_t safe_weight(Float_t w) {
  if (std::isfinite(w) && w >= MIN_WEIGHT && w <= MAX_WEIGHT) return w;
  else return 1.0;
}

class CategoryHistogram {
private:
    // Histogram configuration
    TString histName;
    TString xAxisLabel;
    Int_t nXBins;
    Double_t xMin;
    Double_t xMax;
    Int_t nCategories;
    std::vector<TString> Categories;  // Use vector instead of array
    
    // Unweighted histogram
    TH2D* unweighted_histogram;
    TH1D* unweighted_histogram_mc;
    
    // Histogram universes -- don't need categories for these
    std::map<TString, Bool_t> universe_names;
    std::map<TString, TH1D*> universe_histograms;
    std::map<TString, TH2D*> covariance_matrices;

    // Detector variation histograms
    std::map<TString, TH1D*> detvar_histograms;
    std::map<TString, TH2D*> detvar_covariance_matrices;
    
public:
    CategoryHistogram(const TString& name, const TString& label,
                     Int_t nbins, Double_t min, Double_t max,
                     Int_t ncat, const std::vector<TString>& cat)
        : histName(name), xAxisLabel(label),
          nXBins(nbins), xMin(min), xMax(max),
          nCategories(ncat), Categories(cat)
    {
        // Create title string
        TString title = Form("%s;%s;Counts", xAxisLabel.Data(), xAxisLabel.Data());
        
        // Create TH2 histogram
        unweighted_histogram = new TH2D(histName, title,
                                        nXBins, xMin, xMax,
                                        nCategories, 0, nCategories);
        
        // Set y-axis bin labels
        for (int i = 0; i < nCategories; i++) {
            unweighted_histogram->GetYaxis()->SetBinLabel(i + 1, Categories[i]);
        }
        
        // Disassociate from current open file
        unweighted_histogram->SetDirectory(nullptr);
    }

    // Copy constructor
    CategoryHistogram(const CategoryHistogram& other)
        : histName(other.histName), xAxisLabel(other.xAxisLabel),
          nXBins(other.nXBins), xMin(other.xMin), xMax(other.xMax),
          nCategories(other.nCategories), Categories(other.Categories),
          universe_names(other.universe_names)
    {
        // Clone unweighted histogram
        unweighted_histogram = (TH2D*)other.unweighted_histogram->Clone();
        unweighted_histogram->SetDirectory(nullptr);
        
        // Clone universe histograms
        for (const auto& pair : other.universe_histograms) {
            TH1D* cloned_hist = (TH1D*)pair.second->Clone();
            cloned_hist->SetDirectory(nullptr);
            universe_histograms[pair.first] = cloned_hist;
        }

        // Clone detvar histograms
        for (const auto& pair : other.detvar_histograms) {
            TH1D* cloned_hist = (TH1D*)pair.second->Clone();
            cloned_hist->SetDirectory(nullptr);
            detvar_histograms[pair.first] = cloned_hist;
        }
    }
    
    ~CategoryHistogram() {
        delete unweighted_histogram;
        for (auto& pair : universe_histograms) {
            delete pair.second;
        }
    }

    inline TH2D* GetUnweightedHistogram() { return unweighted_histogram; };
    inline std::map<TString, Bool_t> GetUniverseNames() { return universe_names; };
    inline std::map<TString, TH1D*> GetUniverseHistograms() { return universe_histograms; };

    // Make a new universe histogram
    void CreateUniverseHistogram(TString univ) {

        // Create title string
        TString title = Form("%s %s;%s;Counts", xAxisLabel.Data(), univ.Data(), xAxisLabel.Data());
        
        // Create TH2 histogram
        TH1D* weighted_histogram = new TH1D(histName, title,
                                            nXBins, xMin, xMax);
        
        // Disassociate from current open file
        weighted_histogram->SetDirectory(nullptr);

        // Add to map
        universe_histograms[univ] = weighted_histogram;

    }

    // Make a new detvar histogram
    void CreateDetVarHistogram(TString detvar) {

        // Create title string
        TString title = Form("%s %s;%s;Counts", xAxisLabel.Data(), detvar.Data(), xAxisLabel.Data());
        
        // Create TH2 histogram
        TH1D* weighted_histogram = new TH1D(histName, title,
                                            nXBins, xMin, xMax);
        
        // Disassociate from current open file
        weighted_histogram->SetDirectory(nullptr);

        // Add to map
        detvar_histograms[detvar] = weighted_histogram;

    }

    // Fill unweighted histogram
    void FillUnweighted(Double_t x, Double_t y, Double_t w) {
        unweighted_histogram->Fill(x, y, safe_weight(w));
    }

    // Fill corresponding universe histogram with weight
    void FillUniverse(TString univ, Double_t x, Double_t w) {
        // Create histogram if it doesn't exist yet
        if (universe_histograms.find(univ) == universe_histograms.end()) {
            CreateUniverseHistogram(univ);
        }

        universe_histograms[univ]->Fill(x, safe_weight(w));
    }

    void FillUniversesFromMap(std::map<TString, Universe> weight_map, Double_t x) {
        for (const auto& [universe_name, universe] : weight_map) {
            if (universe_names.find(universe_name) == universe_names.end()) {
                universe_names[universe_name] = universe.average;
            }
            for (Int_t k=0; k < universe.weights_vector.size(); ++k) {
                FillUniverse(Form("%s_%d", universe_name.Data(), k),
                             x, safe_weight(universe.weights_vector[k] * universe.cv_correction));
            }
        }
    }

    // Fill corresponding detvar histogram with weight
    void FillDetVarUniverse(const char* detvar, Double_t x, Double_t w) {
        // Create histogram if it doesn't exist yet
        TString detvar_string(detvar);
        if (detvar_histograms.find(detvar_string) == detvar_histograms.end()) {
            CreateDetVarHistogram(detvar_string);
        }

        detvar_histograms[detvar_string]->Fill(x, safe_weight(w));
    }

    // Clone all histograms and reset them, returning a new CategoryHistogram
    CategoryHistogram CloneAndReset() const {
        CategoryHistogram cloned(*this);
        
        // Reset unweighted histogram
        cloned.unweighted_histogram->Reset();

        // Clear the universe_names map
        cloned.universe_names.clear();
        
        // Reset universe histograms
        for (auto& pair : cloned.universe_histograms) {
            pair.second->Reset();
        }

        // Reset detvar histograms
        for (auto& pair : cloned.detvar_histograms) {
            pair.second->Reset();
        }
        
        return cloned;
    }

    // Scale all histograms by a factor
    void Scale(Double_t scale_factor) {
        // Scale unweighted histogram
        unweighted_histogram->Scale(scale_factor);
        
        // Scale universe histograms
        for (auto& pair : universe_histograms) {
            pair.second->Scale(scale_factor);
        }

        // Scale detvar histograms
        for (auto& pair : detvar_histograms) {
            pair.second->Scale(scale_factor);
        }
    }

    // Add histograms from another CategoryHistogram to this one
    void Add(const CategoryHistogram& other, Double_t scale_factor = 1.0) {
        // Add unweighted histogram
        unweighted_histogram->Add(other.unweighted_histogram, scale_factor);

        // Merge universe_names maps
        for (const auto& [name, value] : other.universe_names) {
            universe_names[name] = value;
        }
        
        // Add universe histograms
        for (const auto& pair : other.universe_histograms) {
            TString universe = pair.first;
            
            // If this universe exists in current object, add to it
            if (universe_histograms.find(universe) != universe_histograms.end()) {
                universe_histograms[universe]->Add(pair.second, scale_factor);
            }
            // Otherwise, create a new histogram by cloning and scaling
            else {
                TH1D* new_hist = (TH1D*)pair.second->Clone();
                new_hist->Scale(scale_factor);
                new_hist->SetDirectory(nullptr);
                universe_histograms[universe] = new_hist;
            }
        }

        // Add detvar histograms
        for (const auto& pair : other.detvar_histograms) {
            TString detvar = pair.first;
            
            // If this detvar exists in current object, add to it
            if (detvar_histograms.find(detvar) != detvar_histograms.end()) {
                detvar_histograms[detvar]->Add(pair.second, scale_factor);
            }
            // Otherwise, create a new histogram by cloning and scaling
            else {
                TH1D* new_hist = (TH1D*)pair.second->Clone();
                new_hist->Scale(scale_factor);
                new_hist->SetDirectory(nullptr);
                detvar_histograms[detvar] = new_hist;
            }
        }
    }

    void MakeProjectionMC() {
        // Project unweighted histogram for MC categories
        // Assumes "EXT" and "Data" are the two last categories
        unweighted_histogram_mc = unweighted_histogram->ProjectionX(Form("%s_MC", histName.Data()),
                                                                    1, nCategories - 2);
    }

    TH2D* InitializeCovarianceMatrix(TString univ_name) {
        TH2D* cov_mat = new TH2D(Form("cov_mat_%s", univ_name.Data()), univ_name,
                                 unweighted_histogram->GetNbinsX(),
                                 0,
                                 unweighted_histogram->GetNbinsX(),
                                 unweighted_histogram->GetNbinsX(),
                                 0,
                                 unweighted_histogram->GetNbinsX());

        cov_mat->Reset();
        return cov_mat;
    }

    void FillCovarianceMatrix(TH2D* cov_mat, TH1D* detvar, TH1D* cv) {
        for (Int_t i = 1; i <= cv->GetNbinsX(); ++i) {

            Float_t cv_i = cv->GetBinContent(i);
            Float_t univ_i = detvar->GetBinContent(i);

            for (Int_t j = 1; j <= cv->GetNbinsX(); ++j) {

                Float_t cv_j = cv->GetBinContent(j);
                Float_t univ_j = detvar->GetBinContent(j);

                Float_t cov = (cv_i - univ_i) * (cv_j - univ_j);

                // The lower bound of each covariance matrix bin is the bin index
                // Filling using the zero-based bin indices and the covariance as
                // weight increments the existing element by current covariance value
                cov_mat->Fill(i-1, j-1, cov);
            }
        }
    }

    // Use unweighted MC histogram if no CV hist provided
    void FillCovarianceMatrix(TH2D* cov_mat, TH1D* univ) {
        FillCovarianceMatrix(cov_mat, univ, unweighted_histogram_mc);
    }

    void MakeAllMatrices() {
        for (const auto& univ : universe_names) {
            TString univ_name = univ.first;
            Int_t num_universes = 0;
            TH2D* cov_mat = InitializeCovarianceMatrix(univ_name);
            for (const auto& hist : universe_histograms) {
                TString hist_name = hist.first;
                if (hist_name.Contains(univ_name)) {
                    FillCovarianceMatrix(cov_mat, hist.second);
                    num_universes++;
                }
            } // end loop over histograms
            covariance_matrices[univ_name] = cov_mat;
            if (univ.second) {
                cov_mat->Scale(1.0 / num_universes);
            }
        } // end loop over universes
    }

    void MakeAllDetVarMatrices(std::map<TString, TString> detvar_map) {

        for (const auto& pair : detvar_map) {
            TString detvar_name = pair.first;
            TString cv_name = pair.second;
            TH2D* cov_mat = InitializeCovarianceMatrix(detvar_name);

            TH1D* detvar_hist = detvar_histograms[detvar_name];
            TH1D* cv_hist = detvar_histograms[cv_name];

            FillCovarianceMatrix(cov_mat, detvar_hist, cv_hist);

            detvar_covariance_matrices[detvar_name] = cov_mat;
        } // end loop over detVars

    }

    // Write all histograms to a TFile in a specified directory
    void Finalize(TFile* file, const TString& directory_name,
                  std::map<TString, TString> detvar_map,
                  Bool_t write_univ_hist=false,
                  Bool_t write_detvar_hist=false) {
        if (!file || !file->IsOpen()) {
            std::cerr << "Error: Invalid or closed TFile provided to Write()" << std::endl;
            return;
        }
        
        // Save current directory
        TDirectory* current_dir = gDirectory;
        
        // Navigate to the file
        file->cd();
        
        // Create directory if it doesn't exist, or navigate to it if it does
        TDirectory* target_dir = file->GetDirectory(directory_name);
        if (!target_dir) {
            target_dir = file->mkdir(directory_name);
        }
        target_dir->cd();
        
        // Write unweighted histogram
        unweighted_histogram->Write("unweighted_histogram_2d");

        // Project unweighted histogram for MC categories
        MakeProjectionMC();
        unweighted_histogram_mc->Write("unweighted_histogram_mc_1d");

        // Write universe histograms
        if (write_univ_hist) {
            for (const auto& pair : universe_histograms) {
                pair.second->Write(pair.first);
            }
        }

        // Write covariance matrices
        MakeAllMatrices();
        for (const auto& pair : covariance_matrices) {
            pair.second->Write(Form("CovMat_%s", pair.first.Data()));
        }

        // Write detvar histograms
        if (write_detvar_hist) {
            for (const auto& pair : detvar_histograms) {
                pair.second->Write(pair.first);
            }
        }

        // Write detVar covariance matrices
        MakeAllDetVarMatrices(detvar_map);
        for (const auto& pair : detvar_covariance_matrices) {
            pair.second->Write(Form("CovMat_%s", pair.first.Data()));
        }
        
        // Restore original directory
        current_dir->cd();
    }

};

// Particle categories
const Int_t N_PARTICLE_CATEGORIES = 10;
const std::vector<TString> PARTICLE_CATEGORIES = {"Cosmic", "e", "#mu", "#gamma", "#pi^{#pm}", "K^{#pm}", "n", "p", "EXT", "Data"};

Int_t pdg_to_category(const Int_t pdgCode) {
    if (pdgCode == 0) {
        return 0;
    } else if (abs(pdgCode) == 11) {
        return 1;
    } else if (abs(pdgCode) == 13) {
        return 2;
    } else if (pdgCode == 22) {
        return 3;
    } else if (abs(pdgCode) == 211) {
        return 4;
    } else if (abs(pdgCode) == 321) {
        return 5;
    } else if (pdgCode == 2112) {
        return 6;
    } else if (pdgCode == 2212) {
        return 7;
    } else {
        return 999;
    }
}

const Int_t EXT_category = 8;
const Int_t Data_category = 9;

/* -------------------------------------------------------------------------- */
/*                  Map matching detVar samples to CV ntuples                 */
/* -------------------------------------------------------------------------- */

std::map<TString, TString> detvar_map = {
                                             {"detVarLYrayl",    "detVarCV"},
                                             {"detVarLYatten",   "detVarCV"},
                                             {"detVarLYdown",    "detVarCVLYDown"},
                                             {"detVarWMAngleXZ", "detVarCV"},
                                             {"detVarWMAngleYZ", "detVarCV"},
                                             {"detVarWMX",       "detVarCV"},
                                             {"detVarWMYZ",      "detVarCV"},
                                             {"detVarWMdEdx",    "detVarCVdEdx"},
                                             {"detVarSCE",       "detVarCVExtra"},
                                             {"detVarRecomb2",   "detVarCVExtra"}
                                        };

/* -------------------------------------------------------------------------- */
/*                             Workhorse function                             */
/* -------------------------------------------------------------------------- */

void FillHistograms(RunFile file, CategoryHistogram& hist,
                    std::map<TString, Universe> weight_map,
                    Float_t x, Float_t y, Float_t w) {
    if (file.is_detvar) {
        hist.FillDetVarUniverse(file.type, x, w);
    } else {
        hist.FillUnweighted(x, y, w);
    }

    if (file.is_reweightable_mc) hist.FillUniversesFromMap(weight_map, x);
}

void selection_plots(const char* inputFileList, const char* outputFile, const char* treeName = "stv_tree") {

    // Create the output file
    TFile *fileOut = TFile::Open(outputFile, "RECREATE");
    if (!fileOut || fileOut->IsZombie()) {
        std::cerr << "Error creating output file" << std::endl;
        fileOut->Close();
        return;
    }

    // Input variable declaration

    // Boolean cuts
    Bool_t reco_vertex_in_FV;
    Bool_t pfp_starts_in_PCV;
    Bool_t topo_cut_passed;
    Bool_t numu_cc;
    Bool_t muon_passed_mom_cuts;
    Bool_t no_reco_showers;
    Bool_t has_p_candidate;
    Bool_t protons_contained;

    // PFParticles
    std::vector<Int_t>* pfp_generation_v = nullptr;
    std::vector<Float_t>* trk_score_v = nullptr;
    std::vector<Float_t>* backtracked_pdg = nullptr;
    std::vector<Float_t>* trk_distance_v = nullptr;
    std::vector<Float_t>* trk_len_v = nullptr;
    std::vector<Float_t>* trk_llr_pid_score_v = nullptr;

    // Weights
    Bool_t is_mc;
    Float_t spline_weight;
    Float_t tuned_cv_weight;
    std::vector<Float_t>* weight_flux_all = nullptr;
    std::vector<Float_t>* weight_reint_all = nullptr;
    std::vector<Float_t>* weight_All_UBGenie = nullptr;
    std::vector<Float_t>* weight_AxFFCCQEshape_UBGenie = nullptr;
    std::vector<Float_t>* weight_DecayAngMEC_UBGenie = nullptr;
    std::vector<Float_t>* weight_NormCCCOH_UBGenie = nullptr;
    std::vector<Float_t>* weight_NormNCCOH_UBGenie = nullptr;
    std::vector<Float_t>* weight_RPA_CCQE_UBGenie = nullptr;
    std::vector<Float_t>* weight_ThetaDelta2NRad_UBGenie = nullptr;
    std::vector<Float_t>* weight_Theta_Delta2Npi_UBGenie = nullptr;
    std::vector<Float_t>* weight_VecFFCCQEshape_UBGenie = nullptr;
    std::vector<Float_t>* weight_XSecShape_CCMEC_UBGenie = nullptr;
    std::vector<Float_t>* weight_xsr_scc_Fa3_SCC = nullptr;
    std::vector<Float_t>* weight_xsr_scc_Fv3_SCC = nullptr;

    // Declare histograms
    CategoryHistogram hTrackScoreMuon("hTrackScoreMuon", "Track Score",
                                      50, 0.0, 1.0,
                                      N_PARTICLE_CATEGORIES, PARTICLE_CATEGORIES);
    CategoryHistogram hTrackDistMuon("hTrackDistMuon", "Track Distance",
                                      50, 0.0, 15.0,
                                      N_PARTICLE_CATEGORIES, PARTICLE_CATEGORIES);
    CategoryHistogram hTrackLenMuon("hTrackLenMuon", "Track Length",
                                      50, 0.0, 100.0,
                                      N_PARTICLE_CATEGORIES, PARTICLE_CATEGORIES);
    CategoryHistogram hTrackPIDMuon("hTrackPIDMuon", "Track PID Score",
                                      50, -1.0, 1.0,
                                      N_PARTICLE_CATEGORIES, PARTICLE_CATEGORIES);

    CategoryHistogram hTrackPIDProton("hTrackPIDProton", "Track PID Score",
                                      50, -1.0, 1.0,
                                      N_PARTICLE_CATEGORIES, PARTICLE_CATEGORIES);

    // Parse all files
    auto files = parseFileList(inputFileList);

    // Create maps for run → POT and triggers
    double total_pot = 0.0;
    std::map<int,double> run_pot;
    std::map<int,int> run_triggers;

    for (auto& f : files) {
        if (strcmp(f.type, "onBNB") == 0) {
            run_pot[f.run] = f.pot;
            run_triggers[f.run] = f.triggers;
            total_pot += f.pot;
        } 
    }

    for (auto& f : files) {
    std::cout << "Run " << f.run
              << ", type " << f.type
              << ", filename " << f.filename
              << ", triggers " << f.triggers
              << ", POT " << f.pot
              << std::endl;
    }

    for (const auto& file : files) {

        // Open the input file
        TFile *fileIn = TFile::Open(file.filename, "READ");
        if (!fileIn || fileIn->IsZombie()) {
            std::cerr << "Error opening input file" << std::endl;
            return;
        }

        // Get the input tree
        TTree *treeIn = (TTree*)fileIn->Get(treeName);
        if (!treeIn) {
            std::cerr << "Error getting tree from input file" << std::endl;
            fileIn->Close();
            return;
        }

        Float_t ScaleFactor = 1.0;
        if (strcmp(file.type, "MC") == 0) {
            Float_t summed_pot = dynamic_cast<TParameter<Float_t>*>(fileIn->Get("summed_pot"))->GetVal();
            ScaleFactor = run_pot[file.run] / summed_pot;
        } else if (strcmp(file.type, "extBNB") == 0) {
            ScaleFactor = static_cast<Double_t>(run_triggers[file.run]) / static_cast<Double_t>(file.triggers);
        } else if (strstr(file.type, "detVar") != nullptr) {
            // Scale detVar to total POT
            Float_t summed_pot = dynamic_cast<TParameter<Float_t>*>(fileIn->Get("summed_pot"))->GetVal();
            ScaleFactor = total_pot / summed_pot;
        }

        std::cout << "    Scale factor: " << ScaleFactor << std::endl;

        // Selection cut booleans
        treeIn->SetBranchAddress((selection_name + "_reco_vertex_in_FV").c_str(), &reco_vertex_in_FV);
        treeIn->SetBranchAddress((selection_name + "_pfp_starts_in_PCV").c_str(), &pfp_starts_in_PCV);
        treeIn->SetBranchAddress((selection_name + "_topo_cut_passed").c_str(),     &topo_cut_passed);

        treeIn->SetBranchAddress((selection_name + "_nu_mu_cc").c_str(),                          &numu_cc);
        treeIn->SetBranchAddress((selection_name + "_muon_passed_mom_cuts").c_str(), &muon_passed_mom_cuts);
        treeIn->SetBranchAddress((selection_name + "_no_reco_showers").c_str(),           &no_reco_showers);
        treeIn->SetBranchAddress((selection_name + "_has_p_candidate").c_str(),           &has_p_candidate);
        treeIn->SetBranchAddress((selection_name + "_protons_contained").c_str(),       &protons_contained);

        // PFParticles branches
        treeIn->SetBranchAddress("pfp_generation_v",       &pfp_generation_v);
        treeIn->SetBranchAddress("trk_score_v",                 &trk_score_v);
        treeIn->SetBranchAddress("backtracked_pdg",         &backtracked_pdg);
        treeIn->SetBranchAddress("trk_distance_v",           &trk_distance_v);
        treeIn->SetBranchAddress("trk_len_v",                     &trk_len_v);
        treeIn->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v);

        // Weights branches
        treeIn->SetBranchAddress("is_mc",                                                        &is_mc);
        treeIn->SetBranchAddress("spline_weight",                                        &spline_weight);
        treeIn->SetBranchAddress("tuned_cv_weight",                                    &tuned_cv_weight);
        if (file.is_reweightable_mc) {
            treeIn->SetBranchAddress("weight_flux_all",                                &weight_flux_all);
            treeIn->SetBranchAddress("weight_reint_all",                              &weight_reint_all);
            treeIn->SetBranchAddress("weight_All_UBGenie",                          &weight_All_UBGenie);
            treeIn->SetBranchAddress("weight_AxFFCCQEshape_UBGenie",      &weight_AxFFCCQEshape_UBGenie);
            treeIn->SetBranchAddress("weight_DecayAngMEC_UBGenie",          &weight_DecayAngMEC_UBGenie);
            treeIn->SetBranchAddress("weight_NormCCCOH_UBGenie",              &weight_NormCCCOH_UBGenie);
            treeIn->SetBranchAddress("weight_NormNCCOH_UBGenie",              &weight_NormNCCOH_UBGenie);
            treeIn->SetBranchAddress("weight_RPA_CCQE_UBGenie",                &weight_RPA_CCQE_UBGenie);
            treeIn->SetBranchAddress("weight_ThetaDelta2NRad_UBGenie",  &weight_ThetaDelta2NRad_UBGenie);
            treeIn->SetBranchAddress("weight_Theta_Delta2Npi_UBGenie",  &weight_Theta_Delta2Npi_UBGenie);
            treeIn->SetBranchAddress("weight_VecFFCCQEshape_UBGenie",    &weight_VecFFCCQEshape_UBGenie);
            treeIn->SetBranchAddress("weight_XSecShape_CCMEC_UBGenie",  &weight_XSecShape_CCMEC_UBGenie);
            treeIn->SetBranchAddress("weight_xsr_scc_Fa3_SCC",                  &weight_xsr_scc_Fa3_SCC);
            treeIn->SetBranchAddress("weight_xsr_scc_Fv3_SCC",                  &weight_xsr_scc_Fv3_SCC);
        }

        // Declare temporary histograms
        CategoryHistogram hTempTrackScoreMuon = hTrackScoreMuon.CloneAndReset();
        CategoryHistogram hTempTrackDistMuon  = hTrackDistMuon.CloneAndReset();
        CategoryHistogram hTempTrackLenMuon   = hTrackLenMuon.CloneAndReset();
        CategoryHistogram hTempTrackPIDMuon   = hTrackPIDMuon.CloneAndReset();
        CategoryHistogram hTempTrackPIDProton = hTrackPIDProton.CloneAndReset();

        // Get number of events
        Long64_t nEvents = treeIn->GetEntries();
        std::cout << "    Total number of events: " << nEvents << std::endl;

        // Progress counter
        Int_t reportEvery = nEvents / 10;
        if (reportEvery == 0) reportEvery = 1;

        // Declare weight map
        std::map<TString, Universe> weight_map;

        // Main event loop
        for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {

            //if (iEvent >= 1000) break;

            // Print progress
            if (iEvent % reportEvery == 0) {
                std::cout << "        Processing event " << iEvent << " (" 
                          << (100.0 * iEvent / nEvents) << "%)" << std::endl;
            }

            // Load the current event
            treeIn->GetEntry(iEvent);

            // Fill weight map for MC only
            if (file.is_reweightable_mc) {
                weight_map.clear();
                weight_map["flux_all"]                = {*weight_flux_all, tuned_cv_weight * spline_weight, true};
                weight_map["reint_all"]               = {*weight_reint_all, tuned_cv_weight * spline_weight, true};
                weight_map["All_UBGenie"]             = {*weight_All_UBGenie, spline_weight, true};
                weight_map["AxFFCCQEshape_UBGenie"]   = {*weight_AxFFCCQEshape_UBGenie, spline_weight, false};
                weight_map["DecayAngMEC_UBGenie"]     = {*weight_DecayAngMEC_UBGenie, spline_weight, false};
                weight_map["NormCCCOH_UBGenie"]       = {*weight_NormCCCOH_UBGenie, spline_weight, false};
                weight_map["NormNCCOH_UBGenie"]       = {*weight_NormNCCOH_UBGenie, spline_weight, false};
                weight_map["RPA_CCQE_UBGenie"]        = {*weight_RPA_CCQE_UBGenie, spline_weight, true};
                weight_map["ThetaDelta2NRad_UBGenie"] = {*weight_ThetaDelta2NRad_UBGenie, spline_weight, false};
                weight_map["Theta_Delta2Npi_UBGenie"] = {*weight_Theta_Delta2Npi_UBGenie, spline_weight, false};
                weight_map["VecFFCCQEshape_UBGenie"]  = {*weight_VecFFCCQEshape_UBGenie, spline_weight, false};
                weight_map["XSecShape_CCMEC_UBGenie"] = {*weight_XSecShape_CCMEC_UBGenie, spline_weight, false};
                weight_map["xsr_scc_Fa3_SCC"]         = {*weight_xsr_scc_Fa3_SCC, tuned_cv_weight * spline_weight, true};
                weight_map["xsr_scc_Fv3_SCC"]         = {*weight_xsr_scc_Fv3_SCC, tuned_cv_weight * spline_weight, true};
            }

            // Weight MC events using the MicroBooNE CV tune
            Float_t weight = 1.0;
            if (is_mc) {
                weight = tuned_cv_weight * spline_weight;
            }
            
            // First selection cut
            if (!(reco_vertex_in_FV && pfp_starts_in_PCV && topo_cut_passed)) continue;

            Int_t nParticles = trk_score_v->size();

            for (Long64_t jParticle = 0; jParticle < nParticles; jParticle++) {

                if (pfp_generation_v->at(jParticle) != 2) continue; // only direct neutrino daughters

                Int_t category = -1;
                if (is_mc) {
                    category = pdg_to_category(backtracked_pdg->at(jParticle));
                } else if (strcmp(file.type, "extBNB") == 0) {
                    category = EXT_category;
                } else {
                    category = Data_category;
                }

                Float_t trk_score = trk_score_v->at(jParticle);
                FillHistograms(file, hTempTrackScoreMuon, weight_map, trk_score, category, weight);

                // Select particles passing the track score cut
                if (!(trk_score >= 0.8)) continue;

                Float_t distance = trk_distance_v->at(jParticle);
                FillHistograms(file, hTempTrackDistMuon, weight_map, distance, category, weight);

                // Select particles with distance to vertex less than cut
                if (!(distance <= 4.0)) continue;

                Float_t length = trk_len_v->at(jParticle);
                FillHistograms(file, hTempTrackLenMuon, weight_map, length, category, weight);

                // Select particles with length larger than cut
                if (!(length >= 10.0)) continue;

                Float_t llr_pid_score = trk_llr_pid_score_v->at(jParticle);
                FillHistograms(file, hTempTrackPIDMuon, weight_map, llr_pid_score, category, weight);
                
            }

            // Second selection cut
            if (!(numu_cc && muon_passed_mom_cuts && no_reco_showers && has_p_candidate && protons_contained)) continue;

            for (Long64_t jParticle = 0; jParticle < nParticles; jParticle++) {

                if (pfp_generation_v->at(jParticle) != 2) continue; // only direct neutrino daughters

                Int_t category = -1;
                if (is_mc) {
                    category = pdg_to_category(backtracked_pdg->at(jParticle));
                } else if (strcmp(file.type, "extBNB") == 0) {
                    category = EXT_category;
                } else {
                    category = Data_category;
                }

                Float_t llr_pid_score = trk_llr_pid_score_v->at(jParticle);
                FillHistograms(file, hTempTrackPIDProton, weight_map, llr_pid_score, category, weight);

            }

        }

        // Scale the histograms
        hTempTrackScoreMuon.Scale(ScaleFactor);
        hTempTrackDistMuon.Scale(ScaleFactor);
        hTempTrackLenMuon.Scale(ScaleFactor);
        hTempTrackPIDMuon.Scale(ScaleFactor);
        hTempTrackPIDProton.Scale(ScaleFactor);

        // Add to the outer histograms
        hTrackScoreMuon.Add(hTempTrackScoreMuon);
        hTrackDistMuon.Add(hTempTrackDistMuon);
        hTrackLenMuon.Add(hTempTrackLenMuon);
        hTrackPIDMuon.Add(hTempTrackPIDMuon);
        hTrackPIDProton.Add(hTempTrackPIDProton);

        // Close input file
        if (fileIn && fileIn->IsOpen()) fileIn->Close();

    }

    fileOut->cd();
    hTrackScoreMuon.Finalize(fileOut, "TrackScoreMuon", detvar_map, false, false);
    hTrackDistMuon.Finalize(fileOut, "TrackDistMuon", detvar_map, false, false);
    hTrackLenMuon.Finalize(fileOut, "TrackLenMuon", detvar_map, false, false);
    hTrackPIDMuon.Finalize(fileOut, "TrackPIDMuon", detvar_map, false, false);
    hTrackPIDProton.Finalize(fileOut, "TrackPIDProton", detvar_map, false, false);
    
    // Close output file
    if (fileOut && fileOut->IsOpen()) fileOut->Close();

}