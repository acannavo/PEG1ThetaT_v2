// ====================================================================
// AnalysisAN.C
// ====================================================================
// Systematic study of analyzing power 
//
// Purpose: Generate events for multiple detector configurations to study
//          how detector misalignment affects analyzing power measurements
//			and to study the dilution  
//
// Study Design:
//   - Fixed right detector at R = 16.2°
//   - Scan left detector from L = 14.0° to 18.0° in 0.2° steps
//   - Generate both spin-up and spin-down events
//   - Support both elastic and inelastic channels
//
// Output: ROOT files with detector angle encoded in filename
//   Elastic:   pC_Elas_<E>MeV_MT_P<pol>_<spin>_<theta>p<decimal>.root
//   Inelastic: pC_Inel443_<E>MeV_MT_P<pol>_<spin>_<theta>p<decimal>.root
//
// Usage:
//   root -l
//   .L sigma_pC_Elas_EventGenerator.C
//   .L AnalysisAN.C
//   RunSystematicStudy()              // Run with default parameters
//   RunSystematicStudy(200, 50000)   // Custom energy and events
//
// Requirements:
//   - DetectorConfig.C loaded (for SetDetectorConfig)
//   - EventGenerator.C loaded (for generation functions)
//   - sigma_pC_Elas_EventGenerator.C loaded (main program)
//
// Author: [Your Name]
// Date: January 2026
// ====================================================================

#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

#include "TString.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "AnalyzingPowerUtils.h"
#include "ElasticScattering.h"
#include "InelasticScattering.h"
#include "Kinematics.h"
#include "DetectorConfig.h"
#include "EventGenerator.h"
#include "AnalyzingPowerUtils.h"


using namespace std;

// ====================================================================
// STUDY CONFIGURATION STRUCTURE
// ====================================================================
struct StudyConfig {
    Double_t energy;           // Beam energy [MeV]
    Int_t events_per_config;   // Events per detector configuration per spin
    Double_t polarization;     // Beam polarization (0.0 to 1.0)
    Bool_t run_elastic;        // Generate elastic events
    Bool_t run_inelastic;      // Generate inelastic events
    vector<Double_t> detector_angles;  // Detector positions to scan [degrees]
    
    // Constructor with defaults
    StudyConfig() {
        energy = 200.0;              // 200 MeV
        events_per_config = 100000;  // 100k events per config per spin
        polarization = 0.80;         // 80% polarization
        run_elastic = true;
        run_inelastic = true;
        
        // Default scan: 14.0° to 18.0° in 0.2° steps (21 positions)
        for (Double_t angle = 14.0; angle <= 18.0; angle += 0.2) {
            detector_angles.push_back(angle);
        }
    }
};

// ====================================================================
// PrintStudyHeader - Display study parameters
// ====================================================================
void PrintStudyHeader(const StudyConfig& config) {
    cout << "\n========================================" << endl;
    cout << "  ANALYZING POWER SYSTEMATIC STUDY" << endl;
    cout << "========================================" << endl;
    cout << "\nStudy Parameters:" << endl;
    cout << "  Beam energy:       " << config.energy << " MeV" << endl;
    cout << "  Polarization:      " << config.polarization * 100 << "%" << endl;
    cout << "  Events/config:     " << config.events_per_config << " per spin state" << endl;
    cout << "  Total events:      " << config.events_per_config * 2 << " per configuration" << endl;
    cout << "\nChannels:" << endl;
    cout << "  Elastic:           " << (config.run_elastic ? "YES" : "NO") << endl;
    cout << "  Inelastic (4.43):  " << (config.run_inelastic ? "YES" : "NO") << endl;
    cout << "\nDetector Scan:" << endl;
    cout << "  Number of positions: " << config.detector_angles.size() << endl;
    cout << "  Range:             " << config.detector_angles.front() << "° to " 
         << config.detector_angles.back() << "°" << endl;
    cout << "  Step size:         " << (config.detector_angles.size() > 1 ? 
         config.detector_angles[1] - config.detector_angles[0] : 0.0) << "°" << endl;
    
    // Calculate total datasets
    Int_t datasets_per_config = 0;
    if (config.run_elastic) datasets_per_config += 2;    // SpinUp + SpinDown
    if (config.run_inelastic) datasets_per_config += 2;  // SpinUp + SpinDown
    Int_t total_datasets = config.detector_angles.size() * datasets_per_config;
    
    cout << "\nExpected Output:" << endl;
    cout << "  Configurations:    " << config.detector_angles.size() << endl;
    cout << "  Datasets/config:   " << datasets_per_config << endl;
    cout << "  Total datasets:    " << total_datasets << endl;
    cout << "  Total events:      " << total_datasets * config.events_per_config << endl;
    
    // Estimate time (rough: ~10k events/sec)
    Double_t estimated_seconds = (total_datasets * config.events_per_config) / 10000.0;
    Int_t hours = (Int_t)(estimated_seconds / 3600);
    Int_t minutes = (Int_t)((estimated_seconds - hours*3600) / 60);
    cout << "\nEstimated time:    ~" << hours << "h " << minutes << "m" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// PrintProgress - Display progress during generation
// ====================================================================
void PrintProgress(Int_t current, Int_t total, Double_t angle, 
                   const char* channel, const char* spin) {
    cout << "\n>>> Configuration " << current << "/" << total << " <<<" << endl;
    cout << "  Detector angle: " << angle << "°" << endl;
    cout << "  Channel:        " << channel << endl;
    cout << "  Spin state:     " << spin << endl;
}

// ====================================================================
// RunSystematicStudy - Main driver function
// ====================================================================
// Generate events for systematic study of detector positioning effects
//
// Parameters:
//   energy          - Beam energy [MeV] (default: 200)
//   events_per_run  - Events per configuration per spin (default: 100000)
//   polarization    - Beam polarization 0-1 (default: 0.80)
//   run_elastic     - Generate elastic events (default: true)
//   run_inelastic   - Generate inelastic events (default: true)
//
// Example usage:
//   RunSystematicStudy()                    // Use all defaults
//   RunSystematicStudy(200, 50000)         // Custom events
//   RunSystematicStudy(200, 100000, 0.70)  // Custom polarization
// ====================================================================
void RunSystematicStudy(Double_t energy = 200.0,
                        Int_t events_per_run = 100000,
                        Double_t polarization = 0.80,
                        Bool_t run_elastic = true,
                        Bool_t run_inelastic = true,
                        vector<Double_t>* custom_angles = nullptr)
{
    // ================================================================
    // 1. Setup configuration
    // ================================================================
    StudyConfig config;
    config.energy = energy;
    config.events_per_config = events_per_run;
    config.polarization = polarization;
    config.run_elastic = run_elastic;
    config.run_inelastic = run_inelastic;
    
    // Use custom angles if provided, otherwise use default range
    if (custom_angles != nullptr && custom_angles->size() > 0) {
        config.detector_angles = *custom_angles;
    }
    
    // Print study information
    PrintStudyHeader(config);
    
    // Start timer
    TStopwatch timer;
    timer.Start();
    
    // ================================================================
    // 2. Main generation loop
    // ================================================================
    Int_t total_configs = config.detector_angles.size();
    Int_t datasets_per_config = 0;
    if (config.run_elastic) datasets_per_config += 2;
    if (config.run_inelastic) datasets_per_config += 2;
    Int_t total_datasets = total_configs * datasets_per_config;
    Int_t current_dataset = 0;
    
    cout << "Starting systematic scan...\n" << endl;
    
    // Loop over all detector angles
    for (Int_t i = 0; i < total_configs; i++) {
        Double_t theta_det = config.detector_angles[i];
        
        cout << "\n========================================" << endl;
        cout << "  Position " << (i+1) << "/" << total_configs 
             << " : Detector = " << theta_det << "°" << endl;
        cout << "========================================" << endl;
        
        // Set detector configuration
        SetDetectorConfig(theta_det);
        
        // ============================================================
        // ELASTIC CHANNEL
        // ============================================================
        if (config.run_elastic) {
            // Spin-Up
            current_dataset++;
            PrintProgress(current_dataset, total_datasets, theta_det, 
                         "Elastic", "UP");
            SingleRunMultithreadPolarized(config.energy, 
                                         config.events_per_config, 
                                         config.polarization, 
                                         +1);  // Spin up
            
            // Spin-Down
            current_dataset++;
            PrintProgress(current_dataset, total_datasets, theta_det, 
                         "Elastic", "DOWN");
            SingleRunMultithreadPolarized(config.energy, 
                                         config.events_per_config, 
                                         config.polarization, 
                                         -1);  // Spin down
        }
        
        // ============================================================
        // INELASTIC CHANNEL
        // ============================================================
        if (config.run_inelastic) {
            // Spin-Up
            current_dataset++;
            PrintProgress(current_dataset, total_datasets, theta_det, 
                         "Inelastic (4.43 MeV)", "UP");
            SingleRunMultithreadInelasticPolarized(config.energy, 
                                                  config.events_per_config, 
                                                  config.polarization, 
                                                  +1);  // Spin up
            
            // Spin-Down
            current_dataset++;
            PrintProgress(current_dataset, total_datasets, theta_det, 
                         "Inelastic (4.43 MeV)", "DOWN");
            SingleRunMultithreadInelasticPolarized(config.energy, 
                                                  config.events_per_config, 
                                                  config.polarization, 
                                                  -1);  // Spin down
        }
        
        cout << "\n>>> Configuration " << (i+1) << "/" << total_configs 
             << " COMPLETE <<<\n" << endl;
    }
    
    // ================================================================
    // 3. Summary
    // ================================================================
    timer.Stop();
    
    cout << "\n========================================" << endl;
    cout << "  SYSTEMATIC STUDY COMPLETE!" << endl;
    cout << "========================================" << endl;
    cout << "\nGeneration Statistics:" << endl;
    cout << "  Configurations:    " << total_configs << endl;
    cout << "  Total datasets:    " << total_datasets << endl;
    cout << "  Total events:      " << total_datasets * config.events_per_config << endl;
    cout << "  Elapsed time:      " << timer.RealTime() << " s" << endl;
    cout << "  Events/second:     " << (total_datasets * config.events_per_config) / timer.RealTime() << endl;
    
    cout << "\nOutput Files:" << endl;
    if (config.run_elastic) {
        cout << "  Elastic files:     " << total_configs * 2 << endl;
        cout << "    Format: pC_Elas_" << (Int_t)config.energy << "MeV_MT_P" 
             << (Int_t)(config.polarization*100) << "_SpinUp/Down_XXpX.root" << endl;
    }
    if (config.run_inelastic) {
        cout << "  Inelastic files:   " << total_configs * 2 << endl;
        cout << "    Format: pC_Inel443_" << (Int_t)config.energy << "MeV_MT_P" 
             << (Int_t)(config.polarization*100) << "_SpinUp/Down_XXpX.root" << endl;
    }
    
    cout << "\nNext Steps:" << endl;
    cout << "  1. Verify files were created (ls -lh *.root)" << endl;
    cout << "  2. Run analysis script to calculate analyzing power" << endl;
    cout << "  3. Plot A_N vs detector position" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// RunQuickTest - Generate only 3 configurations for testing
// ====================================================================
// Quick test with minimal configurations to verify everything works
//
// Generates only 3 detector positions: 14.0°, 16.2°, 18.0°
// Useful for testing before running full study
// ====================================================================
void RunQuickTest(Double_t energy = 200.0,
                  Int_t events_per_run = 100000,
                  Double_t polarization = 0.80)
{
    cout << "\n========================================" << endl;
    cout << "  QUICK TEST MODE" << endl;
    cout << "========================================" << endl;
    cout << "Running with only 3 detector positions:" << endl;
    cout << "  14.0°, 16.2°, 18.0°" << endl;
    cout << "  " << events_per_run << " events per dataset" << endl;
    cout << "========================================\n" << endl;
    
    // Create custom angle list with positions:
    vector<Double_t> test_angles;
    //test_angles.push_back(14.0);
	//test_angles.push_back(14.5);
	//test_angles.push_back(15.0);
	//test_angles.push_back(15.5);
	//test_angles.push_back(16.0);
    test_angles.push_back(16.2);
	//test_angles.push_back(17.0);
	//test_angles.push_back(17.5);
    //test_angles.push_back(18.0);
    
    // Run with custom angles
    RunSystematicStudy(energy, events_per_run, polarization, 
                      true, true, &test_angles);
}

// ====================================================================
// RunElasticOnly - Generate only elastic channel
// ====================================================================
void RunElasticOnly(Double_t energy = 200.0,
                    Int_t events_per_run = 100000,
                    Double_t polarization = 0.80)
{
    cout << "\n=== ELASTIC CHANNEL ONLY ===" << endl;
    RunSystematicStudy(energy, events_per_run, polarization, 
                      true,   // elastic
                      false); // no inelastic
}

// ====================================================================
// RunInelasticOnly - Generate only inelastic channel
// ====================================================================
void RunInelasticOnly(Double_t energy = 200.0,
                      Int_t events_per_run = 100000,
                      Double_t polarization = 0.80)
{
    cout << "\n=== INELASTIC CHANNEL ONLY ===" << endl;
    RunSystematicStudy(energy, events_per_run, polarization, 
                      false,  // no elastic
                      true);  // inelastic
}

// ====================================================================
// RunCustomRange - Generate with custom detector angle range
// ====================================================================
// Example: RunCustomRange(200, 100000, 0.80, 15.0, 17.0, 0.5)
//          Scans from 15.0° to 17.0° in 0.5° steps
// ====================================================================
void RunCustomRange(Double_t energy,
                    Int_t events_per_run,
                    Double_t polarization,
                    Double_t angle_min,
                    Double_t angle_max,
                    Double_t angle_step = 0.2)
{
    cout << "\n=== CUSTOM ANGLE RANGE ===" << endl;
    cout << "Range: " << angle_min << "° to " << angle_max 
         << "° in " << angle_step << "° steps" << endl;
    
    // Build custom angle range
    vector<Double_t> custom_angles;
    for (Double_t angle = angle_min; angle <= angle_max + 0.001; angle += angle_step) {
        custom_angles.push_back(angle);
    }
    
    cout << "Total angles: " << custom_angles.size() << endl;
    
    RunSystematicStudy(energy, events_per_run, polarization, 
                      true, true, &custom_angles);
}

// ====================================================================
// PrintUsageExamples - Show how to use this script
// ====================================================================
void PrintUsageExamples() {
    cout << "\n========================================" << endl;
    cout << "  AnalysisAN.C Usage Examples" << endl;
    cout << "========================================" << endl;
    cout << "\n1. Full systematic study (default settings):" << endl;
    cout << "   RunSystematicStudy()" << endl;
    cout << "   → 200 MeV, 100k events, 80% pol, 21 angles\n" << endl;
    
    cout << "2. Custom energy and events:" << endl;
    cout << "   RunSystematicStudy(200, 50000)" << endl;
    cout << "   → 200 MeV, 50k events per dataset\n" << endl;
    
    cout << "3. Custom polarization:" << endl;
    cout << "   RunSystematicStudy(200, 100000, 0.70)" << endl;
    cout << "   → 70% polarization\n" << endl;
    
    cout << "4. Quick test (only 3 positions):" << endl;
    cout << "   RunQuickTest()" << endl;
    cout << "   → Fast test with 14°, 16.2°, 18°\n" << endl;
    
    cout << "5. Elastic channel only:" << endl;
    cout << "   RunElasticOnly()" << endl;
    cout << "   → Skip inelastic channel\n" << endl;
    
    cout << "6. Inelastic channel only:" << endl;
    cout << "   RunInelasticOnly()" << endl;
    cout << "   → Skip elastic channel\n" << endl;
    
    cout << "7. Custom angle range:" << endl;
    cout << "   RunCustomRange(200, 100000, 0.80, 15.0, 17.0, 0.5)" << endl;
    cout << "   → Scan 15° to 17° in 0.5° steps\n" << endl;
    
    cout << "========================================\n" << endl;
}



// ====================================================================
// ====================================================================
// 		POLARIMETER Misalignemnt Study
// ====================================================================
// ====================================================================




// ====================================================================
// ANALYSIS CODE
// ====================================================================
// This version properly separates right and left detector events
// based on the azimuthal angle φ
// ====================================================================

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TLatex.h"
#include "PParticle.h"
#include <map>
#include <fstream>

// ====================================================================
// DATA INPUT CONFIGURATION
// ====================================================================
const TString DATA_INPUT_PATH = "Data/";  // Path to data files

// Helper function to construct input path
TString MakeInputPath(const char* filename) {
    return DATA_INPUT_PATH + TString(filename);
}

// ====================================================================
// STRUCTURES FOR ANALYSIS
// ====================================================================

/* // Store event counts for one detector configuration
struct DetectorCounts {
    Double_t theta_det;      // Detector angle [degrees]
    Int_t N_right_up;        // Right detector (φ≈0), spin-up
    Int_t N_right_down;      // Right detector (φ≈0), spin-down
    Int_t N_left_up;         // Left detector (φ≈π), spin-up
    Int_t N_left_down;       // Left detector (φ≈π), spin-down
    
    DetectorCounts() : theta_det(0), N_right_up(0), N_right_down(0), 
                       N_left_up(0), N_left_down(0) {}
}; */

/* // Store analyzing power result
struct ANResult {
    Double_t theta_L;        // Left detector angle [degrees]
    Double_t theta_R;        // Right detector angle [degrees]
    Double_t AN;             // Analyzing power
    Double_t AN_error;       // Statistical error
    Int_t N_R_up;            // Right detector, spin-up
    Int_t N_R_down;          // Right detector, spin-down
    Int_t N_L_up;            // Left detector, spin-up
    Int_t N_L_down;          // Left detector, spin-down
    
    ANResult() : theta_L(0), theta_R(0), AN(0), AN_error(0),
                 N_R_up(0), N_R_down(0), N_L_up(0), N_L_down(0) {}
}; */

// ====================================================================
// ExtractDetectorAngle - Parse angle from filename
// ====================================================================
Double_t ExtractDetectorAngle(TString filename) {
    Ssiz_t last_underscore = filename.Last('_');
    Ssiz_t dot_root = filename.Last('.');
    
    if (last_underscore == kNPOS || dot_root == kNPOS) {
        cout << "Warning: Could not parse angle from: " << filename << endl;
        return -999.0;
    }
    
    TString angle_str = filename(last_underscore + 1, dot_root - last_underscore - 1);
    angle_str.ReplaceAll("p", ".");
    
    return angle_str.Atof();
}

/* // ====================================================================
// CountEventsInFileByPhi - Count events separated by φ angle
// ====================================================================
// Separates events into right detector (φ ≈ 0) and left detector (φ ≈ π)
// Uses detector acceptance windows to determine which events go where
//
// Parameters:
//   filename  - ROOT file to analyze
//   phi_window - Azimuthal acceptance window [radians] (default: 0.017)
//   is_spinup - True for spin-up, false for spin-down
//   N_right   - Output: events with φ ≈ 0 (right detector)
//   N_left    - Output: events with φ ≈ π (left detector)
// ====================================================================
void CountEventsInFileByPhi(TString filename, 
                            Double_t phi_window,
                            Bool_t is_spinup,
                            Int_t& N_right, 
                            Int_t& N_left) {
    N_right = 0;
    N_left = 0;
    
    TString fullpath = MakeInputPath(filename.Data());
    
    TFile* file = TFile::Open(fullpath.Data());
    if (!file || file->IsZombie()) {
        cout << "ERROR: Cannot open " << fullpath << endl;
        return;
    }
    
    TTree* tree = (TTree*)file->Get("data");
    if (!tree) {
        cout << "ERROR: No 'data' tree in " << fullpath << endl;
        file->Close();
        return;
    }
    
    // Setup branch reading
    TClonesArray* particles = new TClonesArray("PParticle");
    tree->SetBranchAddress("Particles", &particles);
    
    Int_t nentries = tree->GetEntries();
    
    // Define φ acceptance windows
    const Double_t phi_min_right = -phi_window;
    const Double_t phi_max_right = +phi_window;
    const Double_t phi_min_left = TMath::Pi() - phi_window;
    const Double_t phi_max_left = TMath::Pi() + phi_window;
    
    // Loop over all events
    for (Int_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        
        // Get proton (first particle)
        if (particles->GetEntries() < 1) continue;
        PParticle* proton = (PParticle*)particles->At(0);
        
        // Get φ angle
        Double_t phi = proton->Phi();
        
        // Check which detector this event hit
        Bool_t in_right = (phi >= phi_min_right && phi <= phi_max_right);
        Bool_t in_left = (phi >= phi_min_left && phi <= phi_max_left);
        
        if (in_right) {
            N_right++;
        } else if (in_left) {
            N_left++;
        }
        // Events outside both windows are ignored (shouldn't happen)
    }
    
    file->Close();
    delete file;
    delete particles;
} */

// ====================================================================
// FindDataFiles - Scan directory and count events by detector
// ====================================================================
map<Double_t, DetectorCounts> FindDataFiles(const char* channel, 
                                            Int_t energy, 
                                            Int_t polarization,
                                            Double_t phi_window = 0.017) {
    map<Double_t, DetectorCounts> data_map;
    
    TString pattern = Form("%s_%dMeV_MT_P%d", channel, energy, polarization);
    
    cout << "\n=== Scanning for files: " << pattern << "_*.root ===" << endl;
    cout << "Search path: " << DATA_INPUT_PATH << endl;
    cout << "Phi window: ±" << phi_window << " rad (±" 
         << phi_window*TMath::RadToDeg() << "°)" << endl;
    
    void* dir = gSystem->OpenDirectory(DATA_INPUT_PATH.Data());
    if (!dir) {
        cout << "ERROR: Cannot open directory: " << DATA_INPUT_PATH << endl;
        return data_map;
    }
    
    const char* entry;
    Int_t nfiles_found = 0;
    
    while ((entry = gSystem->GetDirEntry(dir))) {
        TString filename(entry);
        
        if (!filename.BeginsWith(pattern) || !filename.EndsWith(".root")) {
            continue;
        }
        
        Double_t theta_det = ExtractDetectorAngle(filename);
        Bool_t is_spinup = filename.Contains("SpinUp");
        
        if (theta_det < 0) continue;
        
        // Count events by φ (right vs left)
        Int_t N_right, N_left;
        CountEventsInFileByPhi(DATA_INPUT_PATH + filename, phi_window, is_spinup, N_right, N_left);
        
        // Initialize map entry if needed
        if (data_map.find(theta_det) == data_map.end()) {
            data_map[theta_det] = DetectorCounts();
            data_map[theta_det].theta_det = theta_det;
        }
        
        // Store counts
        if (is_spinup) {
            data_map[theta_det].N_right_up += N_right;
            data_map[theta_det].N_left_up += N_left;
            cout << "  ✓ " << filename << " : " 
                 << N_right << " right, " << N_left << " left (spin-up)" << endl;
        } else {
            data_map[theta_det].N_right_down += N_right;
            data_map[theta_det].N_left_down += N_left;
            cout << "  ✓ " << filename << " : " 
                 << N_right << " right, " << N_left << " left (spin-down)" << endl;
        }
        
        nfiles_found++;
    }
    
    gSystem->FreeDirectory(dir);
    
    cout << "\nFound " << nfiles_found << " files" << endl;
    cout << "Configurations: " << data_map.size() << endl;
    
    return data_map;
}


/* 
// ====================================================================
// CalculateAsymmetry - Compute spin asymmetry from correlation
// ====================================================================
// Formula: ε = [√(N_R_up × N_L_down) - √(N_L_up × N_R_down)] /
//              [√(N_R_up × N_L_down) + √(N_L_up × N_R_down)]
//
// This is the ASYMMETRY (ε), not analyzing power (A_N)
// To get A_N, divide by beam polarization: A_N = ε / P
//
// Note: Result stored in result.AN for backward compatibility
// ====================================================================
ANResult CalculateAsymmetry(Double_t theta_R, Double_t theta_L,
                            Int_t N_R_up, Int_t N_R_down,
                            Int_t N_L_up, Int_t N_L_down) {
    ANResult result;
    result.theta_R = theta_R;
    result.theta_L = theta_L;
    result.N_R_up = N_R_up;
    result.N_R_down = N_R_down;
    result.N_L_up = N_L_up;
    result.N_L_down = N_L_down;
    
    if (N_R_up == 0 || N_R_down == 0 || N_L_up == 0 || N_L_down == 0) {
        cout << "WARNING: Zero counts detected!" << endl;
        result.AN = 0;
        result.AN_error = 999.0;
        return result;
    }
    
    // Cast to Double_t BEFORE multiplication to prevent integer overflow
    Double_t term1 = TMath::Sqrt((Double_t)N_R_up * (Double_t)N_L_down);
    Double_t term2 = TMath::Sqrt((Double_t)N_L_up * (Double_t)N_R_down);
    
    Double_t numerator = term1 - term2;
    Double_t denominator = term1 + term2;

    if (denominator == 0 || std::isnan(numerator) || std::isnan(denominator)) {
        cout << "WARNING: Invalid math at theta_L = " << theta_L << endl;
        result.AN = 0;
        result.AN_error = 999;
        return result;
    }
    
    // Calculate asymmetry
    result.AN = numerator / denominator;
    
    // Correct error propagation for asymmetry measurement
    // For ε = (√(a·b) - √(c·d)) / (√(a·b) + √(c·d))
    // where a=N_R_up, b=N_L_down, c=N_L_up, d=N_R_down
    
    // Relative errors on sqrt terms (Poisson statistics)
    Double_t rel_err_term1 = 0.5 * TMath::Sqrt(1.0/N_R_up + 1.0/N_L_down);
    Double_t rel_err_term2 = 0.5 * TMath::Sqrt(1.0/N_L_up + 1.0/N_R_down);

    // Absolute errors on terms
    Double_t err_term1 = term1 * rel_err_term1;
    Double_t err_term2 = term2 * rel_err_term2;

    // Error on numerator (uncorrelated terms)
    Double_t err_numerator = TMath::Sqrt(err_term1*err_term1 + err_term2*err_term2);
    
    // Error on denominator (same terms)
    Double_t err_denominator = err_numerator;

    // Full error propagation for quotient f = (a-b)/(a+b)
    Double_t rel_err = TMath::Sqrt(TMath::Power(err_numerator/numerator, 2) + 
                                    TMath::Power(err_denominator/denominator, 2));
    result.AN_error = TMath::Abs(result.AN) * rel_err;
    
    return result;
}
 */

// ====================================================================
// AnalyzeDetectorSystematics - Main analysis function
// ====================================================================
void AnalyzeDetectorSystematics(const char* channel = "pC_Elas",
                                Int_t energy = 200,
                                Int_t polarization = 80,
                                Double_t theta_R_fixed = 16.2,
                                const char* output_name = "AN_systematic",
                                Double_t phi_window = 0.017) {
    
    cout << "\n========================================" << endl;
    cout << "  ANALYZING POWER SYSTEMATIC ANALYSIS" << endl;
    cout << "========================================" << endl;
    cout << "Channel:           " << channel << endl;
    cout << "Energy:            " << energy << " MeV" << endl;
    cout << "Polarization:      " << polarization << "%" << endl;
    cout << "Right detector:    " << theta_R_fixed << "° (fixed)" << endl;
    cout << "Data path:         " << DATA_INPUT_PATH << endl;
    cout << "========================================" << endl;
    
    // Find and load data files
    map<Double_t, DetectorCounts> data_map = FindDataFiles(channel, energy, 
                                                           polarization, phi_window);
    
    if (data_map.size() == 0) {
        cout << "\nERROR: No data files found!" << endl;
        return;
    }
    
    // Check for right detector data
    if (data_map.find(theta_R_fixed) == data_map.end()) {
        cout << "\nERROR: No data for right detector at " << theta_R_fixed << "°!" << endl;
        cout << "Available angles: ";
        for (auto& pair : data_map) {
            cout << pair.first << "° ";
        }
        cout << endl;
        return;
    }
    
    // Get right detector counts (from RIGHT detector at θ_R)
    DetectorCounts data_at_R = data_map[theta_R_fixed];
    Int_t N_R_up = data_at_R.N_right_up;      // Right detector, spin up
    Int_t N_R_down = data_at_R.N_right_down;  // Right detector, spin down
    
    cout << "\nRight detector (R=" << theta_R_fixed << "°):" << endl;
    cout << "  N_right_up   = " << N_R_up << endl;
    cout << "  N_right_down = " << N_R_down << endl;
    
    // Calculate analyzing power for each left detector position
    cout << "\n=== Computing Asymmetry and Analyzing Power ===" << endl;
    cout << Form("%-8s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s", 
                 "L [deg]", "N_L_up", "N_L_down", "N_R_up", "N_R_down", 
                 "Asym(ε)", "δε", "A_N", "δA_N") << endl;
    cout << "-----------------------------------------------------------------------------------------------------------" << endl;
    
    vector<ANResult> results;
    Double_t P = polarization / 100.0;  // Convert to fraction for A_N calculation
    
    for (auto& pair : data_map) {
        Double_t theta_L = pair.first;
        DetectorCounts data_at_L = pair.second;
        
        // Get left detector counts (from LEFT detector at θ_L)
        Int_t N_L_up = data_at_L.N_left_up;      // Left detector, spin up
        Int_t N_L_down = data_at_L.N_left_down;  // Left detector, spin down
        
        // Calculate asymmetry (stored in result.AN for backward compatibility)
        //ANResult result = CalculateAsymmetry(theta_R_fixed, theta_L,
        //                                     N_R_up, N_R_down,
        //                                    N_L_up, N_L_down);
        //ANResult asymmetry = CalculateAsymmetry(theta_R_fixed, theta_L,
        //                                N_R_up, N_R_down,
        //                               N_L_up, N_L_down);

		//results.push_back(result);
        
        // Calculate analyzing power from asymmetry
        //Double_t analyzing_power, analyzing_power_err;
		//CalculateAnalyzingPower(asymmetry, P, analyzing_power, analyzing_power_err);
		//Double_t analyzing_power = result.AN / P;
        //Double_t analyzing_power_err = result.AN_error / P;
        
		
		
		// Calculate asymmetry (stored in result.AN for backward compatibility)
		ANResult asymmetry = CalculateAsymmetry(theta_R_fixed, theta_L,
												N_R_up, N_R_down,
												N_L_up, N_L_down);

		// Calculate analyzing power from asymmetry
		Double_t P = polarization / 100.0;
		Double_t analyzing_power = AsymmetryToAnalyzingPower(asymmetry.AN, P);
		Double_t analyzing_power_err = AsymmetryErrorToAnalyzingPowerError(asymmetry.AN_error, P);

		// Store in results vector
		results.push_back(asymmetry);  // ← This is the issue - you're trying to use 'result' but it's called 'asymmetry'

		cout << Form("%7.1f  | %10d | %10d | %10d | %10d | %+9.5f | %10.5f | %+9.5f | %10.5f",
					 theta_L, N_L_up, N_L_down, N_R_up, N_R_down, 
					 asymmetry.AN, asymmetry.AN_error,              // Asymmetry
					 analyzing_power, analyzing_power_err)          // A_N
					 << endl;

    }
    
    if (results.size() == 0) {
        cout << "\nERROR: No valid left detector positions found!" << endl;
        return;
    }
    
    // Create output arrays for plotting
    Int_t npoints = results.size();
    vector<Double_t> theta_L_vec, asymmetry_vec, asymmetry_err_vec,
                     AN_vec, AN_err_vec, delta_theta_vec;
    
    for (const auto& result : results) {
        theta_L_vec.push_back(result.theta_L);
        
        // Asymmetry (what the formula directly gives)
        asymmetry_vec.push_back(result.AN);
        asymmetry_err_vec.push_back(result.AN_error);
        
        // Analyzing power (asymmetry / polarization)
        AN_vec.push_back(result.AN / P);
        AN_err_vec.push_back(result.AN_error / P);
        
        delta_theta_vec.push_back(result.theta_L - theta_R_fixed);
    }
    
    // Create plots
    cout << "\n=== Creating Plots ===" << endl;
    
    TCanvas* c1 = new TCanvas("c1", "Asymmetry and Analyzing Power Systematics", 1600, 1000);
    c1->Divide(2, 2);  // 2x2 grid
    
    // ================================================================
    // TOP ROW: ASYMMETRY PLOTS
    // ================================================================
    
    // Plot 1: Asymmetry vs Left Detector Position
    c1->cd(1);
    gPad->SetGrid();
    
    TGraphErrors* gr_asym_1 = new TGraphErrors(npoints, &theta_L_vec[0], 
                                                &asymmetry_vec[0], 
                                                0, &asymmetry_err_vec[0]);
    gr_asym_1->SetTitle("Spin Asymmetry vs Left Detector Position;#theta_{L} [deg];#varepsilon");
    gr_asym_1->SetMarkerStyle(20);
    gr_asym_1->SetMarkerSize(1.2);
    gr_asym_1->SetMarkerColor(kBlue);
    gr_asym_1->SetLineColor(kBlue);
    gr_asym_1->Draw("APE");
    
    TLine* line1 = new TLine(theta_R_fixed, gr_asym_1->GetYaxis()->GetXmin(),
                             theta_R_fixed, gr_asym_1->GetYaxis()->GetXmax());
    line1->SetLineColor(kRed);
    line1->SetLineStyle(2);
    line1->SetLineWidth(2);
    line1->Draw("same");
    
    TLegend* leg1 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg1->AddEntry(gr_asym_1, "Asymmetry #varepsilon", "pe");
    leg1->AddEntry(line1, Form("Right det. (R=%.1f#circ)", theta_R_fixed), "l");
    leg1->Draw();
    
    // Plot 2: Asymmetry vs Detector Separation
    c1->cd(2);
    gPad->SetGrid();
    
    TGraphErrors* gr_asym_2 = new TGraphErrors(npoints, &delta_theta_vec[0], 
                                                &asymmetry_vec[0],
                                                0, &asymmetry_err_vec[0]);
    gr_asym_2->SetTitle("Spin Asymmetry vs Detector Separation;#theta_{L} - #theta_{R} [deg];#varepsilon");
    gr_asym_2->SetMarkerStyle(20);
    gr_asym_2->SetMarkerSize(1.2);
    gr_asym_2->SetMarkerColor(kBlue);
    gr_asym_2->SetLineColor(kBlue);
    gr_asym_2->Draw("APE");
    
    TLine* line2 = new TLine(0, gr_asym_2->GetYaxis()->GetXmin(),
                             0, gr_asym_2->GetYaxis()->GetXmax());
    line2->SetLineColor(kRed);
    line2->SetLineStyle(2);
    line2->SetLineWidth(2);
    line2->Draw("same");
    
    TLegend* leg2 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg2->AddEntry(gr_asym_2, "Asymmetry #varepsilon", "pe");
    leg2->AddEntry(line2, "#theta_{L} = #theta_{R}", "l");
    leg2->Draw();
    
    // ================================================================
    // BOTTOM ROW: ANALYZING POWER PLOTS
    // ================================================================
    
    // Plot 3: Analyzing Power vs Left Detector Position
    c1->cd(3);
    gPad->SetGrid();
    
    TGraphErrors* gr_AN_1 = new TGraphErrors(npoints, &theta_L_vec[0], 
                                              &AN_vec[0], 
                                              0, &AN_err_vec[0]);
    gr_AN_1->SetTitle(Form("Analyzing Power vs Left Detector Position (P=%d%%);#theta_{L} [deg];A_{N}", polarization));
    gr_AN_1->SetMarkerStyle(21);
    gr_AN_1->SetMarkerSize(1.2);
    gr_AN_1->SetMarkerColor(kGreen+2);
    gr_AN_1->SetLineColor(kGreen+2);
    gr_AN_1->Draw("APE");
    
    TLine* line3 = new TLine(theta_R_fixed, gr_AN_1->GetYaxis()->GetXmin(),
                             theta_R_fixed, gr_AN_1->GetYaxis()->GetXmax());
    line3->SetLineColor(kRed);
    line3->SetLineStyle(2);
    line3->SetLineWidth(2);
    line3->Draw("same");
    
    TLegend* leg3 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg3->AddEntry(gr_AN_1, "A_{N} = #varepsilon/P", "pe");
    leg3->AddEntry(line3, Form("Right det. (R=%.1f#circ)", theta_R_fixed), "l");
    leg3->Draw();
    
    // Plot 4: Analyzing Power vs Detector Separation
    c1->cd(4);
    gPad->SetGrid();
    
    TGraphErrors* gr_AN_2 = new TGraphErrors(npoints, &delta_theta_vec[0], 
                                              &AN_vec[0],
                                              0, &AN_err_vec[0]);
    gr_AN_2->SetTitle(Form("Analyzing Power vs Detector Separation (P=%d%%);#theta_{L} - #theta_{R} [deg];A_{N}", polarization));
    gr_AN_2->SetMarkerStyle(21);
    gr_AN_2->SetMarkerSize(1.2);
    gr_AN_2->SetMarkerColor(kGreen+2);
    gr_AN_2->SetLineColor(kGreen+2);
    gr_AN_2->Draw("APE");
    
    TLine* line4 = new TLine(0, gr_AN_2->GetYaxis()->GetXmin(),
                             0, gr_AN_2->GetYaxis()->GetXmax());
    line4->SetLineColor(kRed);
    line4->SetLineStyle(2);
    line4->SetLineWidth(2);
    line4->Draw("same");
    
    TLegend* leg4 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg4->AddEntry(gr_AN_2, "A_{N} = #varepsilon/P", "pe");
    leg4->AddEntry(line4, "#theta_{L} = #theta_{R}", "l");
    leg4->Draw();
    
	c1->Draw();
	c1->Update();
	
    TString pdf_name = Form("%s_%s_%dMeV.pdf", output_name, channel, energy);
    c1->SaveAs(pdf_name.Data());
    cout << "Plot saved: " << pdf_name << endl;
    
    // Save results to text file
    TString txt_name = Form("%s_%s_%dMeV.txt", output_name, channel, energy);
    ofstream outfile(txt_name.Data());
    
    outfile << "========================================" << endl;
    outfile << "ANALYZING POWER SYSTEMATIC STUDY RESULTS" << endl;
    outfile << "========================================" << endl;
    outfile << "Channel:        " << channel << endl;
    outfile << "Energy:         " << energy << " MeV" << endl;
    outfile << "Polarization:   " << polarization << "%" << endl;
    outfile << "Right detector: " << theta_R_fixed << " degrees" << endl;
    outfile << "========================================" << endl;
    outfile << endl;
    
    outfile << Form("%-10s %-10s %-12s %-10s %-10s %-10s %-10s %-12s %-12s %-12s %-12s %-12s",
                    "theta_L", "theta_R", "L-R", "N_L_up", "N_L_down", 
                    "N_R_up", "N_R_down", "Asymmetry", "delta_asym", "A_N", "delta_AN", "Signif") << endl;
    outfile << "----------------------------------------------------------------------------------------------------------------------" << endl;
    
    for (const auto& result : results) {
        Double_t delta_theta = result.theta_L - result.theta_R;
        Double_t analyzing_power = result.AN / P;
        Double_t analyzing_power_err = result.AN_error / P;
        Double_t significance = TMath::Abs(analyzing_power) / analyzing_power_err;
        
        outfile << Form("%10.1f %10.1f %12.2f %10d %10d %10d %10d %12.6f %12.6f %12.6f %12.6f %12.2f",
                        result.theta_L, result.theta_R, delta_theta,
                        result.N_L_up, result.N_L_down,
                        result.N_R_up, result.N_R_down,
                        result.AN, result.AN_error,              // Asymmetry
                        analyzing_power, analyzing_power_err,    // A_N
                        significance) << endl;
    }
    
    outfile.close();
    cout << "Results saved: " << txt_name << endl;
    
    // ================================================================
    // NEW PLOT: Ideal vs Misaligned Analyzing Power Comparison
    // ================================================================
    cout << "\n=== Creating Misalignment Effect Plot ===" << endl;
    
    // Calculate ideal A_N (both detectors at same angle)
    vector<Double_t> AN_ideal_vec, AN_misaligned_deviation_vec, relative_deviation_vec;
    
    // Determine which channel for ideal A_N calculation
    Bool_t is_elastic = (TString(channel).Contains("Elas"));
    Bool_t is_inelastic = (TString(channel).Contains("Inel"));
    
    if (is_inelastic) {
        // Load inelastic data for A_N calculation
        LoadInelasticData();
    }
    
    for (size_t i = 0; i < theta_L_vec.size(); i++) {
        Double_t theta = theta_L_vec[i];  // Use left detector angle as reference
        Double_t AN_ideal = 0.0;
        
        if (is_elastic) {
            // For elastic: use formula directly with lab angle
            AN_ideal = GetElasticAnalyzingPower(energy, theta);
        } else if (is_inelastic) {
            // For inelastic: need to convert Lab → CM
            // Use numerical inversion of kinematics
            Double_t theta_cm_min = 0.0;
            Double_t theta_cm_max = 50.0;
            Double_t tolerance = 0.001;
            
            while (theta_cm_max - theta_cm_min > tolerance) {
                Double_t theta_cm_mid = (theta_cm_min + theta_cm_max) / 2.0;
                Double_t theta_lab_test = ConvertThetaCMtoLab(energy, theta_cm_mid);
                
                if (theta_lab_test < theta) {
                    theta_cm_min = theta_cm_mid;
                } else {
                    theta_cm_max = theta_cm_mid;
                }
            }
            Double_t theta_cm = (theta_cm_min + theta_cm_max) / 2.0;
            AN_ideal = GetInelasticAnalyzingPower(theta_cm);
        }
        
        AN_ideal_vec.push_back(AN_ideal);
        
        // Calculate deviation
        Double_t AN_measured = AN_vec[i];  // From misaligned setup
        Double_t absolute_deviation = AN_measured - AN_ideal;
        Double_t relative_deviation = 0.0;
        if (TMath::Abs(AN_ideal) > 0.001) {
            relative_deviation = 100.0 * absolute_deviation / AN_ideal;  // Percentage
        }
        
        AN_misaligned_deviation_vec.push_back(absolute_deviation);
        relative_deviation_vec.push_back(relative_deviation);
    }
    
    // Create new canvas for misalignment comparison
    TCanvas* c2 = new TCanvas("c2", "Ideal vs Misaligned Analyzing Power", 1600, 900);
    c2->Divide(2, 2);  // 2x2 grid
    
    // ----------------------------------------------------------------
    // Plot 1: Ideal A_N and Measured A_N vs θ_L
    // ----------------------------------------------------------------
    c2->cd(1);
    gPad->SetGrid();
    
    // Ideal A_N (both detectors perfectly aligned)
    TGraph* gr_AN_ideal = new TGraph(npoints, &theta_L_vec[0], &AN_ideal_vec[0]);
    gr_AN_ideal->SetTitle("Ideal vs Measured A_{N};#theta [deg];A_{N}");
    gr_AN_ideal->SetLineColor(kBlack);
    gr_AN_ideal->SetLineWidth(3);
    gr_AN_ideal->SetLineStyle(1);
    gr_AN_ideal->Draw("AL");
    
    // Measured A_N (with misalignment)
    TGraphErrors* gr_AN_meas = new TGraphErrors(npoints, &theta_L_vec[0], 
                                                 &AN_vec[0], 0, &AN_err_vec[0]);
    gr_AN_meas->SetMarkerStyle(20);
    gr_AN_meas->SetMarkerSize(1.0);
    gr_AN_meas->SetMarkerColor(kRed);
    gr_AN_meas->SetLineColor(kRed);
    gr_AN_meas->Draw("PE SAME");
    
    // Reference line at theta_R
    TLine* ref_line1 = new TLine(theta_R_fixed, gr_AN_ideal->GetYaxis()->GetXmin(),
                                 theta_R_fixed, gr_AN_ideal->GetYaxis()->GetXmax());
    ref_line1->SetLineColor(kBlue);
    ref_line1->SetLineStyle(2);
    ref_line1->SetLineWidth(2);
    ref_line1->Draw("same");
    
    TLegend* leg_c2_1 = new TLegend(0.15, 0.70, 0.50, 0.88);
    leg_c2_1->SetHeader("Polarimeter Misalignment Effect", "C");
    leg_c2_1->AddEntry(gr_AN_ideal, "Ideal (both det. at same #theta)", "l");
    leg_c2_1->AddEntry(gr_AN_meas, Form("Measured (R=%.1f#circ, L variable)", theta_R_fixed), "pe");
    leg_c2_1->AddEntry(ref_line1, Form("Fixed det. at R=%.1f#circ", theta_R_fixed), "l");
    leg_c2_1->Draw();
    
    // ----------------------------------------------------------------
    // Plot 2: Absolute Deviation (A_N_measured - A_N_ideal)
    // ----------------------------------------------------------------
    c2->cd(2);
    gPad->SetGrid();
    
    TGraphErrors* gr_deviation = new TGraphErrors(npoints, &theta_L_vec[0], 
                                                   &AN_misaligned_deviation_vec[0],
                                                   0, &AN_err_vec[0]);
    gr_deviation->SetTitle("Absolute Deviation from Ideal;#theta_{L} [deg];A_{N}^{meas} - A_{N}^{ideal}");
    gr_deviation->SetMarkerStyle(21);
    gr_deviation->SetMarkerSize(1.0);
    gr_deviation->SetMarkerColor(kBlue+2);
    gr_deviation->SetLineColor(kBlue+2);
    gr_deviation->Draw("APE");
    
    // Zero line
    TLine* zero_line1 = new TLine(theta_L_vec.front(), 0.0, theta_L_vec.back(), 0.0);
    zero_line1->SetLineColor(kBlack);
    zero_line1->SetLineStyle(2);
    zero_line1->SetLineWidth(2);
    zero_line1->Draw("same");
    
    // Reference line at theta_R
    TLine* ref_line2 = new TLine(theta_R_fixed, gr_deviation->GetYaxis()->GetXmin(),
                                 theta_R_fixed, gr_deviation->GetYaxis()->GetXmax());
    ref_line2->SetLineColor(kRed);
    ref_line2->SetLineStyle(2);
    ref_line2->SetLineWidth(2);
    ref_line2->Draw("same");
    
    TLegend* leg_c2_2 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_c2_2->AddEntry(gr_deviation, "Deviation from ideal", "pe");
    leg_c2_2->AddEntry(zero_line1, "Zero deviation", "l");
    leg_c2_2->AddEntry(ref_line2, Form("R=%.1f#circ", theta_R_fixed), "l");
    leg_c2_2->Draw();
    
    // ----------------------------------------------------------------
    // Plot 3: Relative Deviation (percentage)
    // ----------------------------------------------------------------
    c2->cd(3);
    gPad->SetGrid();
    
    TGraphErrors* gr_rel_deviation = new TGraphErrors(npoints, &theta_L_vec[0], 
                                                       &relative_deviation_vec[0],
                                                       0, 0);  // No y-errors for percentage
    gr_rel_deviation->SetTitle("Relative Deviation from Ideal;#theta_{L} [deg];Relative Deviation [%]");
    gr_rel_deviation->SetMarkerStyle(22);
    gr_rel_deviation->SetMarkerSize(1.2);
    gr_rel_deviation->SetMarkerColor(kMagenta+2);
    gr_rel_deviation->SetLineColor(kMagenta+2);
    gr_rel_deviation->Draw("APE");
    
    // Zero line
    TLine* zero_line2 = new TLine(theta_L_vec.front(), 0.0, theta_L_vec.back(), 0.0);
    zero_line2->SetLineColor(kBlack);
    zero_line2->SetLineStyle(2);
    zero_line2->SetLineWidth(2);
    zero_line2->Draw("same");
    
    // Reference line at theta_R
    TLine* ref_line3 = new TLine(theta_R_fixed, gr_rel_deviation->GetYaxis()->GetXmin(),
                                 theta_R_fixed, gr_rel_deviation->GetYaxis()->GetXmax());
    ref_line3->SetLineColor(kRed);
    ref_line3->SetLineStyle(2);
    ref_line3->SetLineWidth(2);
    ref_line3->Draw("same");
    
    TLegend* leg_c2_3 = new TLegend(0.15, 0.75, 0.50, 0.88);
    leg_c2_3->AddEntry(gr_rel_deviation, "100 #times (A_{N}^{meas}-A_{N}^{ideal})/A_{N}^{ideal}", "pe");
    leg_c2_3->AddEntry(zero_line2, "Zero deviation", "l");
    leg_c2_3->AddEntry(ref_line3, Form("R=%.1f#circ", theta_R_fixed), "l");
    leg_c2_3->Draw();
    
    // ----------------------------------------------------------------
    // Plot 4: Deviation vs Misalignment (ΔΘ = θ_L - θ_R)
    // ----------------------------------------------------------------
    c2->cd(4);
    gPad->SetGrid();
    
    TGraphErrors* gr_dev_vs_misalign = new TGraphErrors(npoints, &delta_theta_vec[0], 
                                                         &AN_misaligned_deviation_vec[0],
                                                         0, &AN_err_vec[0]);
    gr_dev_vs_misalign->SetTitle("Deviation vs Misalignment;#theta_{L} - #theta_{R} [deg];A_{N}^{meas} - A_{N}^{ideal}");
    gr_dev_vs_misalign->SetMarkerStyle(21);
    gr_dev_vs_misalign->SetMarkerSize(1.0);
    gr_dev_vs_misalign->SetMarkerColor(kGreen+3);
    gr_dev_vs_misalign->SetLineColor(kGreen+3);
    gr_dev_vs_misalign->Draw("APE");
    
    // Zero line (perfect alignment)
    TLine* zero_line3 = new TLine(delta_theta_vec.front(), 0.0, delta_theta_vec.back(), 0.0);
    zero_line3->SetLineColor(kBlack);
    zero_line3->SetLineStyle(2);
    zero_line3->SetLineWidth(2);
    zero_line3->Draw("same");
    
    // Vertical line at ΔΘ = 0 (perfect alignment)
    TLine* perfect_align = new TLine(0.0, gr_dev_vs_misalign->GetYaxis()->GetXmin(),
                                     0.0, gr_dev_vs_misalign->GetYaxis()->GetXmax());
    perfect_align->SetLineColor(kRed);
    perfect_align->SetLineStyle(2);
    perfect_align->SetLineWidth(2);
    perfect_align->Draw("same");
    
    TLegend* leg_c2_4 = new TLegend(0.15, 0.75, 0.45, 0.88);
    leg_c2_4->AddEntry(gr_dev_vs_misalign, "Misalignment effect", "pe");
    leg_c2_4->AddEntry(perfect_align, "Perfect alignment", "l");
    leg_c2_4->Draw();
    
    c2->Draw();
    c2->Update();
    
    // Save misalignment comparison plot
    TString misalign_pdf = Form("%s_%s_%dMeV_misalignment.pdf", output_name, channel, energy);
    c2->SaveAs(misalign_pdf.Data());
    cout << "Misalignment plot saved: " << misalign_pdf << endl;
    
    // Print misalignment statistics
    cout << "\n=== Misalignment Effect Statistics ===" << endl;
    
    // Find maximum absolute deviation
    Double_t max_abs_dev = 0.0;
    for (size_t i = 0; i < AN_misaligned_deviation_vec.size(); i++) {
        if (TMath::Abs(AN_misaligned_deviation_vec[i]) > TMath::Abs(max_abs_dev)) {
            max_abs_dev = AN_misaligned_deviation_vec[i];
        }
    }
    
    // Find maximum relative deviation
    Double_t max_rel_dev = 0.0;
    for (size_t i = 0; i < relative_deviation_vec.size(); i++) {
        if (TMath::Abs(relative_deviation_vec[i]) > TMath::Abs(max_rel_dev)) {
            max_rel_dev = relative_deviation_vec[i];
        }
    }
    
    cout << "  Maximum absolute deviation: " << max_abs_dev << endl;
    cout << "  Maximum relative deviation: " << max_rel_dev << " %" << endl;
    
    // Find deviation at perfect alignment (θ_L = θ_R)
    for (size_t i = 0; i < theta_L_vec.size(); i++) {
        if (TMath::Abs(theta_L_vec[i] - theta_R_fixed) < 0.1) {
            Double_t dev_at_alignment = AN_misaligned_deviation_vec[i];
            cout << "  Deviation at θ_L ≈ θ_R = " << theta_R_fixed << "°: " 
                 << dev_at_alignment << " (should be ~0)" << endl;
            break;
        }
    }
    cout << "========================================" << endl;
    
    // Summary statistics
    cout << "\n========================================" << endl;
    cout << "  ANALYSIS COMPLETE" << endl;
    cout << "========================================" << endl;
    cout << "Detector positions analyzed: " << npoints << endl;
    
    cout << "\nAsymmetry Statistics:" << endl;
    Double_t asym_min = *min_element(asymmetry_vec.begin(), asymmetry_vec.end());
    Double_t asym_max = *max_element(asymmetry_vec.begin(), asymmetry_vec.end());
    cout << "  Minimum ε:   " << asym_min << endl;
    cout << "  Maximum ε:   " << asym_max << endl;
    cout << "  Range:       " << (asym_max - asym_min) << endl;
    
    cout << "\nAnalyzing Power Statistics:" << endl;
    Double_t AN_min = *min_element(AN_vec.begin(), AN_vec.end());
    Double_t AN_max = *max_element(AN_vec.begin(), AN_vec.end());
    cout << "  Minimum A_N: " << AN_min << endl;
    cout << "  Maximum A_N: " << AN_max << endl;
    cout << "  Range:       " << (AN_max - AN_min) << endl;
    cout << "  Polarization: " << polarization << "%" << endl;
    
    cout << "\nOutput files:" << endl;
    cout << "  Plot:   " << pdf_name << endl;
    cout << "  Table:  " << txt_name << endl;
    cout << "========================================\n" << endl;
}


// ====================================================================
// AnalyzeBothChannels - Analyze elastic and inelastic together
// ====================================================================
void AnalyzeBothChannels(Int_t energy = 200,
                         Int_t polarization = 80,
                         Double_t theta_R_fixed = 16.2,
                         Double_t phi_window = 0.017) {
    cout << "\n========================================" << endl;
    cout << "  ANALYZING BOTH CHANNELS" << endl;
    cout << "========================================\n" << endl;
    
    cout << "\n>>> ELASTIC CHANNEL <<<" << endl;
    AnalyzeDetectorSystematics("pC_Elas", energy, polarization, 
                              theta_R_fixed, "AN_systematic_elastic", phi_window);
    
    cout << "\n>>> INELASTIC CHANNEL <<<" << endl;
    AnalyzeDetectorSystematics("pC_Inel443", energy, polarization, 
                              theta_R_fixed, "AN_systematic_inelastic", phi_window);
    
    cout << "\n========================================" << endl;
    cout << "  BOTH CHANNELS COMPLETE" << endl;
    cout << "========================================\n" << endl;
}




// ====================================================================
// ====================================================================
// 		DILUTION Study
// ====================================================================
// ====================================================================

// ====================================================================
// AnalyzeAnalyzingPowerDilution - FIXED for compilation
// ====================================================================
void AnalyzeAnalyzingPowerDilution(Int_t energy = 200,
                                   Int_t polarization = 80,
                                   Double_t theta_R = 16.2,
                                   Double_t theta_L = 16.2,
                                   Double_t k_min = 0.0,
                                   Double_t k_max = 0.1,
                                   Double_t phi_window = 0.017)
{
    cout << "\n========================================" << endl;
    cout << "  ANALYZING POWER DILUTION STUDY" << endl;
    cout << "========================================" << endl;
    cout << "Reference: Zelenski et al., Physics of Atomic Nuclei" << endl;
    cout << "           Vol. 76, No. 12, pp. 1502-1509 (2013)" << endl;
    cout << "========================================" << endl;
    cout << "\nParameters:" << endl;
    cout << "  Energy:        " << energy << " MeV" << endl;
    cout << "  Polarization:  " << polarization << "%" << endl;
    cout << "  Right detector: " << theta_R << "°" << endl;
    cout << "  Left detector:  " << theta_L << "°" << endl;
    cout << "  k range:       [" << k_min << ", " << k_max << "]" << endl;
    cout << "========================================\n" << endl;
    
    Double_t P = polarization / 100.0;
    
    // ================================================================
    // Step 1: Load elastic files
    // ================================================================
    cout << "=== Step 1: Loading Elastic Data ===" << endl;
    
    // Format angles for filenames
    TString angle_R_str = Form("%.1f", theta_R);
    angle_R_str.ReplaceAll(".", "p");
    TString angle_L_str = Form("%.1f", theta_L);
    angle_L_str.ReplaceAll(".", "p");
    
    // Construct filenames - FIX: use .Data() for TString in Form()
    TString elas_up_R = Form("%s/pC_Elas_%dMeV_MT_P%d_SpinUp_%s.root", 
                             (const char*)DATA_INPUT_PATH, energy, polarization, angle_R_str.Data());
    TString elas_down_R = Form("%s/pC_Elas_%dMeV_MT_P%d_SpinDown_%s.root", 
                               (const char*)DATA_INPUT_PATH, energy, polarization, angle_R_str.Data());
    
    cout << "Loading elastic files:" << endl;
    cout << "  " << elas_up_R << endl;
    cout << "  " << elas_down_R << endl;
    
    // Count events - FIX: function returns void, just call it
    Int_t N_R_up_elas, N_L_up_elas, N_R_down_elas, N_L_down_elas;
    
    CountEventsInFileByPhi(elas_up_R, phi_window, true, N_R_up_elas, N_L_up_elas);
    CountEventsInFileByPhi(elas_down_R, phi_window, false, N_R_down_elas, N_L_down_elas);
    
    // If left detector at different angle
    Int_t N_R_up_elas_L = N_R_up_elas;
    Int_t N_L_up_elas_L = N_L_up_elas;
    Int_t N_R_down_elas_L = N_R_down_elas;
    Int_t N_L_down_elas_L = N_L_down_elas;
    
    if (TMath::Abs(theta_L - theta_R) > 0.01) {
        TString elas_up_L = Form("%s/pC_Elas_%dMeV_MT_P%d_SpinUp_%s.root", 
                                 (const char*)DATA_INPUT_PATH, energy, polarization, angle_L_str.Data());
        TString elas_down_L = Form("%s/pC_Elas_%dMeV_MT_P%d_SpinDown_%s.root", 
                                   (const char*)DATA_INPUT_PATH, energy, polarization, angle_L_str.Data());
        
        cout << "  " << elas_up_L << endl;
        cout << "  " << elas_down_L << endl;
        
        CountEventsInFileByPhi(elas_up_L, phi_window, true, N_R_up_elas_L, N_L_up_elas_L);
        CountEventsInFileByPhi(elas_down_L, phi_window, false, N_R_down_elas_L, N_L_down_elas_L);
    }
    
    // Calculate elastic asymmetry
    ANResult elas_asymmetry = CalculateAsymmetry(
        theta_R, theta_L,
        N_R_up_elas, N_R_down_elas,
        N_L_up_elas_L, N_L_down_elas_L
    );
    
    Double_t A_elas = AsymmetryToAnalyzingPower(elas_asymmetry.AN, P);
    Double_t A_elas_err = AsymmetryErrorToAnalyzingPowerError(elas_asymmetry.AN_error, P);
    
    cout << "\nElastic Channel Results:" << endl;
    cout << "  θ_R = " << theta_R << "°, θ_L = " << theta_L << "°" << endl;
    cout << "  N_R↑ = " << N_R_up_elas << ", N_R↓ = " << N_R_down_elas << endl;
    cout << "  N_L↑ = " << N_L_up_elas_L << ", N_L↓ = " << N_L_down_elas_L << endl;
    cout << "  Asymmetry ε = " << elas_asymmetry.AN << " ± " << elas_asymmetry.AN_error << endl;
    cout << "  A_N(elastic) = " << A_elas << " ± " << A_elas_err << endl;
    
    // ================================================================
    // Step 2: Load inelastic files
    // ================================================================
    cout << "\n=== Step 2: Loading Inelastic Data ===" << endl;
    
    TString inel_up_R = Form("%s/pC_Inel443_%dMeV_MT_P%d_SpinUp_%s.root", 
                             (const char*)DATA_INPUT_PATH, energy, polarization, angle_R_str.Data());
    TString inel_down_R = Form("%s/pC_Inel443_%dMeV_MT_P%d_SpinDown_%s.root", 
                               (const char*)DATA_INPUT_PATH, energy, polarization, angle_R_str.Data());
    
    cout << "Loading inelastic files:" << endl;
    cout << "  " << inel_up_R << endl;
    cout << "  " << inel_down_R << endl;
    
    Int_t N_R_up_inel, N_L_up_inel, N_R_down_inel, N_L_down_inel;
    
    CountEventsInFileByPhi(inel_up_R, phi_window, true, N_R_up_inel, N_L_up_inel);
    CountEventsInFileByPhi(inel_down_R, phi_window, false, N_R_down_inel, N_L_down_inel);
    
    Int_t N_R_up_inel_L = N_R_up_inel;
    Int_t N_L_up_inel_L = N_L_up_inel;
    Int_t N_R_down_inel_L = N_R_down_inel;
    Int_t N_L_down_inel_L = N_L_down_inel;
    
    if (TMath::Abs(theta_L - theta_R) > 0.01) {
        TString inel_up_L = Form("%s/pC_Inel443_%dMeV_MT_P%d_SpinUp_%s.root", 
                                 (const char*)DATA_INPUT_PATH, energy, polarization, angle_L_str.Data());
        TString inel_down_L = Form("%s/pC_Inel443_%dMeV_MT_P%d_SpinDown_%s.root", 
                                   (const char*)DATA_INPUT_PATH, energy, polarization, angle_L_str.Data());
        
        cout << "  " << inel_up_L << endl;
        cout << "  " << inel_down_L << endl;
        
        CountEventsInFileByPhi(inel_up_L, phi_window, true, N_R_up_inel_L, N_L_up_inel_L);
        CountEventsInFileByPhi(inel_down_L, phi_window, false, N_R_down_inel_L, N_L_down_inel_L);
    }
    
    ANResult inel_asymmetry = CalculateAsymmetry(
        theta_R, theta_L,
        N_R_up_inel, N_R_down_inel,
        N_L_up_inel_L, N_L_down_inel_L
    );
    
    Double_t A_inel = AsymmetryToAnalyzingPower(inel_asymmetry.AN, P);
    Double_t A_inel_err = AsymmetryErrorToAnalyzingPowerError(inel_asymmetry.AN_error, P);
    
    cout << "\nInelastic Channel Results:" << endl;
    cout << "  θ_R = " << theta_R << "°, θ_L = " << theta_L << "°" << endl;
    cout << "  N_R↑ = " << N_R_up_inel << ", N_R↓ = " << N_R_down_inel << endl;
    cout << "  N_L↑ = " << N_L_up_inel_L << ", N_L↓ = " << N_L_down_inel_L << endl;
    cout << "  Asymmetry ε = " << inel_asymmetry.AN << " ± " << inel_asymmetry.AN_error << endl;
    cout << "  A_N(inelastic) = " << A_inel << " ± " << A_inel_err << endl;
    
    // ================================================================
    // Step 3: Calculate A_Σ with proper error propagation
    // ================================================================
    cout << "\n=== Step 3: Calculating Dilution with Error Propagation ===" << endl;
    cout << "Formula: A_Σ = (A_elas + A_inel × k) / (1 + k)" << endl;
    cout << "\nError propagation:" << endl;
    cout << "  δA_Σ = (1/(1+k)) × √(δA_elas² + k² × δA_inel²)" << endl;
    cout << endl;
    
    const Int_t npoints = 100;
    vector<Double_t> k_values, A_sigma_values, A_sigma_errors;
    
    for (Int_t i = 0; i < npoints; i++) {
        Double_t k = k_min + (k_max - k_min) * i / (npoints - 1);
        
        // Zelenski formula
        Double_t A_sigma = (A_elas + A_inel * k) / (1.0 + k);
        
        // Proper error propagation
        Double_t term1 = 1.0 / (1.0 + k);
        Double_t term2 = k / (1.0 + k);
        Double_t A_sigma_err = TMath::Sqrt(term1*term1 * A_elas_err*A_elas_err + 
                                           term2*term2 * A_inel_err*A_inel_err);
        
        k_values.push_back(k);
        A_sigma_values.push_back(A_sigma);
        A_sigma_errors.push_back(A_sigma_err);
        
        if (i % 20 == 0 || i == npoints-1) {
            Double_t change_percent = 100.0 * (A_sigma - A_elas) / A_elas;
            cout << Form("  k = %.4f:  A_Σ = %.4f ± %.4f  (change: %+.2f%%)", 
                        k, A_sigma, A_sigma_err, change_percent) << endl;
        }
    }
    
    // ================================================================
    // Step 4: Create plots
    // ================================================================
    cout << "\n=== Step 4: Creating Dilution Plot ===" << endl;
    
    TCanvas* c = new TCanvas("c_dilution", "Analyzing Power Dilution", 1200, 1000);
    c->Divide(1, 2);
    
    // Panel 1: A_Σ vs k
    c->cd(1);
    gPad->SetLogx(1);
    gPad->SetGrid();
    gPad->SetMargin(0.12, 0.05, 0.12, 0.08);
    
    Double_t y_min = TMath::Min(A_elas, A_inel) * 0.95;
    Double_t y_max = TMath::Max(A_elas, A_inel) * 1.05;
    
    TGraphErrors* gr_sigma = new TGraphErrors(npoints, &k_values[0], 
                                              &A_sigma_values[0], 
                                              0, &A_sigma_errors[0]);
    gr_sigma->SetTitle(Form("Analyzing Power Dilution (E=%d MeV, P=%d%%, #theta=%.1f#circ);k = N_{inel}/N_{elas};A_{#Sigma}", 
                           energy, polarization, theta_R));
    gr_sigma->SetMarkerStyle(20);
    gr_sigma->SetMarkerSize(0.6);
    gr_sigma->SetMarkerColor(kBlue+2);
    gr_sigma->SetLineColor(kBlue+2);
    gr_sigma->SetLineWidth(2);
    gr_sigma->GetYaxis()->SetRangeUser(y_min, y_max);
    gr_sigma->Draw("APE");
    
    Double_t x_min_plot = (k_min > 0) ? k_min : 0.001;
    TLine* line_elas = new TLine(0, A_elas, k_max, A_elas);
    line_elas->SetLineColor(kGreen+2);
    line_elas->SetLineWidth(3);
    line_elas->SetLineStyle(2);
    line_elas->Draw("same");
    
    TLine* line_inel = new TLine(0, A_inel, k_max, A_inel);
    line_inel->SetLineColor(kRed+2);
    line_inel->SetLineWidth(3);
    line_inel->SetLineStyle(2);
    line_inel->Draw("same");
    
    TLegend* leg1 = new TLegend(0.15, 0.15, 0.55, 0.35);
    leg1->SetHeader("Zelenski Formula: A_{#Sigma} = (A_{1} + A_{2}#upoint k)/(1+k)", "C");
    leg1->AddEntry(gr_sigma, "A_{#Sigma}(k) - Mixed detection", "lpe");
    leg1->AddEntry(line_elas, Form("A_{elastic} = %.4f #pm %.4f", A_elas, A_elas_err), "l");
    leg1->AddEntry(line_inel, Form("A_{inelastic} = %.4f #pm %.4f", A_inel, A_inel_err), "l");
    leg1->Draw();
    
    TLatex* latex1 = new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.035);
    latex1->DrawLatex(0.60, 0.85, Form("p + ^{12}C at %d MeV", energy));
    latex1->DrawLatex(0.60, 0.80, Form("#theta_{R} = %.1f#circ, #theta_{L} = %.1f#circ", theta_R, theta_L));
    
    // Panel 2: Relative change
    c->cd(2);
    gPad->SetLogx(1);
    gPad->SetGrid();
    gPad->SetMargin(0.12, 0.05, 0.12, 0.08);
    
    vector<Double_t> change_rel;
    for (size_t i = 0; i < k_values.size(); i++) {
        change_rel.push_back(100.0 * (A_sigma_values[i] - A_elas) / A_elas);
    }
    
    TGraph* gr_change = new TGraph(npoints, &k_values[0], &change_rel[0]);
    gr_change->SetTitle("Relative Change from Pure Elastic;k = N_{inel}/N_{elas};(A_{#Sigma} - A_{elastic})/A_{elastic} [%]");
    gr_change->SetMarkerStyle(21);
    gr_change->SetMarkerSize(0.8);
    gr_change->SetMarkerColor(kMagenta+2);
    gr_change->SetLineColor(kMagenta+2);
    gr_change->SetLineWidth(2);
    
    Double_t change_max = *max_element(change_rel.begin(), change_rel.end());
    Double_t change_min = *min_element(change_rel.begin(), change_rel.end());
    Double_t change_range = TMath::Max(TMath::Abs(change_max), TMath::Abs(change_min));
    gr_change->GetYaxis()->SetRangeUser(-change_range*1.2, change_range*1.2);
    
    gr_change->Draw("APL");
    
    TLine* zero_line = new TLine(x_min_plot, 0.0, k_max, 0.0);
    zero_line->SetLineColor(kBlack);
    zero_line->SetLineStyle(2);
    zero_line->SetLineWidth(2);
    zero_line->Draw("same");
    
    TLegend* leg2 = new TLegend(0.15, 0.75, 0.65, 0.88);
    leg2->AddEntry(gr_change, "Change = (A_{#Sigma} - A_{elas})/A_{elas}", "lp");
    leg2->AddEntry(zero_line, "No change (pure elastic, k=0)", "l");
    leg2->Draw();
    
    c->Update();
    
    TString pdf_name = Form("AN_dilution_%dMeV_P%d_theta%.1f.pdf", 
                           energy, polarization, theta_R);
    pdf_name.ReplaceAll(".", "p");
    c->SaveAs(pdf_name.Data());
    cout << "Plot saved: " << pdf_name << endl;
    
    // ================================================================
    // Step 5: Summary
    // ================================================================
    cout << "\n========================================" << endl;
    cout << "  DILUTION ANALYSIS SUMMARY" << endl;
    cout << "========================================" << endl;
    cout << "\nMeasured Analyzing Powers:" << endl;
    cout << "  A_N(elastic)   = " << A_elas << " ± " << A_elas_err << endl;
    cout << "  A_N(inelastic) = " << A_inel << " ± " << A_inel_err << endl;
    cout << "  Ratio: A_inel/A_elas = " << (A_inel/A_elas) << endl;
    
    cout << "\nChange at Representative k Values:" << endl;
    for (Double_t k_test : {0.001, 0.01, 0.05, 0.10}) {
        if (k_test >= k_min && k_test <= k_max) {
            Double_t A_sig = (A_elas + A_inel * k_test) / (1.0 + k_test);
            Double_t change = 100.0 * (A_sig - A_elas) / A_elas;
            
            Double_t t1 = 1.0 / (1.0 + k_test);
            Double_t t2 = k_test / (1.0 + k_test);
            Double_t A_sig_err = TMath::Sqrt(t1*t1*A_elas_err*A_elas_err + t2*t2*A_inel_err*A_inel_err);
            
            cout << Form("  k = %.3f:  A_Σ = %.4f ± %.4f, change = %+.2f%%", 
                        k_test, A_sig, A_sig_err, change) << endl;
        }
    }
    
    cout << "\nPhysical Interpretation:" << endl;
    if (A_inel > A_elas) {
        cout << "  A_inel > A_elas → Inelastic contamination INCREASES A_Σ" << endl;
        cout << "  This is ENHANCEMENT, not dilution!" << endl;
    } else if (A_inel < A_elas) {
        cout << "  A_inel < A_elas → Inelastic contamination DECREASES A_Σ" << endl;
        cout << "  True dilution effect" << endl;
    }
    cout << "========================================\n" << endl;
}


// ====================================================================
// QuickAnalyzingPowerLookup - Fast A_N extraction at single angle
// ====================================================================
// Convenience function to quickly extract A_N for both channels
// at a specific detector configuration.
//
// Parameters:
//   channel      - "pC_Elas" or "pC_Inel443"
//   energy       - Beam energy [MeV]
//   polarization - Beam polarization [%]
//   theta_R      - Right detector angle [degrees]
//   theta_L      - Left detector angle [degrees]
//   phi_window   - Azimuthal acceptance window [rad]
//
// Returns: Prints A_N value and error to console
// ====================================================================
void QuickAnalyzingPowerLookup(const char* channel,
                               Int_t energy = 200,
                               Int_t polarization = 80,
                               Double_t theta_R = 16.2,
                               Double_t theta_L = 16.2,
                               Double_t phi_window = 0.017)
{
    cout << "\n========================================" << endl;
    cout << "  QUICK A_N LOOKUP" << endl;
    cout << "========================================" << endl;
    cout << "Channel:       " << channel << endl;
    cout << "Energy:        " << energy << " MeV" << endl;
    cout << "Polarization:  " << polarization << "%" << endl;
    cout << "θ_R:           " << theta_R << "°" << endl;
    cout << "θ_L:           " << theta_L << "°" << endl;
    cout << "========================================" << endl;
    
    Double_t P = polarization / 100.0;
    
    // Load data
    map<Double_t, DetectorCounts> data = FindDataFiles(channel, energy, 
                                                       polarization, phi_window);
    
    if (data.find(theta_R) == data.end()) {
        cerr << "ERROR: No data at θ_R = " << theta_R << "°" << endl;
        return;
    }
    
    if (theta_L != theta_R && data.find(theta_L) == data.end()) {
        cerr << "ERROR: No data at θ_L = " << theta_L << "°" << endl;
        return;
    }
    
    // Get counts
    DetectorCounts data_R = data[theta_R];
    DetectorCounts data_L = (theta_L == theta_R) ? data_R : data[theta_L];
    
    // Calculate asymmetry
    ANResult asymmetry = CalculateAsymmetry(
        theta_R, theta_L,
        data_R.N_right_up, data_R.N_right_down,
        data_L.N_left_up, data_L.N_left_down
    );
    
    Double_t AN = AsymmetryToAnalyzingPower(asymmetry.AN, P);
    Double_t AN_err = AsymmetryErrorToAnalyzingPowerError(asymmetry.AN_error, P);
    
    // Results
    cout << "\nResults:" << endl;
    cout << "  Counts: N_R↑=" << data_R.N_right_up << ", N_R↓=" << data_R.N_right_down;
    cout << ", N_L↑=" << data_L.N_left_up << ", N_L↓=" << data_L.N_left_down << endl;
    cout << "  Asymmetry ε = " << asymmetry.AN << " ± " << asymmetry.AN_error << endl;
    cout << "  A_N = " << AN << " ± " << AN_err << endl;
    cout << "  Significance: " << TMath::Abs(AN/AN_err) << " σ" << endl;
    cout << "========================================\n" << endl;
}











// ====================================================================
// PrintAnalysisHelp - Show available analysis functions
// ====================================================================
void PrintAnalysisHelp() {
    cout << "\n========================================" << endl;
    cout << "  ANALYSIS FUNCTIONS GUIDE" << endl;
    cout << "========================================" << endl;
    
    cout << "\n1. Detector Systematics (Misalignment Study):" << endl;
    cout << "   AnalyzeDetectorSystematics(\"pC_Elas\", 200, 80, 16.2)" << endl;
    cout << "   AnalyzeDetectorSystematics(\"pC_Inel443\", 200, 80, 16.2)" << endl;
    cout << "   → Creates AN_systematic_*.pdf (4 panels)" << endl;
    cout << "   → Creates AN_systematic_*_misalignment.pdf (NEW!)" << endl;
    
    cout << "\n2. Both Channels Together:" << endl;
    cout << "   AnalyzeBothChannels(200, 80, 16.2)" << endl;
    cout << "   → Runs elastic + inelastic sequentially" << endl;
    
    cout << "\n3. Analyzing Power Dilution Study (NEW!):" << endl;
    cout << "   AnalyzeAnalyzingPowerDilution(200, 80, 16.2, 16.2)" << endl;
    cout << "   → Studies mixed elastic+inelastic detection" << endl;
    cout << "   → Based on Zelenski formula: A_Σ = (A_1 + A_2·k)/(1+k)" << endl;
    cout << "   → Creates AN_dilution_*.pdf" << endl;
    
    cout << "\n   Custom k range:" << endl;
    cout << "   AnalyzeAnalyzingPowerDilution(200, 80, 16.2, 16.2, 0.0, 0.5)" << endl;
    cout << "   → Studies k from 0 to 0.5 (0% to 50% contamination)" << endl;
    
    cout << "\n4. Quick A_N Lookup:" << endl;
    cout << "   QuickAnalyzingPowerLookup(\"pC_Elas\", 200, 80, 16.2, 16.2)" << endl;
    cout << "   QuickAnalyzingPowerLookup(\"pC_Inel443\", 200, 80, 16.2, 16.2)" << endl;
    cout << "   → Fast extraction of A_N at specific angle" << endl;
    
    cout << "\n========================================" << endl;
    cout << "PARAMETERS GUIDE" << endl;
    cout << "========================================" << endl;
    cout << "  energy:        Beam energy [MeV] (150-240)" << endl;
    cout << "  polarization:  Beam polarization [%] (0-100)" << endl;
    cout << "  theta_R:       Right detector angle [degrees]" << endl;
    cout << "  theta_L:       Left detector angle [degrees]" << endl;
    cout << "  k_min, k_max:  Contamination ratio range" << endl;
    cout << "  phi_window:    Azimuthal window [rad] (default: 0.017)" << endl;
    
    cout << "\n========================================" << endl;
    cout << "DILUTION STUDY DETAILS" << endl;
    cout << "========================================" << endl;
    cout << "Physical meaning of k:" << endl;
    cout << "  k = N_inel / N_elas  (ratio of event rates)" << endl;
    cout << "\nTypical scenarios:" << endl;
    cout << "  k = 0.00  → Pure elastic (no contamination)" << endl;
    cout << "  k = 0.01  → 1% inelastic contamination" << endl;
    cout << "  k = 0.05  → 5% inelastic contamination" << endl;
    cout << "  k = 0.10  → 10% inelastic contamination" << endl;
    cout << "  k = 1.00  → Equal elastic and inelastic" << endl;
    
    cout << "\nInterpretation:" << endl;
    cout << "  If A_inel > A_elas: contamination increases A_Σ" << endl;
    cout << "  If A_inel < A_elas: contamination decreases A_Σ (dilution)" << endl;
    
    cout << "\nReference:" << endl;
    cout << "  Zelenski et al., Physics of Atomic Nuclei" << endl;
    cout << "  Vol. 76, No. 12, pp. 1502-1509 (2013)" << endl;
    cout << "  https://link.springer.com/article/10.1134/S1063778813120156" << endl;
    
    cout << "\n========================================" << endl;
    cout << "KEY TECHNICAL NOTE" << endl;
    cout << "========================================" << endl;
    cout << "Events are separated by φ angle:" << endl;
    cout << "  Right detector: φ ≈ 0   (|φ| < φ_window)" << endl;
    cout << "  Left detector:  φ ≈ π   (|φ-π| < φ_window)" << endl;
    cout << "This ensures proper correlation measurement!" << endl;
    cout << "========================================\n" << endl;
}