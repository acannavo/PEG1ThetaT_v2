// InelasticScattering.C
// Implementation of inelastic scattering data loading and access

#include "InelasticScattering.h"
#include "TSystem.h"
#include "TH1D.h"           
#include "TLorentzVector.h" 
#include "TMath.h"          
#include "PhysicalConstants.h" 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

using namespace std;

// ====================================================================
// Global TGraph objects (definition)
// ====================================================================
TGraph* g_inelastic_xs = nullptr;
TGraph* g_inelastic_AN = nullptr;

// ====================================================================
// LoadInelasticData: Load digitized inelastic scattering data
// ====================================================================
void LoadInelasticData()
{
    cout << "\n========================================" << endl;
    cout << "Loading Inelastic Data" << endl;
    cout << "========================================" << endl;
    
    // ================================================================
    // Load Cross Section
    // ================================================================
    const char* xs_file = "inelastic_crosssection_443_200MeV.csv";
    ifstream file_xs(xs_file);
    
    if (!file_xs.is_open()) {
        cerr << "ERROR: Cannot open " << xs_file << endl;
        cerr << "Current directory: " << gSystem->pwd() << endl;
        cerr << "Make sure the CSV file is in the same directory!" << endl;
        exit(1);
    }
    
    vector<double> theta_xs, sigma_xs;
    string line;
    
    while (getline(file_xs, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        double theta, sigma;
        if (sscanf(line.c_str(), "%lf,%lf", &theta, &sigma) == 2 ||
            sscanf(line.c_str(), "%lf, %lf", &theta, &sigma) == 2) {
            theta_xs.push_back(theta);
            sigma_xs.push_back(sigma);
        }
    }
    file_xs.close();
    
    g_inelastic_xs = new TGraph(theta_xs.size(), &theta_xs[0], &sigma_xs[0]);
    g_inelastic_xs->SetName("g_inelastic_xs");
    cout << "✓ Cross section: " << theta_xs.size() << " points" << endl;
    
    // ================================================================
    // Load Analyzing Power
    // ================================================================
    const char* an_file = "inelastic_analyzingpower_443_200MeV.csv";
    ifstream file_an(an_file);
    
    if (!file_an.is_open()) {
        cerr << "ERROR: Cannot open " << an_file << endl;
        cerr << "Current directory: " << gSystem->pwd() << endl;
        cerr << "Make sure the CSV file is in the same directory!" << endl;
        exit(1);
    }
    
    vector<double> theta_an, an_values;
    
    while (getline(file_an, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        double theta, an;
        if (sscanf(line.c_str(), "%lf,%lf", &theta, &an) == 2 ||
            sscanf(line.c_str(), "%lf, %lf", &theta, &an) == 2) {
            theta_an.push_back(theta);
            an_values.push_back(an);
        }
    }
    file_an.close();
    
    g_inelastic_AN = new TGraph(theta_an.size(), &theta_an[0], &an_values[0]);
    g_inelastic_AN->SetName("g_inelastic_AN");
    cout << "✓ Analyzing power: " << theta_an.size() << " points" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// GetInelasticCrossSection: Interpolate cross section at given angle
// ====================================================================
Double_t GetInelasticCrossSection(Double_t theta_cm_deg)
{
    if (!g_inelastic_xs) {
        cerr << "ERROR: Call LoadInelasticData() first!" << endl;
        exit(1);
    }
    return g_inelastic_xs->Eval(theta_cm_deg);
}

// ====================================================================
// GetInelasticAnalyzingPower: Interpolate A_N at given angle
// ====================================================================
Double_t GetInelasticAnalyzingPower(Double_t theta_cm_deg)
{
    if (!g_inelastic_AN) {
        cerr << "ERROR: Call LoadInelasticData() first!" << endl;
        exit(1);
    }
    return g_inelastic_AN->Eval(theta_cm_deg);
}



// ====================================================================
// CalculateCMMomentumInelastic: Calculate CM momentum for inelastic
// ====================================================================
Double_t CalculateCMMomentumInelastic(Double_t ekin, Double_t E_ex)
{
    // Convert to GeV
    Double_t ekinGeV = ekin / 1000.;
    
    // Initial proton momentum
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    
    // Initial 4-momentum: beam proton + target carbon at rest
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    
    // Mandelstam s (total energy squared in CM)
    Double_t s = iState.Mag2();
    
    // For inelastic: final state has carbon with excitation energy
    // Effective mass of excited carbon: m_C* = m_C + E_ex
    Double_t mC_star = mC + E_ex;
    
    // CM momentum formula with excited final state
    Double_t pcm_inel = TMath::Sqrt((s - (mp + mC_star)*(mp + mC_star)) * 
                                     (s - (mp - mC_star)*(mp - mC_star)) / (4.*s));
    
    return pcm_inel;
}

// ====================================================================
// CreateInelasticSamplingHistogram: Create histogram for importance sampling
// ====================================================================
TH1D* CreateInelasticSamplingHistogram(Double_t theta_cm_min, Double_t theta_cm_max, 
                                       Int_t nbins)
{
    if (!g_inelastic_xs) {
        cerr << "ERROR: Call LoadInelasticData() first!" << endl;
        return nullptr;
    }
    
    TH1D* h = new TH1D("h_inelastic_sampling", 
                       "Inelastic cross section sampling", 
                       nbins, theta_cm_min, theta_cm_max);
    
    // Fill histogram with interpolated cross section values
    for (Int_t ibin = 1; ibin <= nbins; ibin++) {
        Double_t theta_cm = h->GetBinCenter(ibin);
        Double_t xs = GetInelasticCrossSection(theta_cm);
        h->SetBinContent(ibin, xs);
    }
    
    return h;
}