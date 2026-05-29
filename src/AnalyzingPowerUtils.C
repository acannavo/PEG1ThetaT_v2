// ====================================================================
// AnalyzingPowerUtils.C
// ====================================================================
// Implementation of analyzing power utility functions
// ====================================================================

#include "AnalyzingPowerUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "PParticle.h"
#include "TMath.h"
#include <iostream>

using namespace std;

// ====================================================================
// CountEventsInFileByPhi - Count events separated by φ angle
// ====================================================================
// Separates events into right detector (φ ≈ 0) and left detector (φ ≈ π)
// Uses detector acceptance windows to determine which events go where
//
// Parameters:
//   filename   - ROOT file to analyze
//   phi_window - Azimuthal acceptance window [radians]
//   is_spinup  - True for spin-up, false for spin-down (currently unused but kept for compatibility)
//   N_right    - Output: events with φ ≈ 0 (right detector)
//   N_left     - Output: events with φ ≈ π (left detector)
// ====================================================================
void CountEventsInFileByPhi(const TString& filename, 
                            Double_t phi_window,
                            Bool_t is_spinup,
                            Int_t& N_right, 
                            Int_t& N_left) {
    N_right = 0;
    N_left = 0;
    
    TFile* file = TFile::Open(filename.Data());
    if (!file || file->IsZombie()) {
        cerr << "ERROR: Cannot open " << filename << endl;
        return;
    }
    
    TTree* tree = (TTree*)file->Get("data");
    if (!tree) {
        cerr << "ERROR: No 'data' tree in " << filename << endl;
        file->Close();
        delete file;
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
        // Events outside both windows are ignored
    }
    
    file->Close();
    delete file;
    delete particles;
}

// ====================================================================
// CalculateAsymmetry - Compute spin asymmetry from correlation
// ====================================================================
// Formula: ε = [√(N_R_up × N_L_down) - √(N_L_up × N_R_down)] /
//              [√(N_R_up × N_L_down) + √(N_L_up × N_R_down)]
//
// This calculates the SPIN ASYMMETRY (ε), NOT the analyzing power (A_N)
// To get A_N, divide by beam polarization: A_N = ε / P
//
// The result is stored in result.AN for backward compatibility, but
// represents asymmetry, not analyzing power!
//
// Returns: ANResult with asymmetry in the AN field
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
    
    // Check for zero counts
    if (N_R_up == 0 || N_R_down == 0 || N_L_up == 0 || N_L_down == 0) {
        cerr << "WARNING: Zero counts detected in CalculateAsymmetry!" << endl;
        cerr << "  N_R_up=" << N_R_up << ", N_R_down=" << N_R_down 
             << ", N_L_up=" << N_L_up << ", N_L_down=" << N_L_down << endl;
        result.AN = 0;
        result.AN_error = 999.0;
        return result;
    }
    
    // Cast to Double_t BEFORE multiplication to prevent integer overflow
    // This is critical for large event counts (>46340 can overflow Int_t)
    Double_t term1 = TMath::Sqrt((Double_t)N_R_up * (Double_t)N_L_down);
    Double_t term2 = TMath::Sqrt((Double_t)N_L_up * (Double_t)N_R_down);
    
    Double_t numerator = term1 - term2;
    Double_t denominator = term1 + term2;

    // Check for invalid math
    if (denominator == 0 || std::isnan(numerator) || std::isnan(denominator)) {
        cerr << "WARNING: Invalid math in CalculateAsymmetry at theta_L = " << theta_L << endl;
        result.AN = 0;
        result.AN_error = 999;
        return result;
    }
    
    // Calculate asymmetry (stored in result.AN)
    result.AN = numerator / denominator;
    
    // ================================================================
    // Correct error propagation for asymmetry measurement
    // ================================================================
    // For ε = (√(a·b) - √(c·d)) / (√(a·b) + √(c·d))
    // where a=N_R_up, b=N_L_down, c=N_L_up, d=N_R_down
    
    // Relative errors on sqrt terms (Poisson statistics: δN = √N)
    // For √(a·b): relative error = 0.5 × √(1/a + 1/b)
    Double_t rel_err_term1 = 0.5 * TMath::Sqrt(1.0/N_R_up + 1.0/N_L_down);
    Double_t rel_err_term2 = 0.5 * TMath::Sqrt(1.0/N_L_up + 1.0/N_R_down);

    // Absolute errors on terms
    Double_t err_term1 = term1 * rel_err_term1;
    Double_t err_term2 = term2 * rel_err_term2;

    // Error on numerator (uncorrelated terms add in quadrature)
    Double_t err_numerator = TMath::Sqrt(err_term1*err_term1 + err_term2*err_term2);
    
    // Error on denominator (same terms as numerator)
    Double_t err_denominator = err_numerator;

    // Full error propagation for quotient f = (a-b)/(a+b)
    // Using: δf = √[(δa/a)² + (δb/b)²] × |f|
    Double_t rel_err = TMath::Sqrt(TMath::Power(err_numerator/numerator, 2) + 
                                    TMath::Power(err_denominator/denominator, 2));
    result.AN_error = TMath::Abs(result.AN) * rel_err;
    
    return result;
}

// ====================================================================
// AsymmetryToAnalyzingPower - Convert asymmetry to analyzing power
// ====================================================================
// Simple conversion: A_N = ε / P
//
// Parameters:
//   asymmetry    - Measured spin asymmetry ε
//   polarization - Beam polarization (0.0 to 1.0, e.g., 0.80 for 80%)
//
// Returns: Analyzing power A_N
// ====================================================================
Double_t AsymmetryToAnalyzingPower(Double_t asymmetry, Double_t polarization) {
    if (polarization == 0) {
        cerr << "ERROR: Polarization is zero in AsymmetryToAnalyzingPower!" << endl;
        return 0;
    }
    return asymmetry / polarization;
}

// ====================================================================
// AsymmetryErrorToAnalyzingPowerError - Convert error
// ====================================================================
// Error propagation: δ(A_N) = δε / P
//
// Parameters:
//   asymmetry_error - Error on asymmetry measurement
//   polarization    - Beam polarization (0.0 to 1.0)
//
// Returns: Error on analyzing power
// ====================================================================
Double_t AsymmetryErrorToAnalyzingPowerError(Double_t asymmetry_error, Double_t polarization) {
    if (polarization == 0) {
        cerr << "ERROR: Polarization is zero in AsymmetryErrorToAnalyzingPowerError!" << endl;
        return 999.0;
    }
    return asymmetry_error / polarization;
}

// ====================================================================
// CalculateAnalyzingPower - Full calculation with error propagation
// ====================================================================
// Converts an asymmetry result to analyzing power with proper error propagation
//
// Parameters:
//   asymmetry_result       - Result from CalculateAsymmetry()
//   polarization           - Beam polarization (0.0 to 1.0)
//   analyzing_power        - Output: A_N value
//   analyzing_power_error  - Output: Error on A_N
// ====================================================================
void CalculateAnalyzingPower(const ANResult& asymmetry_result,
                             Double_t polarization,
                             Double_t& analyzing_power,
                             Double_t& analyzing_power_error) {
    analyzing_power = AsymmetryToAnalyzingPower(asymmetry_result.AN, polarization);
    analyzing_power_error = AsymmetryErrorToAnalyzingPowerError(asymmetry_result.AN_error, polarization);
}