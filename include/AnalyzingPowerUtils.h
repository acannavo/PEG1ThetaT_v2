// ====================================================================
// AnalyzingPowerUtils.h
// ====================================================================
// Utility functions and data structures for analyzing power measurements
// using spin-correlation method
//
// This module provides:
// - Data structures for storing results
// - Event counting by azimuthal angle (phi)
// - Spin asymmetry calculation
// - Conversion from asymmetry to analyzing power
//
// Usage:
//   #include "AnalyzingPowerUtils.h"
//
// Author: PEG Project
// Date: January 2026
// ====================================================================

#ifndef ANALYZINGPOWERUTILS_H
#define ANALYZINGPOWERUTILS_H

#include "TString.h"
#include "TMath.h"

// ====================================================================
// Data Structures
// ====================================================================

// Store event counts for one detector configuration
struct DetectorCounts {
    Double_t theta_det;      // Detector angle [degrees]
    Int_t N_right_up;        // Right detector (φ≈0), spin-up
    Int_t N_right_down;      // Right detector (φ≈0), spin-down
    Int_t N_left_up;         // Left detector (φ≈π), spin-up
    Int_t N_left_down;       // Left detector (φ≈π), spin-down
    
    DetectorCounts() : theta_det(0), N_right_up(0), N_right_down(0), 
                       N_left_up(0), N_left_down(0) {}
};

// Store analyzing power result
struct ANResult {
    Double_t theta_L;        // Left detector angle [degrees]
    Double_t theta_R;        // Right detector angle [degrees]
    Double_t AN;             // Analyzing power (or asymmetry - see comments)
    Double_t AN_error;       // Statistical error
    Int_t N_R_up;            // Right detector, spin-up
    Int_t N_R_down;          // Right detector, spin-down
    Int_t N_L_up;            // Left detector, spin-up
    Int_t N_L_down;          // Left detector, spin-down
    
    ANResult() : theta_L(0), theta_R(0), AN(0), AN_error(0),
                 N_R_up(0), N_R_down(0), N_L_up(0), N_L_down(0) {}
};

// ====================================================================
// Function Declarations
// ====================================================================

// Count events separated by φ angle
void CountEventsInFileByPhi(const TString& filename, 
                            Double_t phi_window,
                            Bool_t is_spinup,
                            Int_t& N_right, 
                            Int_t& N_left);

// Calculate spin asymmetry from event counts
ANResult CalculateAsymmetry(Double_t theta_R, Double_t theta_L,
                            Int_t N_R_up, Int_t N_R_down,
                            Int_t N_L_up, Int_t N_L_down);

// Convert asymmetry to analyzing power
Double_t AsymmetryToAnalyzingPower(Double_t asymmetry, Double_t polarization);

// Convert asymmetry error to analyzing power error
Double_t AsymmetryErrorToAnalyzingPowerError(Double_t asymmetry_error, Double_t polarization);

// Full calculation: asymmetry → analyzing power with error propagation
void CalculateAnalyzingPower(const ANResult& asymmetry_result,
                             Double_t polarization,
                             Double_t& analyzing_power,
                             Double_t& analyzing_power_error);

#endif // ANALYZINGPOWERUTILS_H