// ElasticScattering.h
// Elastic scattering cross section data for p+C

#ifndef ELASTIC_SCATTERING_H
#define ELASTIC_SCATTERING_H

#include "TMath.h"

// ====================================================================
// Include optical model cross section data (log scale)
// These files contain parameterizations of dσ/dΩ from optical model
// calculations for p+C elastic scattering at different beam energies
// ====================================================================
#include "XSLog150MeV.C"
#include "XSLog160MeV.C"
#include "XSLog170MeV.C"
#include "XSLog180MeV.C"
#include "XSLog190MeV.C"
#include "XSLog200MeV.C"
#include "XSLog210MeV.C"
#include "XSLog220MeV.C"
#include "XSLog230MeV.C"
#include "XSLog240MeV.C"

// ====================================================================
// t-based cross section and analyzing power splines
// ====================================================================
#include "XStSpline.C"
#include "APSpline.C"

// ====================================================================
// t-based cross section and analyzing power splines
// These provide energy-independent parameterizations as functions of
// momentum transfer t (instead of theta_CM at specific energies)
// ====================================================================
double XStSpline(double t);  // Cross section vs t [(MeV/c)²]
double APSpline(double t);   // Analyzing power vs t [(MeV/c)²]

// ====================================================================
// Convert log cross section to linear scale
// The input files provide log10(σ), so we convert: σ = 10^(log10(σ))
// ====================================================================
double fXS150Op(double x);
double fXS160Op(double x);
double fXS170Op(double x);
double fXS180Op(double x);
double fXS190Op(double x);
double fXS200Op(double x);
double fXS210Op(double x);
double fXS220Op(double x);
double fXS230Op(double x);
double fXS240Op(double x);

// ====================================================================
// Helper function: Get cross section formula string for given energy
// Returns the appropriate function name based on beam kinetic energy
// ====================================================================
const char* GetXSFormula(Double_t ekin);




// ====================================================================
// Analyzing Power Parameters for Elastic Scattering
// Based on: https://pos.sissa.it/324/024/pdf
// Wissink et al., Phys. Rev. C 45, R504 (1992)
// ====================================================================
const Double_t T0_AN = 189.0;        // MeV - energy at maximum A_N
const Double_t theta0_AN = 17.3;     // degrees - angle at maximum A_N (lab frame)
const Double_t alpha_AN = 1.21e-4;   // MeV^-2
const Double_t beta_AN = 1.61e-3;    // MeV^-1 deg^-1
const Double_t gamma_AN = 1.00e-2;   // deg^-2

// ====================================================================
// GetElasticAnalyzingPower: Calculate A_N for elastic scattering
// ====================================================================
// Formula: A_N(T,θ) = 1 - α(T-T₀)² - β(T-T₀)(θ-θ₀) - γ(θ-θ₀)²
// 
// Input:
//   - ekin: beam kinetic energy [MeV]
//   - theta_lab_deg: scattering angle in lab frame [degrees]
// 
// Output:
//   - analyzing power A_N (dimensionless, typically 0 to 1)
// ====================================================================
Double_t GetElasticAnalyzingPower(Double_t ekin, Double_t theta_lab_deg);




#endif // ELASTIC_SCATTERING_H