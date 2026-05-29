// ElasticScattering.C
// Implementation of elastic scattering cross section conversions

#include "ElasticScattering.h"
#include <iostream>

using namespace std;

// ====================================================================
// Convert log cross section to linear scale
// The input files provide log10(σ), so we convert: σ = 10^(log10(σ))
// ====================================================================

double fXS150Op(double x) { 
    return TMath::Power(10, XSLog150MeV(x)); 
}

double fXS160Op(double x) { 
    return TMath::Power(10, XSLog160MeV(x)); 
}

double fXS170Op(double x) { 
    return TMath::Power(10, XSLog170MeV(x)); 
}

double fXS180Op(double x) { 
    return TMath::Power(10, XSLog180MeV(x)); 
}

double fXS190Op(double x) { 
    return TMath::Power(10, XSLog190MeV(x)); 
}

double fXS200Op(double x) { 
    return TMath::Power(10, XSLog200MeV(x)); 
}

double fXS210Op(double x) { 
    return TMath::Power(10, XSLog210MeV(x)); 
}

double fXS220Op(double x) { 
    return TMath::Power(10, XSLog220MeV(x)); 
}

double fXS230Op(double x) { 
    return TMath::Power(10, XSLog230MeV(x)); 
}

double fXS240Op(double x) { 
    return TMath::Power(10, XSLog240MeV(x)); 
}

// ====================================================================
// Helper function: Get cross section formula string for given energy
// Returns the appropriate function name based on beam kinetic energy
// ====================================================================
const char* GetXSFormula(Double_t ekin)
{
    if      (ekin == 150) return "fXS150Op(x)";
    else if (ekin == 160) return "fXS160Op(x)";
    else if (ekin == 170) return "fXS170Op(x)";
    else if (ekin == 180) return "fXS180Op(x)";
    else if (ekin == 190) return "fXS190Op(x)";
    else if (ekin == 200) return "fXS200Op(x)";
    else if (ekin == 210) return "fXS210Op(x)";
    else if (ekin == 220) return "fXS220Op(x)";
    else if (ekin == 230) return "fXS230Op(x)";
    else if (ekin == 240) return "fXS240Op(x)";
    else {
        std::cout << "Warning: Energy " << ekin << " not in optical model data." << std::endl;
        std::cout << "Using cos(theta) distribution instead." << std::endl;
        return "cos(x)";
    }
}


// ====================================================================
// GetElasticAnalyzingPower: Calculate A_N for elastic scattering
// ====================================================================
Double_t GetElasticAnalyzingPower(Double_t ekin, Double_t theta_lab_deg)
{
    Double_t dT = ekin - T0_AN;
    Double_t dTheta = theta_lab_deg - theta0_AN;
    
    Double_t AN = 1.0 - alpha_AN*dT*dT - beta_AN*dT*dTheta - gamma_AN*dTheta*dTheta;
    
    // Clamp to physical range [-1, 1]
    if (AN > 1.0) AN = 1.0;
    if (AN < -1.0) AN = -1.0;
    
    return AN;
}