// TestKinematicsConversions.C
// Unit tests for t ↔ θ_CM conversion functions

#include "Kinematics.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void TestKinematicsConversions()
{
    cout << "\n╔════════════════════════════════════════════════════════╗" << endl;
    cout << "║  Testing t ↔ θ_CM Conversion Functions                ║" << endl;
    cout << "╚════════════════════════════════════════════════════════╝\n" << endl;
    
    Double_t energy = 200.0;  // MeV
    
    // Test 1: Calculate CM momentum
    cout << "Test 1: CM Momentum Calculation" << endl;
    cout << "────────────────────────────────────────────────────────" << endl;
    Double_t pcm = CalculateCMMomentum(energy);
    cout << "  Beam energy: " << energy << " MeV" << endl;
    cout << "  CM momentum: " << pcm << " MeV/c" << endl;
    cout << "  Expected:    ~585.67 MeV/c" << endl;
    cout << "  Match:       " << (abs(pcm - 585.67) < 1.0 ? "✓ PASS" : "✗ FAIL") << endl;
    cout << endl;
    
    // Test 2: Round-trip conversion θ → t → θ
    cout << "Test 2: Round-Trip Conversion (θ → t → θ)" << endl;
    cout << "────────────────────────────────────────────────────────" << endl;
    cout << setw(12) << "θ_CM (deg)" 
         << setw(15) << "t (GeV²)" 
         << setw(15) << "θ_back (deg)" 
         << setw(15) << "Error (deg)" 
         << setw(10) << "Status" << endl;
    cout << "────────────────────────────────────────────────────────" << endl;
    
    Bool_t all_passed = true;
    for (Double_t theta = 0.0; theta <= 50.0; theta += 5.0) {
        Double_t t = ConvertThetaCMtoT(theta, energy);
        Double_t theta_back = ConvertTtoThetaCM(t, energy);
        Double_t error = abs(theta - theta_back);
        Bool_t passed = (error < 1.0e-6);  // 1 micro-degree tolerance
        
        cout << fixed << setprecision(6);
        cout << setw(12) << theta 
             << setw(15) << t 
             << setw(15) << theta_back 
             << setw(15) << error 
             << setw(10) << (passed ? "✓" : "✗") << endl;
        
        if (!passed) all_passed = false;
    }
    cout << "────────────────────────────────────────────────────────" << endl;
    cout << "Result: " << (all_passed ? "✓ ALL TESTS PASSED" : "✗ SOME TESTS FAILED") << endl;
    cout << endl;
    
    // Test 3: Round-trip conversion t → θ → t
    cout << "Test 3: Round-Trip Conversion (t → θ → t)" << endl;
    cout << "────────────────────────────────────────────────────────" << endl;
    cout << setw(15) << "t (GeV²)" 
         << setw(15) << "θ_CM (deg)" 
         << setw(15) << "t_back (GeV²)" 
         << setw(15) << "Error (GeV²)" 
         << setw(10) << "Status" << endl;
    cout << "────────────────────────────────────────────────────────" << endl;
    
    all_passed = true;
    for (Double_t t = -0.1; t <= 0.0; t += 0.01) {
        Double_t theta = ConvertTtoThetaCM(t, energy);
        Double_t t_back = ConvertThetaCMtoT(theta, energy);
        Double_t error = abs(t - t_back);
        Bool_t passed = (error < 1.0e-9);  // nano-GeV² tolerance
        
        cout << fixed << setprecision(6);
        cout << setw(15) << t 
             << setw(15) << theta 
             << setw(15) << t_back 
             << setw(15) << scientific << error << fixed
             << setw(10) << (passed ? "✓" : "✗") << endl;
        
        if (!passed) all_passed = false;
    }
    cout << "────────────────────────────────────────────────────────" << endl;
    cout << "Result: " << (all_passed ? "✓ ALL TESTS PASSED" : "✗ SOME TESTS FAILED") << endl;
    cout << endl;
    
    // Test 4: Physical values at detector acceptance
    cout << "Test 4: Physical Values at Detector Acceptance" << endl;
    cout << "────────────────────────────────────────────────────────" << endl;
    cout << "For θ_Lab = 16.2° detector:" << endl;
    
    // Approximate CM angles from your Monte Carlo
    Double_t theta_cm_min = 19.4;  // degrees
    Double_t theta_cm_max = 20.6;  // degrees
    
    Double_t t_at_min = ConvertThetaCMtoT(theta_cm_min, energy);
    Double_t t_at_max = ConvertThetaCMtoT(theta_cm_max, energy);
    
    cout << "  θ_CM range:     [" << theta_cm_min << "°, " << theta_cm_max << "°]" << endl;
    cout << "  t range:        [" << t_at_max << ", " << t_at_min << "] GeV²" << endl;
    cout << "  Expected:       ~[-0.0131, -0.0128] GeV²" << endl;
    cout << "  Δt:             " << (t_at_min - t_at_max) << " GeV²" << endl;
    cout << endl;
    
    // Test 5: Boundary conditions
    cout << "Test 5: Boundary Conditions" << endl;
    cout << "────────────────────────────────────────────────────────" << endl;
    
    // Forward scattering: θ = 0° → t = 0
    Double_t t_forward = ConvertThetaCMtoT(0.0, energy);
    cout << "  θ_CM = 0°:   t = " << t_forward << " GeV²";
    cout << "  (expected: 0.0)" << endl;
    
    // Backward scattering: θ = 180° → t = -4p²
    Double_t t_backward = ConvertThetaCMtoT(180.0, energy);
    Double_t t_backward_expected = -4.0 * (pcm/1000.0) * (pcm/1000.0);
    cout << "  θ_CM = 180°: t = " << t_backward << " GeV²";
    cout << "  (expected: " << t_backward_expected << ")" << endl;
    
    cout << "────────────────────────────────────────────────────────" << endl;
    cout << endl;
    
    cout << "╔════════════════════════════════════════════════════════╗" << endl;
    cout << "║  All Conversion Tests Complete                         ║" << endl;
    cout << "╚════════════════════════════════════════════════════════╝\n" << endl;
}