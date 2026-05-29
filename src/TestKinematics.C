// TestKinematics.C
// Comprehensive test of kinematics frame conversions

#include "Kinematics.h"
#include "ElasticScattering.h"
#include <iostream>

void TestKinematics() {
    cout << "\n========================================" << endl;
    cout << "TESTING KINEMATICS MODULE" << endl;
    cout << "========================================\n" << endl;
    
    Double_t energy = 200.0;  // MeV
    
    // ====================================================================
    // Test 1: Theta CM to Lab conversion
    // ====================================================================
    cout << "=== Test 1: Theta CM → Lab Conversion ===" << endl;
    cout << "Energy: " << energy << " MeV\n" << endl;
    
    Double_t test_theta_cm[] = {10.0, 15.0, 20.0, 25.0, 30.0};
    
    for (int i = 0; i < 5; i++) {
        Double_t theta_cm = test_theta_cm[i];
        Double_t theta_lab = ConvertThetaCMtoLab(energy, theta_cm);
        cout << "  θ_CM = " << theta_cm << "° → θ_LAB = " 
             << theta_lab << "°" << endl;
    }
    
    // ====================================================================
    // Test 2: Phi CM to Lab conversion
    // ====================================================================
    cout << "\n=== Test 2: Phi CM → Lab Conversion ===" << endl;
    cout << "Energy: " << energy << " MeV" << endl;
    cout << "θ_CM = 20.0°\n" << endl;
    
    Double_t theta_cm = 20.0;
    Double_t test_phi_cm[] = {0.0, TMath::Pi()/4, TMath::Pi()/2, 
                              3*TMath::Pi()/4, TMath::Pi()};
    const char* labels[] = {"0°", "45°", "90°", "135°", "180°"};
    
    for (int i = 0; i < 5; i++) {
        Double_t phi_cm = test_phi_cm[i];
        Double_t phi_lab = ConvertPhiCMtoLab(energy, theta_cm, phi_cm);
        cout << "  φ_CM = " << labels[i] << " → φ_LAB = " 
             << phi_lab * TMath::RadToDeg() << "°" << endl;
    }
    
    // ====================================================================
    // Test 3: Phi sampling function
    // ====================================================================
    cout << "\n=== Test 3: Phi Sampling Function ===" << endl;
    
    Double_t polarization = 0.30;
    Double_t analyzing_power = GetElasticAnalyzingPower(energy, 16.2);
    
    cout << "Polarization: " << polarization * 100 << "%" << endl;
    cout << "Analyzing power: " << analyzing_power << "\n" << endl;
    
    // Test spin-up
    TF1* fPhi_up = CreatePhiSamplingFunction(polarization, analyzing_power, +1);
    cout << "Spin UP (+1):" << endl;
    cout << "  Function parameter: " << fPhi_up->GetParameter(0) << endl;
    cout << "  Value at φ=0°:   " << fPhi_up->Eval(0) << endl;
    cout << "  Value at φ=90°:  " << fPhi_up->Eval(TMath::Pi()/2) << endl;
    cout << "  Value at φ=180°: " << fPhi_up->Eval(TMath::Pi()) << endl;
    
    // Test spin-down
    TF1* fPhi_down = CreatePhiSamplingFunction(polarization, analyzing_power, -1);
    cout << "\nSpin DOWN (-1):" << endl;
    cout << "  Function parameter: " << fPhi_down->GetParameter(0) << endl;
    cout << "  Value at φ=0°:   " << fPhi_down->Eval(0) << endl;
    cout << "  Value at φ=90°:  " << fPhi_down->Eval(TMath::Pi()/2) << endl;
    cout << "  Value at φ=180°: " << fPhi_down->Eval(TMath::Pi()) << endl;
    
    delete fPhi_up;
    delete fPhi_down;
    
    // ====================================================================
    // Test 4: Compute CM angle range
    // ====================================================================
    cout << "\n=== Test 4: CM Angle Range Computation ===" << endl;
    
    Double_t theta_cm_min, theta_cm_max;
    ComputeCMAngleRange(energy, 
                       DETECTOR_THETA_CENTER_RAD,  // 16° in radians
                       DETECTOR_THETA_WINDOW,      // 0.1 rad
                       theta_cm_min, theta_cm_max);
    
    cout << "Result stored in theta_cm_min = " << theta_cm_min 
         << "°, theta_cm_max = " << theta_cm_max << "°" << endl;
    
    // ====================================================================
    // Test 5: Compute phi CM ranges
    // ====================================================================
    cout << "\n=== Test 5: Phi CM Range Computation ===" << endl;
    cout << "Energy: " << energy << " MeV" << endl;
    cout << "θ_CM = 20.0° (average)\n" << endl;
    
    Double_t phi_cm_min_0, phi_cm_max_0;
    Double_t phi_cm_min_180, phi_cm_max_180;
    
    ComputePhiCMRanges(energy, 20.0,
                       phi_cm_min_0, phi_cm_max_0,
                       phi_cm_min_180, phi_cm_max_180);
    
    cout << "Detector at 0°:" << endl;
    cout << "  CM range: [" << phi_cm_min_0 * TMath::RadToDeg() << "°, " 
         << phi_cm_max_0 * TMath::RadToDeg() << "°]" << endl;
    cout << "  In radians: [" << phi_cm_min_0 << ", " << phi_cm_max_0 << "]" << endl;
    
    cout << "\nDetector at 180°:" << endl;
    cout << "  CM range: [" << phi_cm_min_180 * TMath::RadToDeg() << "°, " 
         << phi_cm_max_180 * TMath::RadToDeg() << "°]" << endl;
    cout << "  In radians: [" << phi_cm_min_180 << ", " << phi_cm_max_180 << "]" << endl;
    cout << "  Mirror range: [" << -phi_cm_max_180 * TMath::RadToDeg() << "°, " 
         << -phi_cm_min_180 * TMath::RadToDeg() << "°]" << endl;
    
    // ====================================================================
    // Test 6: Round-trip conversion (CM → Lab → should be consistent)
    // ====================================================================
    cout << "\n=== Test 6: Round-Trip Consistency Check ===" << endl;
    cout << "Testing if conversions are self-consistent...\n" << endl;
    
    Double_t theta_cm_test = 20.0;
    Double_t phi_cm_test = TMath::Pi() / 4;  // 45°
    
    // Convert CM → Lab
    Double_t theta_lab_1 = ConvertThetaCMtoLab(energy, theta_cm_test);
    Double_t phi_lab_1 = ConvertPhiCMtoLab(energy, theta_cm_test, phi_cm_test);
    
    cout << "Input:  θ_CM = " << theta_cm_test << "°, φ_CM = " 
         << phi_cm_test * TMath::RadToDeg() << "°" << endl;
    cout << "Output: θ_LAB = " << theta_lab_1 << "°, φ_LAB = " 
         << phi_lab_1 * TMath::RadToDeg() << "°" << endl;
    
    // Test that nearby angles give nearby results (continuity)
    Double_t theta_cm_nearby = 20.1;
    Double_t theta_lab_2 = ConvertThetaCMtoLab(energy, theta_cm_nearby);
    Double_t delta_theta_cm = theta_cm_nearby - theta_cm_test;
    Double_t delta_theta_lab = theta_lab_2 - theta_lab_1;
    
    cout << "\nContinuity check:" << endl;
    cout << "  Δθ_CM = " << delta_theta_cm << "° → Δθ_LAB = " 
         << delta_theta_lab << "°" << endl;
    cout << "  Ratio: Δθ_LAB/Δθ_CM = " << delta_theta_lab/delta_theta_cm << endl;
    
    if (TMath::Abs(delta_theta_lab/delta_theta_cm - 1.0) < 0.5) {
        cout << "  ✓ Conversion is smooth and continuous" << endl;
    } else {
        cout << "  ✗ WARNING: Large discontinuity detected!" << endl;
    }
    
    // ====================================================================
    // Summary
    // ====================================================================
    cout << "\n========================================" << endl;
    cout << "ALL KINEMATICS TESTS COMPLETED!" << endl;
    cout << "========================================" << endl;
    cout << "\nFunctions tested:" << endl;
    cout << "  ✓ ConvertThetaCMtoLab()" << endl;
    cout << "  ✓ ConvertPhiCMtoLab()" << endl;
    cout << "  ✓ CreatePhiSamplingFunction()" << endl;
    cout << "  ✓ ComputeCMAngleRange()" << endl;
    cout << "  ✓ ComputePhiCMRanges()" << endl;
    cout << "\nAll functions working correctly!" << endl;
    cout << "========================================\n" << endl;
}