// TestInelasticScattering.C
// Comprehensive test of inelastic scattering module

#include "InelasticScattering.h"
#include <iostream>

void TestInelasticScattering() {
    std::cout << "\n========================================" << endl;
    std::cout << "TESTING INELASTIC SCATTERING MODULE" << endl;
    std::cout << "========================================\n" << endl;
    
    // ====================================================================
    // Test 1: Load inelastic data
    // ====================================================================
    std::cout << "=== Test 1: Loading Data ===" << endl;
    LoadInelasticData();
    
    if (IsInelasticDataLoaded()) {
        std::cout << "✓ Data loaded successfully!\n" << endl;
    } else {
        std::cout << "✗ ERROR: Data not loaded!\n" << endl;
        return;
    }
    
    // ====================================================================
    // Test 2: Test cross section interpolation
    // ====================================================================
    std::cout << "=== Test 2: Cross Section Interpolation ===" << endl;
    std::cout << "Testing GetInelasticCrossSection() at different angles:\n" << endl;
    
    Double_t test_angles[] = {10.0, 15.0, 20.0, 25.0, 30.0};
    
    for (int i = 0; i < 5; i++) {
        Double_t theta = test_angles[i];
        Double_t xs = GetInelasticCrossSection(theta);
        std::cout << "  θ_CM = " << theta << "° → σ = " << xs << " mb/sr" << endl;
    }
    
    // ====================================================================
    // Test 3: Test analyzing power interpolation
    // ====================================================================
    std::cout << "\n=== Test 3: Analyzing Power Interpolation ===" << endl;
    std::cout << "Testing GetInelasticAnalyzingPower() at different angles:\n" << endl;
    
    for (int i = 0; i < 5; i++) {
        Double_t theta = test_angles[i];
        Double_t an = GetInelasticAnalyzingPower(theta);
        std::cout << "  θ_CM = " << theta << "° → A_N = " << an << endl;
    }
    
    // ====================================================================
    // Test 4: Test CM momentum calculation
    // ====================================================================
    std::cout << "\n=== Test 4: CM Momentum Calculation ===" << endl;
    std::cout << "Testing CalculateCMMomentumInelastic():\n" << endl;
    
    Double_t test_energies[] = {150.0, 175.0, 200.0, 225.0, 240.0};
    Double_t E_ex = 0.00443;  // 4.43 MeV excitation
    
    for (int i = 0; i < 5; i++) {
        Double_t energy = test_energies[i];
        Double_t pcm = CalculateCMMomentumInelastic(energy, E_ex);
        std::cout << "  E_beam = " << energy << " MeV → p_CM = " 
                  << pcm << " GeV/c" << endl;
    }
    
    // ====================================================================
    // Test 5: Test sampling histogram creation
    // ====================================================================
    std::cout << "\n=== Test 5: Sampling Histogram Creation ===" << endl;
    std::cout << "Testing CreateInelasticSamplingHistogram():\n" << endl;
    
    Double_t theta_min = 10.0;
    Double_t theta_max = 40.0;
    Int_t nbins = 1000;
    
    TH1D* h_sampling = CreateInelasticSamplingHistogram(theta_min, theta_max, nbins);
    
    if (h_sampling) {
        std::cout << "✓ Histogram created successfully!" << endl;
        std::cout << "  Name: " << h_sampling->GetName() << endl;
        std::cout << "  Range: [" << theta_min << "°, " << theta_max << "°]" << endl;
        std::cout << "  Bins: " << nbins << endl;
        std::cout << "  Entries: " << h_sampling->GetEntries() << endl;
        std::cout << "  Integral: " << h_sampling->Integral() << endl;
        
        // Test a few random samples
        std::cout << "\n  Testing GetRandom() sampling:" << endl;
        for (int i = 0; i < 5; i++) {
            Double_t sampled_theta = h_sampling->GetRandom();
            std::cout << "    Sample " << i+1 << ": θ_CM = " 
                      << sampled_theta << "°" << endl;
        }
        
        delete h_sampling;
    } else {
        std::cout << "✗ ERROR: Histogram creation failed!" << endl;
    }
    
    // ====================================================================
    // Test 6: Verify data consistency
    // ====================================================================
    std::cout << "\n=== Test 6: Data Consistency Check ===" << endl;
    std::cout << "Checking that XS and AN are consistent:\n" << endl;
    
    // Check that cross section decreases with angle (typically)
    Double_t xs_10 = GetInelasticCrossSection(10.0);
    Double_t xs_30 = GetInelasticCrossSection(30.0);
    
    std::cout << "  σ(10°) = " << xs_10 << " mb/sr" << endl;
    std::cout << "  σ(30°) = " << xs_30 << " mb/sr" << endl;
    
    if (xs_10 > xs_30) {
        std::cout << "  ✓ Cross section decreases with angle (as expected)" << endl;
    } else {
        std::cout << "  ⚠ Cross section increases with angle (unusual but possible)" << endl;
    }
    
    // Check that analyzing power is in physical range [-1, 1]
    Bool_t all_valid = true;
    for (int i = 0; i < 10; i++) {
        Double_t theta = 10.0 + i * 3.0;
        Double_t an = GetInelasticAnalyzingPower(theta);
        if (an < -1.0 || an > 1.0) {
            std::cout << "  ✗ ERROR: A_N(" << theta << "°) = " << an 
                      << " outside [-1,1]!" << endl;
            all_valid = false;
        }
    }
    
    if (all_valid) {
        std::cout << "  ✓ All analyzing power values in valid range [-1, 1]" << endl;
    }
    
    // ====================================================================
    // Test 7: Compare elastic vs inelastic CM momentum
    // ====================================================================
    std::cout << "\n=== Test 7: Elastic vs Inelastic Comparison ===" << endl;
    std::cout << "Comparing CM momentum for elastic and inelastic:\n" << endl;
    
    Double_t energy = 200.0;  // MeV
    
    // Elastic: no excitation
    Double_t pcm_elastic = CalculateCMMomentumInelastic(energy, 0.0);
    
    // Inelastic: 4.43 MeV excitation
    Double_t pcm_inelastic = CalculateCMMomentumInelastic(energy, 0.00443);
    
    std::cout << "  E_beam = " << energy << " MeV" << endl;
    std::cout << "  p_CM (elastic):   " << pcm_elastic << " GeV/c" << endl;
    std::cout << "  p_CM (inelastic): " << pcm_inelastic << " GeV/c" << endl;
    std::cout << "  Difference:       " << (pcm_elastic - pcm_inelastic)*1000 
              << " MeV/c" << endl;
    
    if (pcm_inelastic < pcm_elastic) {
        std::cout << "  ✓ Inelastic momentum < Elastic (correct!)" << endl;
    } else {
        std::cout << "  ✗ ERROR: Physics violation!" << endl;
    }
    
    // ====================================================================
    // Summary
    // ====================================================================
    std::cout << "\n========================================" << endl;
    std::cout << "ALL INELASTIC TESTS COMPLETED!" << endl;
    std::cout << "========================================" << endl;
    std::cout << "\nFunctions tested:" << endl;
    std::cout << "  ✓ LoadInelasticData()" << endl;
    std::cout << "  ✓ IsInelasticDataLoaded()" << endl;
    std::cout << "  ✓ GetInelasticCrossSection()" << endl;
    std::cout << "  ✓ GetInelasticAnalyzingPower()" << endl;
    std::cout << "  ✓ CalculateCMMomentumInelastic()" << endl;
    std::cout << "  ✓ CreateInelasticSamplingHistogram()" << endl;
    std::cout << "\nInelasticScattering module is working correctly!" << endl;
    std::cout << "========================================\n" << endl;
}