// ====================================================================
// TestDetectorConfig.C
// ====================================================================
// Test script for DetectorConfig module
//
// Purpose: Verify that detector configuration can be changed dynamically
//          and that all derived constants are calculated correctly
//
// Usage:
//   root -l
//   .L DetectorConfig.C
//   .L TestDetectorConfig.C
//   TestDetectorConfig()
//
// Expected behavior:
//   1. Show default configuration (16.2°)
//   2. Change to different angles
//   3. Verify derived constants are correct
//   4. Reset to defaults
// ====================================================================

void TestDetectorConfig() {
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "  Testing DetectorConfig Module" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // ================================================================
    // Test 1: Show default configuration
    // ================================================================
    std::cout << "\n>>> Test 1: Default Configuration <<<" << std::endl;
    PrintDetectorConfig();
    
    // Verify default values
    std::cout << "Verification:" << std::endl;
    std::cout << "  Expected theta center: 16.2 deg" << std::endl;
    std::cout << "  Actual theta center:   " << DETECTOR_THETA_CENTER << " deg" << std::endl;
    
    if (TMath::Abs(DETECTOR_THETA_CENTER - 16.2) < 0.001) {
        std::cout << "  ✓ Default configuration correct!" << std::endl;
    } else {
        std::cout << "  ✗ ERROR: Default configuration wrong!" << std::endl;
    }
    
    // ================================================================
    // Test 2: Change to 14.0 degrees
    // ================================================================
    std::cout << "\n>>> Test 2: Set to 14.0 degrees <<<" << std::endl;
    SetDetectorConfig(14.0);
    
    // Verify change
    std::cout << "Verification:" << std::endl;
    std::cout << "  Expected theta center: 14.0 deg" << std::endl;
    std::cout << "  Actual theta center:   " << DETECTOR_THETA_CENTER << " deg" << std::endl;
    
    if (TMath::Abs(DETECTOR_THETA_CENTER - 14.0) < 0.001) {
        std::cout << "  ✓ Configuration changed correctly!" << std::endl;
    } else {
        std::cout << "  ✗ ERROR: Configuration not changed!" << std::endl;
    }
    
    // Verify derived constants updated
    Double_t expected_rad = 14.0 * TMath::DegToRad();
    if (TMath::Abs(DETECTOR_THETA_CENTER_RAD - expected_rad) < 0.0001) {
        std::cout << "  ✓ Derived constants updated correctly!" << std::endl;
    } else {
        std::cout << "  ✗ ERROR: Derived constants not updated!" << std::endl;
    }
    
    // ================================================================
    // Test 3: Change to 18.0 degrees with custom window
    // ================================================================
    std::cout << "\n>>> Test 3: Set to 18.0 degrees with custom window <<<" << std::endl;
    SetDetectorConfig(18.0, 0.050);  // Smaller window
    
    std::cout << "Verification:" << std::endl;
    if (TMath::Abs(DETECTOR_THETA_CENTER - 18.0) < 0.001 &&
        TMath::Abs(DETECTOR_THETA_WINDOW - 0.050) < 0.0001) {
        std::cout << "  ✓ Custom configuration set correctly!" << std::endl;
    } else {
        std::cout << "  ✗ ERROR: Custom configuration wrong!" << std::endl;
    }
    
    // ================================================================
    // Test 4: Test a sequence (simulating systematic study)
    // ================================================================
    std::cout << "\n>>> Test 4: Simulate systematic scan <<<" << std::endl;
    
    Double_t angles[] = {14.0, 15.0, 16.0, 16.2, 17.0, 18.0};
    Int_t nangles = 6;
    
    std::cout << "\nScanning through detector positions:" << std::endl;
    std::cout << "  Angle [deg]  |  Center [rad]  |  Min [rad]  |  Max [rad]" << std::endl;
    std::cout << "  -------------------------------------------------------------" << std::endl;
    
    for (Int_t i = 0; i < nangles; i++) {
        SetDetectorConfig(angles[i]);
        std::cout << "     " << Form("%5.1f", angles[i]);
        std::cout << "      |    " << Form("%6.4f", DETECTOR_THETA_CENTER_RAD);
        std::cout << "      |   " << Form("%6.4f", DETECTOR_THETA_MIN);
        std::cout << "    |   " << Form("%6.4f", DETECTOR_THETA_MAX) << std::endl;
    }
    
    std::cout << "\n  ✓ Systematic scan completed!" << std::endl;
    
    // ================================================================
    // Test 5: Reset to defaults
    // ================================================================
    std::cout << "\n>>> Test 5: Reset to defaults <<<" << std::endl;
    ResetDetectorConfig();
    
    if (TMath::Abs(DETECTOR_THETA_CENTER - 16.2) < 0.001) {
        std::cout << "  ✓ Successfully reset to defaults!" << std::endl;
    } else {
        std::cout << "  ✗ ERROR: Reset failed!" << std::endl;
    }
    
    // ================================================================
    // Summary
    // ================================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "  All Tests Completed!" << std::endl;
    std::cout << "========================================" << std::endl;

}