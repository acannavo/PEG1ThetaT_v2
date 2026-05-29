{
    // Load ROOT libraries that PLUTO needs
    gSystem->Load("libEG");
    gSystem->Load("libPhysics");
    gSystem->Load("libMatrix");
    gSystem->Load("libMathCore");
    // Suppress VLA extension warning from Pluto headers
	gSystem->AddIncludePath("-Wno-vla-cxx-extension");

    // Include paths
    gROOT->ProcessLine(".include $PLUTOSYS/include");
    gROOT->ProcessLine(".include /home/acannavo/pluto_v6.02/PEG1ThetaT_v2/include");
    
    // Load PLUTO
    gSystem->Load("$PLUTOSYS/lib/libPluto.so");
    
    // Load modules FIRST (with +)
    gROOT->ProcessLine(".L src/DetectorConfig.C+");
    gROOT->ProcessLine(".L src/ElasticScattering.C+");
    gROOT->ProcessLine(".L src/InelasticScattering.C+");
    gROOT->ProcessLine(".L src/Kinematics.C+");
    gROOT->ProcessLine(".L src/AnalyzingPowerUtils.C+");  // ← NEW utility module
    gROOT->ProcessLine(".L src/MeyerScattering.C+");      // ← TEMPORARILY DISABLED FOR TESTING
    gROOT->ProcessLine(".L src/EventGenerator.C+");
	cout << "  \n- EventGenerator compiled here" << endl;
	
    // Load test files (no +)
    
	gROOT->ProcessLine(".L src/TestDetectorConfig.C");
	cout << "  \n- TestDetectorConfig loaded here" << endl;
	gROOT->ProcessLine(".L AnalysisAN.C+");
    gROOT->ProcessLine(".L ComputeNormalization.C");
	
    cout << "\n=== PEG Environment Loaded ===" << endl;
    cout << "Modules loaded:" << endl;
    cout << "  - DetectorConfig" << endl;
    cout << "  - ElasticScattering" << endl;
    cout << "  - InelasticScattering" << endl;
    cout << "  - Kinematics" << endl;
    cout << "  - AnalyzingPowerUtils" << endl;
    cout << "  - MeyerScattering (NEW)" << endl;
    cout << "  - EventGenerator" << endl;
	cout << "  - AnalysisAN.C " << endl;
	cout << "  - ComputeNormalization.C " << endl;
    cout << "\nTestDetectorConfig() is ready to use!" << endl;
    cout << "MeyerScattering_Info() is available for testing" << endl;
    
    // Load the main event generator macro orchestrator LAST
    gROOT->ProcessLine(".L sigma_pC_Elas_EventGenerator.C");
}