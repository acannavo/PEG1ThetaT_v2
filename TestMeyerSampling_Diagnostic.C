// TestMeyerSampling_Diagnostic.C
// Diagnostic version to understand the sampling

void TestMeyerSampling_Diagnostic()
{
    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "    Meyer Sampling Diagnostic" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;
    
    // Test with SMALLER range first to understand
    Double_t t_min_MeV2 = 6e3;   // (MeV/c)²
    Double_t t_max_MeV2 = 6e5;   // (MeV/c)²
    
    Double_t t_min_GeV2 = -t_max_MeV2 / 1.0e6;
    Double_t t_max_GeV2 = -t_min_MeV2 / 1.0e6;
    
    Double_t ekin = 180.0;
    
    cout << "Creating envelope..." << endl;
    TH1D* envelope = CreateMeyerEnvelope_Elastic(t_min_GeV2, t_max_GeV2);
    
    // Check envelope properties
    cout << "\nEnvelope properties:" << endl;
    cout << "  Total bins: " << envelope->GetNbinsX() << endl;
    cout << "  Bin width: " << envelope->GetBinWidth(1) << " GeV²" << endl;
    
    // Find where XS is maximum
    Double_t max_xs = -1;
    Double_t t_at_max = 0;
    for (int i = 1; i <= envelope->GetNbinsX(); i++) {
        Double_t xs = envelope->GetBinContent(i);
        Double_t t = envelope->GetBinCenter(i);
        if (xs > max_xs) {
            max_xs = xs;
            t_at_max = t;
        }
    }
    
    cout << "  Maximum XS = " << max_xs << " mb/GeV²" << endl;
    cout << "  at t = " << t_at_max << " GeV²" << endl;
    cout << "  at |t| = " << TMath::Abs(t_at_max)*1e6 << " (MeV/c)²" << endl;
    
    // Sample a few events
    cout << "\nSampling 20 test events:" << endl;
    for (int i = 0; i < 20; i++) {
        Double_t t = SampleT_Meyer_Elastic(ekin, envelope);
        Double_t t_MeV2 = TMath::Abs(t) * 1e6;
        cout << "  Event " << i+1 << ": |t| = " << t_MeV2 << " (MeV/c)²" << endl;
    }
    
    // Now check XS values across range
    cout << "\nCross section at different |t|:" << endl;
    Double_t test_t_values[] = {6e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5};
    for (int i = 0; i < 7; i++) {
        Double_t t_MeV2 = test_t_values[i];
        Double_t t_GeV2 = -t_MeV2 / 1e6;
        Double_t xs = MeyerXS_Elastic(t_GeV2, ekin);
        cout << "  |t| = " << t_MeV2 << " (MeV/c)²  →  XS = " << xs << " mb/GeV²" << endl;
    }
    
    delete envelope;
}