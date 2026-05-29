// TestMeyerSampling_Optimized.C
// Fixed version that doesn't clone 50k bin histograms

void TestMeyerSampling_Optimized()
{
    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "    Meyer Envelope Sampling Test" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;
    
    // Test parameters
    Double_t t_min = -0.15;  // GeV²
    Double_t t_max = -0.01;  // GeV²
    Double_t ekin = 180.0;   // MeV (midpoint between 160-200)
    Int_t n_samples = 10000; // Number of samples to test
    
    cout << "Test configuration:" << endl;
    cout << "  t range: " << t_min << " to " << t_max << " GeV²" << endl;
    cout << "  Energy: " << ekin << " MeV" << endl;
    cout << "  Samples: " << n_samples << endl;
    cout << endl;
    
    // ================================================================
    // 1. Create envelope
    // ================================================================
    cout << "Creating Meyer envelope (elastic)..." << endl;
    TH1D* envelope = CreateMeyerEnvelope_Elastic(t_min, t_max);  // 50k bins
    cout << "  Envelope created with " << envelope->GetNbinsX() << " bins" << endl;
    cout << endl;
    
    // ================================================================
    // 2. Sample t values and collect statistics
    // ================================================================
    cout << "Sampling " << n_samples << " events..." << endl;
    TH1D* h_samples = new TH1D("h_samples", "Sampled t distribution", 100, t_min, t_max);
    
    TStopwatch timer;
    timer.Start();
    
    for (Int_t i = 0; i < n_samples; i++) {
        Double_t t = SampleT_Meyer_Elastic(ekin, envelope);
        h_samples->Fill(t);
        
        if (i == 0) {
            cout << "  First sampled t = " << t << " GeV²" << endl;
        }
    }
    
    timer.Stop();
    
    cout << "  Done!" << endl;
    cout << "  Time: " << timer.RealTime() << " seconds" << endl;
    cout << "  Rate: " << n_samples / timer.RealTime() << " samples/sec" << endl;
    cout << endl;
    
    // ================================================================
    // 3. Compare sampled distribution to expected XS
    // ================================================================
    cout << "Comparing sampled distribution to expected XS..." << endl;
    
    // Create histogram with expected distribution
    TH1D* h_expected = new TH1D("h_expected", "Expected XS", 100, t_min, t_max);
    for (Int_t ibin = 1; ibin <= 100; ibin++) {
        Double_t t = h_expected->GetBinCenter(ibin);
        Double_t xs = MeyerXS_Elastic(t, ekin);
        h_expected->SetBinContent(ibin, xs);
    }
    
    // Normalize both to same area for comparison
    h_samples->Scale(1.0 / h_samples->Integral());
    h_expected->Scale(1.0 / h_expected->Integral());
    
    cout << "  Sampled mean t: " << h_samples->GetMean() << " GeV²" << endl;
    cout << "  Expected mean t: " << h_expected->GetMean() << " GeV²" << endl;
    cout << "  Difference: " << TMath::Abs(h_samples->GetMean() - h_expected->GetMean()) << " GeV²" << endl;
    cout << endl;
    
    // ================================================================
    // 4. Plot results (using FRESH histograms, not clones!)
    // ================================================================
    cout << "Creating comparison plot..." << endl;
    
    TCanvas* c1 = new TCanvas("c1", "Meyer Sampling Test", 1200, 600);
    c1->Divide(2, 1);
    
    // Left: Envelope vs actual XS
    c1->cd(1);
    gPad->SetLogy();
    gPad->SetGrid();
    
    // Create FRESH 100-bin histograms for plotting (not clones of 50k!)
    TH1D* h_xs_env = new TH1D("h_xs_env", "Envelope", 100, t_min, t_max);
    TH1D* h_xs_160 = new TH1D("h_xs_160", "XS 160", 100, t_min, t_max);
    TH1D* h_xs_200 = new TH1D("h_xs_200", "XS 200", 100, t_min, t_max);
    TH1D* h_xs_actual = new TH1D("h_xs_actual", "XS actual", 100, t_min, t_max);
    
    // Fill with XS values (only 100 bins, fast!)
    for (Int_t ibin = 1; ibin <= 100; ibin++) {
        Double_t t = h_xs_env->GetBinCenter(ibin);
        
        // Get envelope value from the 50k bin histogram
        Int_t env_bin = envelope->FindBin(t);
        Double_t env_val = envelope->GetBinContent(env_bin);
        h_xs_env->SetBinContent(ibin, env_val);
        
        // Calculate XS at different energies
        h_xs_160->SetBinContent(ibin, MeyerXS_Elastic(t, 160.0));
        h_xs_200->SetBinContent(ibin, MeyerXS_Elastic(t, 200.0));
        h_xs_actual->SetBinContent(ibin, MeyerXS_Elastic(t, ekin));
    }
    
    // Draw plots
    h_xs_env->SetLineColor(kBlack);
    h_xs_env->SetLineWidth(2);
    h_xs_env->SetLineStyle(2);
    h_xs_env->SetTitle("Meyer Envelope Method;|t| [GeV^{2}];d#sigma/dt [mb/GeV^{2}]");
    h_xs_env->GetXaxis()->SetTitleSize(0.045);
    h_xs_env->GetYaxis()->SetTitleSize(0.045);
    h_xs_env->GetXaxis()->SetLabelSize(0.04);
    h_xs_env->GetYaxis()->SetLabelSize(0.04);
    h_xs_env->SetMinimum(0.01);
    h_xs_env->Draw("HIST");
    
    h_xs_160->SetLineColor(kBlue);
    h_xs_160->SetLineWidth(2);
    h_xs_160->Draw("HIST SAME");
    
    h_xs_200->SetLineColor(kRed);
    h_xs_200->SetLineWidth(2);
    h_xs_200->Draw("HIST SAME");
    
    h_xs_actual->SetLineColor(kGreen+2);
    h_xs_actual->SetLineWidth(3);
    h_xs_actual->Draw("HIST SAME");
    
    TLegend* leg1 = new TLegend(0.6, 0.6, 0.85, 0.85);
    leg1->AddEntry(h_xs_env, "Envelope (MAX)", "l");
    leg1->AddEntry(h_xs_160, "XS @ 160 MeV", "l");
    leg1->AddEntry(h_xs_200, "XS @ 200 MeV", "l");
    leg1->AddEntry(h_xs_actual, Form("XS @ %.0f MeV", ekin), "l");
    leg1->Draw();
    
    // Right: Sampled vs expected
    c1->cd(2);
    gPad->SetGrid();
    
    h_expected->SetLineColor(kRed);
    h_expected->SetLineWidth(2);
    h_expected->SetTitle("Sampling Validation;|t| [GeV^{2}];Normalized counts");
    h_expected->GetXaxis()->SetTitleSize(0.045);
    h_expected->GetYaxis()->SetTitleSize(0.045);
    h_expected->GetXaxis()->SetLabelSize(0.04);
    h_expected->GetYaxis()->SetLabelSize(0.04);
    h_expected->Draw("HIST");
    
    h_samples->SetMarkerStyle(20);
    h_samples->SetMarkerSize(0.8);
    h_samples->SetMarkerColor(kBlue);
    h_samples->Draw("P SAME");
    
    TLegend* leg2 = new TLegend(0.6, 0.7, 0.85, 0.85);
    leg2->AddEntry(h_expected, "Expected XS", "l");
    leg2->AddEntry(h_samples, "Sampled", "p");
    leg2->Draw();
    
    c1->Update();
    
    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "✓ Test complete! Check the plots." << endl;
    cout << "════════════════════════════════════════════════\n" << endl;
    
    // Cleanup
    delete envelope;
}