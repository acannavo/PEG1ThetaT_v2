// TestMeyer_LinearBins.C
// Test with LINEAR bins to see if sampling is correct

void TestMeyer_LinearBins()
{
    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "    Meyer Test - LINEAR BINS" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;
    
    // Use SMALLER range to test
    Double_t t_min_MeV2 = 6e3;
    Double_t t_max_MeV2 = 1e5;  // Smaller range first!
    
    Double_t t_min_GeV2 = -t_max_MeV2 / 1.0e6;
    Double_t t_max_GeV2 = -t_min_MeV2 / 1.0e6;
    
    Double_t ekin = 180.0;
    Int_t n_samples = 100000;
    
    cout << "Configuration:" << endl;
    cout << "  |t| range: " << t_min_MeV2 << " to " << t_max_MeV2 << " (MeV/c)²" << endl;
    cout << "  Energy: " << ekin << " MeV" << endl;
    cout << endl;
    
    // Create envelope
    TH1D* envelope = CreateMeyerEnvelope_Elastic(t_min_GeV2, t_max_GeV2);
    
    // Sample events in LINEAR bins
    TH1D* h_samples = new TH1D("h_samp", "Sampled", 100, t_min_MeV2, t_max_MeV2);
    TH1D* h_expected = new TH1D("h_exp", "Expected", 100, t_min_MeV2, t_max_MeV2);
    
    cout << "Sampling..." << endl;
    for (Int_t i = 0; i < n_samples; i++) {
        Double_t t_GeV2 = SampleT_Meyer_Elastic(ekin, envelope);
        Double_t t_MeV2 = TMath::Abs(t_GeV2) * 1e6;
        h_samples->Fill(t_MeV2);
    }
    cout << "Done!" << endl;
    cout << endl;
    
    // Fill expected
    for (Int_t i = 1; i <= 100; i++) {
        Double_t t_MeV2 = h_expected->GetBinCenter(i);
        Double_t t_GeV2 = -t_MeV2 / 1e6;
        Double_t xs = MeyerXS_Elastic(t_GeV2, ekin);
        h_expected->SetBinContent(i, xs);
    }
    
    // Scale to compare shapes
    h_samples->Scale(1.0 / h_samples->Integral());
    h_expected->Scale(1.0 / h_expected->Integral());
    
    cout << "Results:" << endl;
    cout << "  Sampled mean: " << h_samples->GetMean() << " (MeV/c)²" << endl;
    cout << "  Expected mean: " << h_expected->GetMean() << " (MeV/c)²" << endl;
    cout << "  XS at minimum |t|: " << MeyerXS_Elastic(t_max_GeV2, ekin) << " mb/GeV²" << endl;
    cout << "  XS at maximum |t|: " << MeyerXS_Elastic(t_min_GeV2, ekin) << " mb/GeV²" << endl;
    cout << endl;
    
    // Plot
    TCanvas* c1 = new TCanvas("c1", "Linear Bins Test", 800, 600);
    gPad->SetLogy();
    
    h_expected->SetLineColor(kRed);
    h_expected->SetLineWidth(2);
    h_expected->SetTitle("Sampling Test (Linear bins);|t| [(MeV/c)^{2}];Normalized d#sigma/dt");
    h_expected->Draw("HIST");
    
    h_samples->SetMarkerStyle(20);
    h_samples->SetMarkerColor(kBlue);
    h_samples->Draw("P SAME");
    
    TLegend* leg = new TLegend(0.6, 0.7, 0.85, 0.85);
    leg->AddEntry(h_expected, "Expected (XS)", "l");
    leg->AddEntry(h_samples, "Sampled events", "p");
    leg->Draw();
    
    TLatex* text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.03);
    text->DrawLatex(0.15, 0.85, Form("Range: %.0e - %.0e (MeV/c)^{2}", t_min_MeV2, t_max_MeV2));
    text->DrawLatex(0.15, 0.80, "Should match if sampling is correct");
    
    c1->Update();
    
    delete envelope;
}