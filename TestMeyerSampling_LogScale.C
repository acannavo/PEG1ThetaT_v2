// TestMeyerSampling_LogScale_Corrected.C
// Proper display: dN/d(log|t|) for log-scale histograms

void TestMeyerSampling_LogScale_Corrected()
{
    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "    Meyer Sampling Test (Corrected)" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;
    
    Double_t t_min_MeV2 = 6e3;
    Double_t t_max_MeV2 = 6e5;
    
    Double_t t_min_GeV2 = -t_max_MeV2 / 1.0e6;
    Double_t t_max_GeV2 = -t_min_MeV2 / 1.0e6;
    
    Double_t ekin = 180.0;
    Int_t n_samples = 100000;
    
    cout << "Test configuration:" << endl;
    cout << "  |t| range: " << t_min_MeV2 << " to " << t_max_MeV2 << " (MeV/c)²" << endl;
    cout << "  Energy: " << ekin << " MeV" << endl;
    cout << "  Samples: " << n_samples << endl;
    cout << endl;
    
    // Create envelope
    cout << "Creating envelope..." << endl;
    TH1D* envelope = CreateMeyerEnvelope_Elastic(t_min_GeV2, t_max_GeV2);
    cout << "  Envelope: " << envelope->GetNbinsX() << " bins" << endl;
    cout << endl;
    
    // Sample events - store in LINEAR bins first to see true distribution
    cout << "Sampling events..." << endl;
    TH1D* h_samples_linear = new TH1D("h_lin", "Linear", 200, t_min_GeV2, t_max_GeV2);
    
    TStopwatch timer;
    timer.Start();
    
    for (Int_t i = 0; i < n_samples; i++) {
        Double_t t_GeV2 = SampleT_Meyer_Elastic(ekin, envelope);
        h_samples_linear->Fill(t_GeV2);
    }
    
    timer.Stop();
    cout << "  Done! Rate: " << n_samples/timer.RealTime() << " events/sec" << endl;
    cout << endl;
    
    // Create LOG-scale histograms for display
    const Int_t nbins_log = 100;
    Double_t log_edges[nbins_log + 1];
    Double_t log_min = TMath::Log10(t_min_MeV2);
    Double_t log_max = TMath::Log10(t_max_MeV2);
    
    for (Int_t i = 0; i <= nbins_log; i++) {
        Double_t log_t = log_min + i * (log_max - log_min) / nbins_log;
        log_edges[i] = TMath::Power(10, log_t);
    }
    
    // Fill log histogram from sampled events
    TH1D* h_samples_log = new TH1D("h_log", "Sampled", nbins_log, log_edges);
    TH1D* h_expected_log = new TH1D("h_exp", "Expected", nbins_log, log_edges);
    
    // Transfer from linear to log histogram
    for (Int_t i = 1; i <= h_samples_linear->GetNbinsX(); i++) {
        Double_t t_GeV2 = h_samples_linear->GetBinCenter(i);
        Double_t t_MeV2 = TMath::Abs(t_GeV2) * 1e6;
        Double_t counts = h_samples_linear->GetBinContent(i);
        h_samples_log->Fill(t_MeV2, counts);
    }
    
    // Fill expected: XS × bin_width
    for (Int_t i = 1; i <= nbins_log; i++) {
        Double_t t_MeV2 = h_expected_log->GetBinCenter(i);
        Double_t t_GeV2 = -t_MeV2 / 1e6;
        Double_t xs = MeyerXS_Elastic(t_GeV2, ekin);
        Double_t bin_width = h_expected_log->GetBinWidth(i) / 1e6;  // Convert to GeV²
        h_expected_log->SetBinContent(i, xs * bin_width);
    }
    
    // Normalize
    h_samples_log->Scale(1.0 / h_samples_log->Integral());
    h_expected_log->Scale(1.0 / h_expected_log->Integral());
    
    cout << "Statistics:" << endl;
    cout << "  Sampled mean |t|: " << h_samples_log->GetMean() << " (MeV/c)²" << endl;
    cout << "  Expected mean |t|: " << h_expected_log->GetMean() << " (MeV/c)²" << endl;
    cout << endl;
    
    // Plot
    TCanvas* c1 = new TCanvas("c1", "Meyer Sampling", 1400, 600);
    c1->Divide(2, 1);
    
    // Left: Linear scale to see where events really are
    c1->cd(1);
    gPad->SetLogy();
    h_samples_linear->SetLineColor(kBlue);
    h_samples_linear->SetTitle("Event Distribution (Linear t scale);t [GeV^{2}];Events / bin");
    h_samples_linear->Draw("HIST");
    
    TLatex* text1 = new TLatex();
    text1->SetNDC();
    text1->SetTextSize(0.04);
    text1->DrawLatex(0.15, 0.85, "Events should cluster at t ~ 0");
    text1->DrawLatex(0.15, 0.80, "(smallest |t|, highest XS)");
    
    // Right: Log scale comparison
    c1->cd(2);
    gPad->SetLogx();
    gPad->SetGrid();
    
    h_expected_log->SetLineColor(kRed);
    h_expected_log->SetLineWidth(2);
    h_expected_log->SetTitle("Sampling Validation (Log scale);|t| [(MeV/c)^{2}];Normalized counts");
    h_expected_log->GetXaxis()->SetMoreLogLabels();
    h_expected_log->Draw("HIST");
    
    h_samples_log->SetMarkerStyle(20);
    h_samples_log->SetMarkerColor(kBlue);
    h_samples_log->SetMarkerSize(0.8);
    h_samples_log->Draw("P SAME");
    
    TLegend* leg = new TLegend(0.6, 0.7, 0.85, 0.85);
    leg->AddEntry(h_expected_log, "Expected", "l");
    leg->AddEntry(h_samples_log, "Sampled", "p");
    leg->Draw();
    
    c1->Update();
    
    cout << "✓ Check LEFT plot: events should peak near t=0 (low |t|)" << endl;
    cout << "✓ Check RIGHT plot: log-scale should match" << endl;
    cout << endl;
    
    delete envelope;
}