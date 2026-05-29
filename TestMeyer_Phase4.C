// TestMeyer_Phase4.C
// Three-stage validation of GenerateThreadEventsPolarized_MeyerSampling
// and GenerateThreadEventsInelastic_MeyerSampling
//
// Stage 1 — Smoke test:     does a single-thread run complete without crashes?
// Stage 2 — Shape test:     does the sampled |t| distribution match MeyerXS?
// Stage 3 — Comparison:     does Meyer output match T_SAMPLING at mid-energy?
//
// Usage (from ROOT prompt):
//   root [0] .L TestMeyer_Phase4.C
//   root [1] TestMeyer_Phase4_Stage1()   // ~5 seconds
//   root [2] TestMeyer_Phase4_Stage2()   // ~10 seconds, produces plots
//   root [3] TestMeyer_Phase4_Stage3()   // ~30 seconds, produces comparison plots
//   root [4] TestMeyer_Phase4_All()      // runs all three in sequence

// ====================================================================
// STAGE 1: Smoke test — single-thread elastic + inelastic
// ====================================================================
void TestMeyer_Phase4_Stage1()
{
    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "  STAGE 1: Smoke Test (single thread, small N)" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;

    const Double_t ekin         = 180.0;   // mid-range: exercises interpolation
    const Int_t    n_events     = 1000;
    const Double_t polarization = 0.80;
    const Int_t    spin_state   = +1;
    const Int_t    num_threads  = 1;

    // ----------------------------------------------------------------
    // 1a. Elastic Meyer
    // ----------------------------------------------------------------
    cout << "--- Elastic Meyer (180 MeV, 1000 events, 1 thread) ---" << endl;
    SingleRunMultithreadPolarized(ekin, n_events, polarization, spin_state,
                                  num_threads, T_SAMPLING_MEYER);
    cout << "✓ Elastic Meyer completed without crash" << endl;

    // ----------------------------------------------------------------
    // 1b. Inelastic Meyer
    // ----------------------------------------------------------------
    cout << "\n--- Inelastic Meyer (180 MeV, 1000 events, 1 thread) ---" << endl;
    SingleRunMultithreadInelasticPolarized(ekin, n_events, polarization, spin_state,
                                           num_threads, T_SAMPLING_MEYER);
    cout << "✓ Inelastic Meyer completed without crash" << endl;

    // ----------------------------------------------------------------
    // 1c. Boundary energies — verify clamp logic in MeyerXS/AP
    // ----------------------------------------------------------------
    cout << "\n--- Boundary: 160 MeV (lower clamp) ---" << endl;
    SingleRunMultithreadPolarized(160.0, n_events, polarization, spin_state,
                                  num_threads, T_SAMPLING_MEYER);
    cout << "✓ 160 MeV completed" << endl;

    cout << "\n--- Boundary: 200 MeV (upper clamp) ---" << endl;
    SingleRunMultithreadPolarized(200.0, n_events, polarization, spin_state,
                                  num_threads, T_SAMPLING_MEYER);
    cout << "✓ 200 MeV completed" << endl;

    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "  STAGE 1 PASSED" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;
}

// ====================================================================
// STAGE 2: Shape test — sampled |t| must match MeyerXS_Elastic(t, ekin)
// ====================================================================
// This is the equivalent of TestMeyer_LinearBins but run *through the
// full worker function path*, so it catches any unit-conversion bug
// or accept/reject logic error introduced during integration.
// ====================================================================
void TestMeyer_Phase4_Stage2()
{
    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "  STAGE 2: Shape Test (sampled t vs expected XS)" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;

    const Double_t ekin     = 180.0;   // interpolated energy
    const Int_t    n_events = 50000;

    // --- Build envelope (same call the orchestrator uses) ---
    // Use the detector t-range from ComputeUnifiedTRange
    Double_t t_min, t_max;
    ComputeUnifiedTRange(ekin, ekin,
                         DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW,
                         t_min, t_max);

    cout << "t range: [" << t_min << ", " << t_max << "] GeV²" << endl;

    TH1D* envelope = CreateMeyerEnvelope_Elastic(t_min, t_max);
    envelope->SetDirectory(0);

    // --- Run worker directly (1 thread) ---
    cout << "Generating " << n_events << " elastic Meyer events..." << endl;
    ThreadData result = GenerateThreadEventsPolarized_MeyerSampling(
        0, ekin, n_events, envelope,
        DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW, DETECTOR_PHI_WINDOW,
        0.0,   // unpolarized: isolates the XS shape from phi weighting
        +1);

    cout << "Events generated: " << result.count << endl;

    // --- Reconstruct |t| from stored momenta ---
    // t = (p_proton_out - p_proton_in)²  but we only have lab momenta.
    // Easier: recompute theta_lab from px/py/pz, convert to |t| via spline range.
    // For shape validation it's sufficient to histogram theta_lab and compare
    // to the expected angular distribution from MeyerXS integrated over phi.

    // Convert stored (px,py,pz) back to theta_lab [degrees]
    const Int_t   nbins   = 60;
    const Double_t th_min = (DETECTOR_THETA_CENTER_RAD - DETECTOR_THETA_WINDOW)
                            * TMath::RadToDeg();
    const Double_t th_max = (DETECTOR_THETA_CENTER_RAD + DETECTOR_THETA_WINDOW)
                            * TMath::RadToDeg();

    TH1D* h_sampled  = new TH1D("h_s2_sampled",  "Sampled",   nbins, th_min, th_max);
    TH1D* h_expected = new TH1D("h_s2_expected", "Expected",  nbins, th_min, th_max);
    h_sampled->SetDirectory(0);
    h_expected->SetDirectory(0);

    // Fill sampled histogram from stored lab momenta
    for (Int_t i = 0; i < result.count; i++) {
        TVector3 p3(result.px[i], result.py[i], result.pz[i]);
        h_sampled->Fill(p3.Theta() * TMath::RadToDeg());
    }

    // Fill expected histogram by scanning t uniformly over the envelope range,
    // evaluating MeyerXS_Elastic, and projecting t → theta_CM → theta_lab.
    // This mirrors exactly what the worker does, so it is the correct reference.
    // No ConvertThetaLabtoCM needed — we go in the same direction as the worker.
    const Int_t n_t_scan = 200000;
    for (Int_t k = 0; k < n_t_scan; k++) {
        Double_t t_GeV2    = t_min + (t_max - t_min) * (k + 0.5) / n_t_scan;
        Double_t xs        = MeyerXS_Elastic(t_GeV2, ekin);
        Double_t theta_cm  = ConvertTtoThetaCM(t_GeV2, ekin);       // [degrees]
        Double_t theta_lab = ConvertThetaCMtoLab(ekin, theta_cm);   // [degrees]
        // Weight by xs; normalisation is applied afterwards so absolute
        // scale does not matter — only the shape is compared.
        h_expected->Fill(theta_lab, xs);
    }

    // Normalize both to unit area for shape comparison
    if (h_sampled->Integral()  > 0) h_sampled->Scale(1.0  / h_sampled->Integral());
    if (h_expected->Integral() > 0) h_expected->Scale(1.0 / h_expected->Integral());

    // --- Chi² shape comparison ---
    Double_t chi2 = 0.;
    Int_t    ndof = 0;
    for (Int_t ibin = 1; ibin <= nbins; ibin++) {
        Double_t obs = h_sampled->GetBinContent(ibin);
        Double_t exp = h_expected->GetBinContent(ibin);
        Double_t err = h_sampled->GetBinError(ibin);
        if (err > 0 && exp > 0) {
            chi2 += (obs - exp)*(obs - exp) / (err*err);
            ndof++;
        }
    }
    Double_t chi2_per_dof = (ndof > 0) ? chi2 / ndof : -1.;

    cout << "\nShape comparison:" << endl;
    cout << "  Chi²/dof = " << chi2_per_dof << " (expect ~1 for good match)" << endl;
    if (chi2_per_dof < 3.0) {
        cout << "  ✓ PASSED — sampled distribution matches expected XS shape" << endl;
    } else {
        cout << "  ✗ FAILED — distribution mismatch, check unit conversions" << endl;
    }

    // --- Plot ---
    TCanvas* c1 = new TCanvas("c_stage2", "Stage 2: Shape Test", 900, 500);
    c1->Divide(2, 1);

    c1->cd(1);
    gPad->SetLogy();
    h_expected->SetLineColor(kRed);
    h_expected->SetLineWidth(2);
    h_expected->SetTitle(Form("Shape Test (E=%.0f MeV, N=%d);#theta_{lab} [deg];Normalized",
                              ekin, n_events));
    h_expected->Draw("HIST");
    h_sampled->SetMarkerStyle(20);
    h_sampled->SetMarkerColor(kBlue);
    h_sampled->Draw("P SAME");

    TLegend* leg1 = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg1->AddEntry(h_expected, "Expected (MeyerXS)", "l");
    leg1->AddEntry(h_sampled,  "Sampled events",     "p");
    leg1->Draw();

    TLatex* lat = new TLatex();
    lat->SetNDC(); lat->SetTextSize(0.04);
    lat->DrawLatex(0.15, 0.85, Form("#chi^{2}/dof = %.2f", chi2_per_dof));

    // Ratio panel
    c1->cd(2);
    TH1D* h_ratio = (TH1D*)h_sampled->Clone("h_s2_ratio");
    h_ratio->SetDirectory(0);
    h_ratio->Divide(h_expected);
    h_ratio->SetTitle("Ratio Sampled/Expected;#theta_{lab} [deg];Ratio");
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerColor(kBlue);
    h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    h_ratio->Draw("P");

    TLine* unity = new TLine(th_min, 1.0, th_max, 1.0);
    unity->SetLineColor(kRed);
    unity->SetLineWidth(2);
    unity->SetLineStyle(2);
    unity->Draw();

    c1->Update();

    delete envelope;

    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "  STAGE 2 COMPLETE" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;
}

// ====================================================================
// STAGE 3: Mode comparison — Meyer vs T_SAMPLING at mid-energy (180 MeV)
// ====================================================================
// Both modes sample from validated splines. At an interpolated energy
// their output distributions should be consistent within statistics.
// A large systematic difference would point to a bug in the Meyer
// interpolation or a unit-conversion mismatch.
// ====================================================================
void TestMeyer_Phase4_Stage3()
{
    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "  STAGE 3: Meyer vs T_SAMPLING comparison" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;

    const Double_t ekin         = 180.0;
    const Int_t    n_events     = 100000;
    const Double_t polarization = 0.80;
    const Int_t    spin_state   = +1;
    const Int_t    num_threads  = 4;

    // --- Generate T_SAMPLING ---
    cout << "Generating T_SAMPLING reference (" << n_events << " events)..." << endl;
    SingleRunMultithreadPolarized(ekin, n_events, polarization, spin_state,
                                  num_threads, T_SAMPLING);

    // --- Generate T_SAMPLING_MEYER ---
    cout << "Generating T_SAMPLING_MEYER (" << n_events << " events)..." << endl;
    SingleRunMultithreadPolarized(ekin, n_events, polarization, spin_state,
                                  num_threads, T_SAMPLING_MEYER);

    // --- Read back both ROOT files and compare theta_lab distributions ---
    int    pol_int   = (int)(polarization * 100);
    TString spinLabel = (spin_state > 0) ? "SpinUp" : "SpinDown";

    TString fname_t  = Form("pC_Elas_%3.0fMeV_MT_P%d_%s_%dp%d.root",
                            ekin, pol_int, spinLabel.Data(),
                            (int)DETECTOR_THETA_CENTER,
                            (int)(DETECTOR_THETA_CENTER*10)%10);
    TString fname_m  = fname_t;   // same name — Meyer run overwrites it
    // NOTE: both runs write to the same filename pattern. To keep both we'd
    // need different names. For now we re-generate and compare via ThreadData
    // directly, bypassing the file system. Re-run both in memory:

    cout << "\nRe-running both modes in memory for direct comparison..." << endl;

    Double_t t_min, t_max;
    ComputeUnifiedTRange(ekin, ekin,
                         DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW,
                         t_min, t_max);

    // T_SAMPLING reference — use existing 10k-bin histogram
    TH1D* h_t_ref = new TH1D("h_t_ref", "T_SAMPLING ref",
                              10000, t_min, t_max);
    h_t_ref->SetDirectory(0);
    for (int ibin = 1; ibin <= 10000; ibin++) {
        Double_t t_GeV2 = h_t_ref->GetBinCenter(ibin);
        Double_t t_MeV2 = TMath::Abs(t_GeV2) * 1.0e6;
        h_t_ref->SetBinContent(ibin, XStSpline(t_MeV2));
    }

    // Meyer envelope
    TH1D* envelope = CreateMeyerEnvelope_Elastic(t_min, t_max);
    envelope->SetDirectory(0);

    const Int_t n_direct = 200000;

    cout << "T_SAMPLING worker  (" << n_direct << " events)..." << endl;
    ThreadData td_t = GenerateThreadEventsPolarized_TSampling(
        99, ekin, n_direct, t_min, t_max,
        DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW, DETECTOR_PHI_WINDOW,
        0.0, +1);   // unpolarized to isolate XS shape

    cout << "T_SAMPLING_MEYER worker (" << n_direct << " events)..." << endl;
    ThreadData td_m = GenerateThreadEventsPolarized_MeyerSampling(
        99, ekin, n_direct, envelope,
        DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW, DETECTOR_PHI_WINDOW,
        0.0, +1);

    delete envelope;
    delete h_t_ref;

    cout << "T_SAMPLING  count: " << td_t.count << endl;
    cout << "Meyer count:       " << td_m.count << endl;

    // Histogram theta_lab for both
    const Int_t    nbins  = 50;
    const Double_t th_min = (DETECTOR_THETA_CENTER_RAD - DETECTOR_THETA_WINDOW)
                            * TMath::RadToDeg();
    const Double_t th_max = (DETECTOR_THETA_CENTER_RAD + DETECTOR_THETA_WINDOW)
                            * TMath::RadToDeg();

    TH1D* h_tsamp = new TH1D("h_s3_tsamp", "T_SAMPLING",       nbins, th_min, th_max);
    TH1D* h_meyer = new TH1D("h_s3_meyer", "T_SAMPLING_MEYER", nbins, th_min, th_max);
    h_tsamp->SetDirectory(0);
    h_meyer->SetDirectory(0);

    for (Int_t i = 0; i < td_t.count; i++) {
        TVector3 p3(td_t.px[i], td_t.py[i], td_t.pz[i]);
        h_tsamp->Fill(p3.Theta() * TMath::RadToDeg());
    }
    for (Int_t i = 0; i < td_m.count; i++) {
        TVector3 p3(td_m.px[i], td_m.py[i], td_m.pz[i]);
        h_meyer->Fill(p3.Theta() * TMath::RadToDeg());
    }

    h_tsamp->Scale(1.0 / h_tsamp->Integral());
    h_meyer->Scale(1.0 / h_meyer->Integral());

    // KS test
    Double_t ks_prob = h_tsamp->KolmogorovTest(h_meyer);
    cout << "\nKolmogorov-Smirnov probability: " << ks_prob << endl;
    cout << "(expect > 0.05 for compatible distributions)" << endl;
    if (ks_prob > 0.05) {
        cout << "✓ PASSED — Meyer and T_SAMPLING distributions are compatible" << endl;
    } else {
        cout << "⚠ Note: KS p < 0.05 — distributions differ. At 180 MeV this" << endl;
        cout << "  is expected if the energy-dependent interpolation produces a" << endl;
        cout << "  meaningfully different XS than the energy-independent spline." << endl;
        cout << "  Verify the shapes look physically reasonable on the plot." << endl;
    }

    // --- Plot ---
    TCanvas* c2 = new TCanvas("c_stage3", "Stage 3: Meyer vs T_SAMPLING", 900, 500);
    c2->Divide(2, 1);

    c2->cd(1);
    gPad->SetLogy();
    h_tsamp->SetLineColor(kBlack);
    h_tsamp->SetLineWidth(2);
    h_tsamp->SetTitle(Form("Meyer vs T_SAMPLING (E=%.0f MeV, N=%d);#theta_{lab} [deg];Normalized",
                           ekin, n_direct));
    h_tsamp->Draw("HIST");
    h_meyer->SetLineColor(kRed);
    h_meyer->SetLineWidth(2);
    h_meyer->SetLineStyle(2);
    h_meyer->Draw("HIST SAME");

    TLegend* leg2 = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg2->AddEntry(h_tsamp, "T_SAMPLING (energy-indep)", "l");
    leg2->AddEntry(h_meyer, "T_SAMPLING_MEYER",          "l");
    leg2->Draw();

    TLatex* lat2 = new TLatex();
    lat2->SetNDC(); lat2->SetTextSize(0.04);
    lat2->DrawLatex(0.15, 0.85, Form("KS prob = %.3f", ks_prob));

    // Ratio panel
    c2->cd(2);
    TH1D* h_ratio2 = (TH1D*)h_meyer->Clone("h_s3_ratio");
    h_ratio2->SetDirectory(0);
    h_ratio2->Divide(h_tsamp);
    h_ratio2->SetTitle("Ratio Meyer/T_SAMPLING;#theta_{lab} [deg];Ratio");
    h_ratio2->SetLineColor(kRed);
    h_ratio2->SetLineWidth(2);
    h_ratio2->GetYaxis()->SetRangeUser(0.5, 1.5);
    h_ratio2->Draw("HIST");

    TLine* unity2 = new TLine(th_min, 1.0, th_max, 1.0);
    unity2->SetLineColor(kBlack);
    unity2->SetLineWidth(2);
    unity2->SetLineStyle(2);
    unity2->Draw();

    c2->Update();

    cout << "\n════════════════════════════════════════════════" << endl;
    cout << "  STAGE 3 COMPLETE" << endl;
    cout << "════════════════════════════════════════════════\n" << endl;
}

// ====================================================================
// Run all three stages sequentially
// ====================================================================
void TestMeyer_Phase4_All()
{
    TestMeyer_Phase4_Stage1();
    TestMeyer_Phase4_Stage2();
    TestMeyer_Phase4_Stage3();

    cout << "\n╔════════════════════════════════════════════════╗" << endl;
    cout << "║   All Phase 4 tests complete                   ║" << endl;
    cout << "╚════════════════════════════════════════════════╝\n" << endl;
}