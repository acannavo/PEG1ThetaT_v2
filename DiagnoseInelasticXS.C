// DiagnoseInelasticXS.C
// ====================================================================
// Diagnostic: compare the two inelastic cross section AND analyzing
// power sources used by THETA_CM and Meyer sampling at 200 MeV.
//
//   THETA_CM source: GetInelasticCrossSection / GetInelasticAnalyzingPower
//                    (digitized CSV data, theta_CM in degrees)
//   Meyer source:    MeyerXS_Inelastic / MeyerAP_Inelastic
//                    (spline data converted from t to theta_CM)
//
// Usage:
//   root [0] .L DiagnoseInelasticXS.C
//   root [1] DiagnoseInelasticXS()
// ====================================================================

void DiagnoseInelasticXS()
{
    const Double_t ekin = 200.0;  // MeV

    // Load inelastic data if not already loaded
    if (!g_inelastic_xs || !g_inelastic_AN) LoadInelasticData();

    // ================================================================
    // Kinematic setup
    // ================================================================
    Double_t theta_cm_min, theta_cm_max;
    ComputeCMAngleRange(ekin, DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW,
                        theta_cm_min, theta_cm_max);
    cout << Form("\nDetector maps to theta_CM: [%.4f, %.4f] deg", theta_cm_min, theta_cm_max) << endl;

    const Double_t theta_plot_min = 5.0;
    const Double_t theta_plot_max = 35.0;
    const Int_t    npts           = 500;

    // ================================================================
    // Build all four curves: XS and AP for both sources
    // ================================================================
    TGraph *g_xs_csv   = new TGraph(npts);  // THETA_CM: digitized CSV
    TGraph *g_xs_meyer = new TGraph(npts);  // Meyer spline
    TGraph *g_ap_csv   = new TGraph(npts);  // THETA_CM: digitized CSV
    TGraph *g_ap_meyer = new TGraph(npts);  // Meyer spline

    for (Int_t i = 0; i < npts; i++) {
        Double_t theta_cm_deg = theta_plot_min
                              + (theta_plot_max - theta_plot_min) * i / (npts - 1.);

        // Convert theta_CM → t for Meyer
        Double_t t_GeV2 = ConvertThetaCMtoT(theta_cm_deg, ekin);

        // XS
        Double_t xs_csv   = GetInelasticCrossSection(theta_cm_deg);
        Double_t xs_meyer = MeyerXS_Inelastic(t_GeV2, ekin);

        // AP
        Double_t ap_csv   = GetInelasticAnalyzingPower(theta_cm_deg);
        Double_t ap_meyer = MeyerAP_Inelastic(t_GeV2, ekin);

        g_xs_csv  ->SetPoint(i, theta_cm_deg, xs_csv);
        g_xs_meyer->SetPoint(i, theta_cm_deg, xs_meyer);
        g_ap_csv  ->SetPoint(i, theta_cm_deg, ap_csv);
        g_ap_meyer->SetPoint(i, theta_cm_deg, ap_meyer);
    }

    // ================================================================
    // Print values at detector centre
    // ================================================================
    Double_t theta_det = 0.5 * (theta_cm_min + theta_cm_max);
    Double_t t_det     = ConvertThetaCMtoT(theta_det, ekin);
    Double_t xs_csv_det   = GetInelasticCrossSection(theta_det);
    Double_t xs_meyer_det = MeyerXS_Inelastic(t_det, ekin);
    Double_t ap_csv_det   = GetInelasticAnalyzingPower(theta_det);
    Double_t ap_meyer_det = MeyerAP_Inelastic(t_det, ekin);

    cout << "\n══════════════════════════════════════════════════════════" << endl;
    cout << Form("  INELASTIC at detector centre (theta_CM = %.2f deg):", theta_det) << endl;
    cout << "══════════════════════════════════════════════════════════" << endl;
    cout << "  Cross section:" << endl;
    cout << Form("    CSV  (THETA_CM): %.4f mb/sr", xs_csv_det) << endl;
    cout << Form("    Meyer:           %.4f mb/sr", xs_meyer_det) << endl;
    cout << Form("    Ratio Meyer/CSV: %.4f", xs_meyer_det / xs_csv_det) << endl;
    cout << "  Analyzing power:" << endl;
    cout << Form("    CSV  (THETA_CM): %.4f", ap_csv_det) << endl;
    cout << Form("    Meyer:           %.4f", ap_meyer_det) << endl;
    cout << Form("    Difference:      %.4f", ap_meyer_det - ap_csv_det) << endl;
    cout << "══════════════════════════════════════════════════════════\n" << endl;

    // ================================================================
    // Canvas: 2 rows × 3 columns
    //   Row 1: XS full range, XS zoom, XS ratio
    //   Row 2: AP full range, AP zoom, AP ratio
    // ================================================================
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c_inel_diag", "Inelastic XS & AP Comparison", 1800, 900);
    c->Divide(3, 2);

    auto StyleGraph = [](TGraph* g, Color_t col, Style_t ls, Width_t lw) {
        g->SetLineColor(col); g->SetLineStyle(ls); g->SetLineWidth(lw);
        g->SetMarkerColor(col);
    };

    StyleGraph(g_xs_csv,   kBlue,     1, 2);
    StyleGraph(g_xs_meyer, kOrange+7, 2, 2);
    StyleGraph(g_ap_csv,   kBlue,     1, 2);
    StyleGraph(g_ap_meyer, kOrange+7, 2, 2);

    auto DrawDetWindow = [&](Double_t ylo, Double_t yhi) {
        TBox* b = new TBox(theta_cm_min, ylo, theta_cm_max, yhi);
        b->SetFillColorAlpha(kYellow, 0.30);
        b->SetLineColor(kYellow+2);
        b->Draw("SAME");
    };

    auto MakeLegend = [](Double_t x1, Double_t y1, Double_t x2, Double_t y2,
                         TGraph* g1, const char* l1,
                         TGraph* g2, const char* l2) {
        TLegend* leg = new TLegend(x1, y1, x2, y2);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042);
        leg->AddEntry(g1, l1, "l");
        leg->AddEntry(g2, l2, "l");
        return leg;
    };

    // ── Panel 1: XS full range ──────────────────────────────────────
    c->cd(1);
    gPad->SetLogy(); gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    TH1F* fr1 = gPad->DrawFrame(theta_plot_min, 1e-3, theta_plot_max, 20.);
    fr1->SetTitle("Inelastic d#sigma/d#Omega;#theta_{CM} [deg];d#sigma/d#Omega [mb/sr]");
    fr1->GetYaxis()->SetTitleOffset(1.3);
    DrawDetWindow(1e-3, 20.);
    g_xs_csv  ->Draw("L SAME");
    g_xs_meyer->Draw("L SAME");
    MakeLegend(0.40, 0.75, 0.92, 0.88,
               g_xs_csv,   "CSV data (THETA_{CM} sampling)",
               g_xs_meyer, "Meyer spline")->Draw();
    TLatex* lat = new TLatex();
    lat->SetNDC(); lat->SetTextSize(0.038); lat->SetTextColor(kGray+2);
    lat->DrawLatex(0.15, 0.92, Form("E_{kin} = %.0f MeV", ekin));

    // ── Panel 2: XS zoom on detector window ────────────────────────
    c->cd(2);
    gPad->SetLogy(); gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    Double_t zm = theta_cm_min - 3., zM = theta_cm_max + 3.;
    // Find y-range in zoom window
    Double_t ylo_z = 1e9, yhi_z = 0;
    for (Int_t i = 0; i < npts; i++) {
        Double_t x, y; g_xs_csv->GetPoint(i, x, y);
        if (x >= zm && x <= zM) { ylo_z = TMath::Min(ylo_z, y); yhi_z = TMath::Max(yhi_z, y); }
        g_xs_meyer->GetPoint(i, x, y);
        if (x >= zm && x <= zM) { ylo_z = TMath::Min(ylo_z, y); yhi_z = TMath::Max(yhi_z, y); }
    }
    ylo_z *= 0.5; yhi_z *= 2.;
    TH1F* fr2 = gPad->DrawFrame(zm, ylo_z, zM, yhi_z);
    fr2->SetTitle("Zoom: detector acceptance;#theta_{CM} [deg];d#sigma/d#Omega [mb/sr]");
    fr2->GetYaxis()->SetTitleOffset(1.3);
    DrawDetWindow(ylo_z, yhi_z);
    g_xs_csv  ->Draw("L SAME");
    g_xs_meyer->Draw("L SAME");
    TLatex* lat2 = new TLatex();
    lat2->SetNDC(); lat2->SetTextSize(0.038);
    lat2->SetTextColor(kBlue);
    lat2->DrawLatex(0.15, 0.87, Form("CSV:   %.3f mb/sr", xs_csv_det));
    lat2->SetTextColor(kOrange+7);
    lat2->DrawLatex(0.15, 0.80, Form("Meyer: %.3f mb/sr", xs_meyer_det));
    lat2->SetTextColor(kGray+2);
    lat2->DrawLatex(0.15, 0.73, Form("Ratio: %.3f", xs_meyer_det/xs_csv_det));

    // ── Panel 3: XS ratio Meyer/CSV ────────────────────────────────
    c->cd(3);
    gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    TGraph* g_xs_ratio = new TGraph(npts);
    for (Int_t i = 0; i < npts; i++) {
        Double_t x, y_csv, y_mey;
        g_xs_csv  ->GetPoint(i, x, y_csv);
        g_xs_meyer->GetPoint(i, x, y_mey);
        if (y_csv > 0) g_xs_ratio->SetPoint(i, x, y_mey / y_csv);
        else           g_xs_ratio->SetPoint(i, x, 0);
    }
    g_xs_ratio->SetLineColor(kRed+1); g_xs_ratio->SetLineWidth(2);
    TH1F* fr3 = gPad->DrawFrame(theta_plot_min, 0., theta_plot_max, 3.0);
    fr3->SetTitle("XS Ratio Meyer / CSV;#theta_{CM} [deg];Meyer / CSV");
    fr3->GetYaxis()->SetTitleOffset(1.3);
    DrawDetWindow(0., 3.0);
    g_xs_ratio->Draw("L SAME");
    TLine* u1 = new TLine(theta_plot_min,1,theta_plot_max,1);
    u1->SetLineColor(kGray+2); u1->SetLineWidth(2); u1->SetLineStyle(2); u1->Draw();
    TLatex* lat3 = new TLatex();
    lat3->SetNDC(); lat3->SetTextSize(0.038); lat3->SetTextColor(kRed+1);
    lat3->DrawLatex(0.15, 0.87, Form("At det. centre: %.3f", xs_meyer_det/xs_csv_det));

    // ── Panel 4: AP full range ──────────────────────────────────────
    c->cd(4);
    gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    TH1F* fr4 = gPad->DrawFrame(theta_plot_min, -1.1, theta_plot_max, 1.1);
    fr4->SetTitle("Inelastic A_{N};#theta_{CM} [deg];A_{N}");
    fr4->GetYaxis()->SetTitleOffset(1.3);
    DrawDetWindow(-1.1, 1.1);
    g_ap_csv  ->Draw("L SAME");
    g_ap_meyer->Draw("L SAME");
    TLine* zero = new TLine(theta_plot_min,0,theta_plot_max,0);
    zero->SetLineColor(kGray); zero->SetLineStyle(2); zero->Draw();
    MakeLegend(0.40, 0.15, 0.92, 0.28,
               g_ap_csv,   "CSV data (THETA_{CM} sampling)",
               g_ap_meyer, "Meyer spline")->Draw();

    // ── Panel 5: AP zoom on detector window ────────────────────────
    c->cd(5);
    gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    Double_t ap_lo = -1.1, ap_hi = 1.1;
    TH1F* fr5 = gPad->DrawFrame(zm, ap_lo, zM, ap_hi);
    fr5->SetTitle("Zoom: A_{N} at detector;#theta_{CM} [deg];A_{N}");
    fr5->GetYaxis()->SetTitleOffset(1.3);
    DrawDetWindow(ap_lo, ap_hi);
    g_ap_csv  ->Draw("L SAME");
    g_ap_meyer->Draw("L SAME");
    TLine* zero2 = new TLine(zm,0,zM,0);
    zero2->SetLineColor(kGray); zero2->SetLineStyle(2); zero2->Draw();
    TLatex* lat5 = new TLatex();
    lat5->SetNDC(); lat5->SetTextSize(0.038);
    lat5->SetTextColor(kBlue);
    lat5->DrawLatex(0.15, 0.87, Form("CSV:   A_N = %.4f", ap_csv_det));
    lat5->SetTextColor(kOrange+7);
    lat5->DrawLatex(0.15, 0.80, Form("Meyer: A_N = %.4f", ap_meyer_det));
    lat5->SetTextColor(kGray+2);
    lat5->DrawLatex(0.15, 0.73, Form("#DeltaA_N = %.4f", ap_meyer_det - ap_csv_det));

    // ── Panel 6: AP difference Meyer - CSV ─────────────────────────
    c->cd(6);
    gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    TGraph* g_ap_diff = new TGraph(npts);
    for (Int_t i = 0; i < npts; i++) {
        Double_t x, y_csv, y_mey;
        g_ap_csv  ->GetPoint(i, x, y_csv);
        g_ap_meyer->GetPoint(i, x, y_mey);
        g_ap_diff ->SetPoint(i, x, y_mey - y_csv);
    }
    g_ap_diff->SetLineColor(kRed+1); g_ap_diff->SetLineWidth(2);
    // Find range
    Double_t diff_max = 0;
    for (Int_t i = 0; i < npts; i++) {
        Double_t x, y; g_ap_diff->GetPoint(i, x, y);
        diff_max = TMath::Max(diff_max, TMath::Abs(y));
    }
    diff_max *= 1.2;
    TH1F* fr6 = gPad->DrawFrame(theta_plot_min, -diff_max, theta_plot_max, diff_max);
    fr6->SetTitle("A_{N} Difference (Meyer #minus CSV);#theta_{CM} [deg];#Delta A_{N}");
    fr6->GetYaxis()->SetTitleOffset(1.3);
    DrawDetWindow(-diff_max, diff_max);
    g_ap_diff->Draw("L SAME");
    TLine* u6 = new TLine(theta_plot_min,0,theta_plot_max,0);
    u6->SetLineColor(kGray+2); u6->SetLineWidth(2); u6->SetLineStyle(2); u6->Draw();
    TLatex* lat6 = new TLatex();
    lat6->SetNDC(); lat6->SetTextSize(0.038); lat6->SetTextColor(kRed+1);
    lat6->DrawLatex(0.15, 0.87, Form("At det. centre: #Delta A_N = %.4f", ap_meyer_det - ap_csv_det));

    c->Update();
    c->SaveAs("DiagnoseInelasticXS_200MeV.pdf");
    cout << "Saved: DiagnoseInelasticXS_200MeV.pdf" << endl;

    // ================================================================
    // Verdict
    // ================================================================
    cout << "\n══════════════════════════════════════════════════════════" << endl;
    cout << "  DIAGNOSIS SUMMARY (INELASTIC)" << endl;
    cout << "══════════════════════════════════════════════════════════" << endl;
    Bool_t xs_mismatch = TMath::Abs(xs_meyer_det/xs_csv_det - 1.0) > 0.05;
    Bool_t ap_mismatch = TMath::Abs(ap_meyer_det - ap_csv_det) > 0.05;
    cout << Form("  XS at det. centre: CSV=%.4f, Meyer=%.4f, ratio=%.4f  %s",
                 xs_csv_det, xs_meyer_det, xs_meyer_det/xs_csv_det,
                 xs_mismatch ? "<-- SIGNIFICANT DIFFERENCE" : "OK") << endl;
    cout << Form("  AP at det. centre: CSV=%.4f, Meyer=%.4f, diff=%.4f  %s",
                 ap_csv_det, ap_meyer_det, ap_meyer_det-ap_csv_det,
                 ap_mismatch ? "<-- SIGNIFICANT DIFFERENCE" : "OK") << endl;
    if (xs_mismatch || ap_mismatch) {
        cout << "\n  The two inelastic sources disagree — this explains the" << endl;
        cout << "  theta_lab and asymmetry differences between the two modes." << endl;
        cout << "  Both implementations are technically correct; the" << endl;
        cout << "  discrepancy reflects different underlying data sources." << endl;
    } else {
        cout << "\n  Both sources agree within 5% — discrepancy is elsewhere." << endl;
    }
    cout << "══════════════════════════════════════════════════════════\n" << endl;
}