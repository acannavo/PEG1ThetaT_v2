// DiagnoseXSComparison.C
// ====================================================================
// Diagnostic: compare the two cross section sources used by
// THETA_CM sampling and Meyer sampling at 200 MeV.
//
// This macro exposes a suspected unit mismatch in how GetXSFormula
// is evaluated inside the THETA_CM worker.
//
// Usage:
//   root [0] .L DiagnoseXSComparison.C
//   root [1] DiagnoseXSComparison()
// ====================================================================

void DiagnoseXSComparison()
{
    const Double_t ekin = 200.0;  // MeV

    // ================================================================
    // Kinematic setup — same as both workers
    // ================================================================
    Double_t ekinGeV = ekin / 1000.;
    Double_t mom     = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    Double_t s   = iState.Mag2();
    Double_t pcm = TMath::Sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));

    // theta_CM range for detector at ~16.2° lab
    Double_t theta_cm_min, theta_cm_max;
    ComputeCMAngleRange(ekin, DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW,
                        theta_cm_min, theta_cm_max);
    cout << Form("\nDetector maps to theta_CM: [%.4f, %.4f] degrees", theta_cm_min, theta_cm_max) << endl;

    // Wider plotting range to show full picture
    const Double_t theta_plot_min = 1.0;   // degrees
    const Double_t theta_plot_max = 40.0;  // degrees
    const Int_t    npts           = 500;

    // ================================================================
    // Curve 1: XSLog200MeV evaluated AS THE THETA_CM WORKER DOES IT
    //          i.e. x = theta_CM_in_radians passed to fXS200Op
    //          This is what actually happens at runtime.
    // ================================================================
    TGraph* g_theta_actual = new TGraph(npts);
    g_theta_actual->SetName("g_theta_actual");

    // ================================================================
    // Curve 2: XSLog200MeV evaluated CORRECTLY
    //          i.e. x = theta_CM_in_degrees passed to fXS200Op
    //          This is what the spline intends.
    // ================================================================
    TGraph* g_theta_correct = new TGraph(npts);
    g_theta_correct->SetName("g_theta_correct");

    // ================================================================
    // Curve 3: MeyerXS_Elastic — evaluated correctly via t
    // ================================================================
    TGraph* g_meyer = new TGraph(npts);
    g_meyer->SetName("g_meyer");

    for (Int_t i = 0; i < npts; i++) {
        Double_t theta_cm_deg = theta_plot_min
                              + (theta_plot_max - theta_plot_min) * i / (npts - 1.);
        Double_t theta_cm_rad = theta_cm_deg * TMath::DegToRad();

        // Curve 1: ACTUAL runtime behaviour — x passed as radians to a spline
        //          that expects degrees. This is the bug candidate.
        Double_t xs_actual = fXS200Op(theta_cm_rad);   // x = 0.017 to 0.698 rad
                                                        // but spline expects 1 to 80 deg

        // Curve 2: CORRECT evaluation — pass degrees
        Double_t xs_correct = fXS200Op(theta_cm_deg);  // x = 1 to 40 degrees

        // Curve 3: Meyer — convert theta_CM → t → Meyer spline
        Double_t t_GeV2    = ConvertThetaCMtoT(theta_cm_deg, ekin);
        Double_t xs_meyer  = MeyerXS_Elastic(t_GeV2, ekin);

        g_theta_actual ->SetPoint(i, theta_cm_deg, xs_actual);
        g_theta_correct->SetPoint(i, theta_cm_deg, xs_correct);
        g_meyer        ->SetPoint(i, theta_cm_deg, xs_meyer);
    }

    // ================================================================
    // Print values at the detector angle for direct comparison
    // ================================================================
    Double_t theta_det = 0.5 * (theta_cm_min + theta_cm_max);  // ~17.8 deg
    Double_t xs_act    = fXS200Op(theta_det * TMath::DegToRad());
    Double_t xs_cor    = fXS200Op(theta_det);
    Double_t t_det     = ConvertThetaCMtoT(theta_det, ekin);
    Double_t xs_mey    = MeyerXS_Elastic(t_det, ekin);

    cout << "\n══════════════════════════════════════════════════════════" << endl;
    cout << Form("  Cross section at detector centre (theta_CM = %.2f deg):", theta_det) << endl;
    cout << "══════════════════════════════════════════════════════════" << endl;
    cout << Form("  fXS200Op(%.4f rad) [AS WORKER USES IT]: %.4f mb/sr",
                 theta_det*TMath::DegToRad(), xs_act) << endl;
    cout << Form("  fXS200Op(%.2f deg) [CORRECT EVALUATION]: %.4f mb/sr",
                 theta_det, xs_cor) << endl;
    cout << Form("  MeyerXS_Elastic(t, 200 MeV)            : %.4f mb/sr",
                 xs_mey) << endl;
    cout << Form("\n  Ratio actual/correct : %.4f  (1.0 = no bug)", xs_act/xs_cor) << endl;
    cout << Form("  Ratio Meyer/correct  : %.4f  (1.0 = perfect agreement)", xs_mey/xs_cor) << endl;
    cout << "══════════════════════════════════════════════════════════" << endl;

    // Are they the same value at different points?
    cout << "\n  If the bug exists, fXS200Op(0.31 rad) = fXS200Op(0.31 deg)" << endl;
    cout << "  which means THETA_CM samples at theta_CM = 0.31 deg NOT 17.8 deg" << endl;
    cout << Form("  fXS200Op(0.31) = %.6f  vs  fXS200Op(0.31 rad = 17.75 deg) = %.6f",
                 fXS200Op(0.31), fXS200Op(17.75)) << endl;

    // ================================================================
    // Canvas: 3 panels
    //   Left:  full range 1-40 deg, log scale
    //   Middle: zoom on detector window
    //   Right:  ratio Meyer/correct at each angle
    // ================================================================
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c_diag", "XS Source Comparison", 1800, 600);
    c->Divide(3, 1);

    // --- Panel 1: Full angular range ---
    c->cd(1);
    gPad->SetLogy();
    gPad->SetGrid();
    gPad->SetLeftMargin(0.14);

    g_theta_actual->SetLineColor(kBlue);
    g_theta_actual->SetLineWidth(2);
    g_theta_correct->SetLineColor(kGreen+2);
    g_theta_correct->SetLineWidth(2);
    g_theta_correct->SetLineStyle(2);
    g_meyer->SetLineColor(kOrange+7);
    g_meyer->SetLineWidth(2);
    g_meyer->SetLineStyle(3);

    // Draw frame first
    TH1F* frame1 = gPad->DrawFrame(theta_plot_min, 1e-4, theta_plot_max, 1e4);
    frame1->SetTitle("d#sigma/d#Omega vs #theta_{CM};#theta_{CM} [deg];d#sigma/d#Omega [mb/sr]");
    frame1->GetYaxis()->SetTitleOffset(1.3);

    g_theta_actual ->Draw("L SAME");
    g_theta_correct->Draw("L SAME");
    g_meyer        ->Draw("L SAME");

    // Detector window band
    TBox* box1 = new TBox(theta_cm_min, 1e-4, theta_cm_max, 1e4);
    box1->SetFillColorAlpha(kYellow, 0.3);
    box1->SetLineColor(kYellow+2);
    box1->Draw("SAME");

    TLegend* leg1 = new TLegend(0.35, 0.70, 0.92, 0.88);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.038);
    leg1->AddEntry(g_theta_actual,  "fXS200Op(#theta_{CM} in RAD) #color[600]{[AS WORKER CALLS IT]}", "l");
    leg1->AddEntry(g_theta_correct, "fXS200Op(#theta_{CM} in DEG) [CORRECT]", "l");
    leg1->AddEntry(g_meyer,         "MeyerXS_Elastic(t, 200 MeV)", "l");
    leg1->AddEntry(box1,            "Detector acceptance window", "f");
    leg1->Draw();

    TLatex* lat1 = new TLatex();
    lat1->SetNDC(); lat1->SetTextSize(0.038); lat1->SetTextColor(kBlue);
    lat1->DrawLatex(0.15, 0.92, "E_{kin} = 200 MeV");

    // --- Panel 2: Zoom on detector window ---
    c->cd(2);
    gPad->SetLogy();
    gPad->SetGrid();
    gPad->SetLeftMargin(0.14);

    Double_t zoom_lo = theta_cm_min - 2.0;
    Double_t zoom_hi = theta_cm_max + 2.0;

    TH1F* frame2 = gPad->DrawFrame(zoom_lo, 0.01, zoom_hi, 1000.);
    frame2->SetTitle("Zoom: detector acceptance;#theta_{CM} [deg];d#sigma/d#Omega [mb/sr]");
    frame2->GetYaxis()->SetTitleOffset(1.3);

    g_theta_actual ->Draw("L SAME");
    g_theta_correct->Draw("L SAME");
    g_meyer        ->Draw("L SAME");

    TBox* box2 = new TBox(theta_cm_min, 0.01, theta_cm_max, 1000.);
    box2->SetFillColorAlpha(kYellow, 0.3);
    box2->SetLineColor(kYellow+2);
    box2->Draw("SAME");

    // Annotate the values at detector centre
    TLatex* lat2 = new TLatex();
    lat2->SetNDC(); lat2->SetTextSize(0.038);
    lat2->SetTextColor(kBlue);
    lat2->DrawLatex(0.15, 0.87, Form("fXS(rad): %.2f mb/sr", xs_act));
    lat2->SetTextColor(kGreen+2);
    lat2->DrawLatex(0.15, 0.81, Form("fXS(deg): %.2f mb/sr", xs_cor));
    lat2->SetTextColor(kOrange+7);
    lat2->DrawLatex(0.15, 0.75, Form("Meyer:    %.2f mb/sr", xs_mey));

    // --- Panel 3: Ratio of each curve to fXS200Op(correct) ---
    c->cd(3);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.14);

    TGraph* g_ratio_actual  = new TGraph(npts);
    TGraph* g_ratio_meyer   = new TGraph(npts);

    for (Int_t i = 0; i < npts; i++) {
        Double_t theta_cm_deg = theta_plot_min
                              + (theta_plot_max - theta_plot_min) * i / (npts - 1.);
        Double_t xs_cor_i  = fXS200Op(theta_cm_deg);
        if (xs_cor_i <= 0.) continue;
        Double_t xs_act_i  = fXS200Op(theta_cm_deg * TMath::DegToRad());
        Double_t t_i       = ConvertThetaCMtoT(theta_cm_deg, ekin);
        Double_t xs_mey_i  = MeyerXS_Elastic(t_i, ekin);
        g_ratio_actual->SetPoint(i, theta_cm_deg, xs_act_i / xs_cor_i);
        g_ratio_meyer ->SetPoint(i, theta_cm_deg, xs_mey_i / xs_cor_i);
    }

    TH1F* frame3 = gPad->DrawFrame(theta_plot_min, 0.0, theta_plot_max, 3.0);
    frame3->SetTitle("Ratio to fXS200Op(correct);#theta_{CM} [deg];XS / fXS200Op(correct)");
    frame3->GetYaxis()->SetTitleOffset(1.3);

    g_ratio_actual->SetLineColor(kBlue);
    g_ratio_actual->SetLineWidth(2);
    g_ratio_meyer ->SetLineColor(kOrange+7);
    g_ratio_meyer ->SetLineWidth(2);
    g_ratio_meyer ->SetLineStyle(3);

    g_ratio_actual->Draw("L SAME");
    g_ratio_meyer ->Draw("L SAME");

    TLine* unity = new TLine(theta_plot_min, 1.0, theta_plot_max, 1.0);
    unity->SetLineColor(kGreen+2);
    unity->SetLineWidth(2);
    unity->SetLineStyle(2);
    unity->Draw();

    TBox* box3 = new TBox(theta_cm_min, 0.0, theta_cm_max, 3.0);
    box3->SetFillColorAlpha(kYellow, 0.3);
    box3->SetLineColor(kYellow+2);
    box3->Draw("SAME");

    TLegend* leg3 = new TLegend(0.35, 0.75, 0.92, 0.88);
    leg3->SetBorderSize(0); leg3->SetFillStyle(0); leg3->SetTextSize(0.038);
    leg3->AddEntry(g_ratio_actual, "fXS200Op(rad) / fXS200Op(deg)", "l");
    leg3->AddEntry(g_ratio_meyer,  "MeyerXS / fXS200Op(deg)",       "l");
    leg3->AddEntry(unity,          "Ratio = 1 (perfect agreement)",  "l");
    leg3->Draw();

    c->Update();
    c->SaveAs("DiagnoseXSComparison_200MeV.pdf");
    cout << "\nSaved: DiagnoseXSComparison_200MeV.pdf" << endl;

    // ================================================================
    // Final verdict printed clearly
    // ================================================================
    cout << "\n══════════════════════════════════════════════════════════" << endl;
    cout << "  DIAGNOSIS SUMMARY" << endl;
    cout << "══════════════════════════════════════════════════════════" << endl;
    Bool_t bug_present = (TMath::Abs(xs_act/xs_cor - 1.0) > 0.01);
    if (bug_present) {
        cout << "  UNIT MISMATCH DETECTED in THETA_CM sampling:" << endl;
        cout << "  fXS200Op is called with theta_CM in RADIANS" << endl;
        cout << "  but XSLog200MeV expects theta_CM in DEGREES." << endl;
        cout << Form("  At detector centre: worker uses XS = %.4f mb/sr", xs_act) << endl;
        cout << Form("  Correct value should be:          %.4f mb/sr", xs_cor) << endl;
        cout << Form("  Meyer value:                      %.4f mb/sr", xs_mey) << endl;
        cout << "  The THETA_CM worker evaluates the XS at the WRONG angle." << endl;
    } else {
        cout << "  No unit mismatch detected — fXS200Op is angle-invariant" << endl;
        cout << "  in the narrow detector window. Difference is physical." << endl;
    }
    cout << Form("\n  Meyer vs correct optical model: ratio = %.4f", xs_mey/xs_cor) << endl;
    cout << "══════════════════════════════════════════════════════════\n" << endl;
}