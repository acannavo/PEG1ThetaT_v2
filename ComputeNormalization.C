// ComputeNormalization.C
// ====================================================================
// Compute the inelastic/elastic cross section ratio integrated over
// the detector acceptance window, using Meyer energy-dependent splines.
//
// This ratio is the weight you need to apply to the inelastic event
// sample when combining elastic and inelastic samples that were each
// generated with the same number of events N:
//
//   N_inel_physical = N_generated * R_inel_elas
//
// where R_inel_elas = sigma_inel_det / sigma_elas_det
//
// Physics:
//   sigma_det = 2*pi * integral[ dσ/dΩ(θ_CM) * sin(θ_CM) ] dθ_CM
//
//   The integral runs over the θ_CM range that maps into the detector
//   acceptance window (computed via ComputeCMAngleRange).
//   The factor 2*pi comes from integrating over all phi (isotropic
//   in phi for the purpose of total rate — polarisation affects the
//   phi distribution but not the total rate).
//
// Usage:
//   root [0] .L ComputeNormalization.C
//   root [1] ComputeNormalization()           // uses DetectorConfig defaults
//   root [1] ComputeNormalization(200., 0.80) // explicit energy and polarization
//
// Prerequisites: all PEG modules must be loaded (rootlogon.C)
// ====================================================================

#include <sstream>

// ====================================================================
// Helper: suppress the ComputeCMAngleRange printout
// ====================================================================
void GetCMRange(Double_t ekin, Double_t& theta_cm_min, Double_t& theta_cm_max)
{
    std::streambuf* old = std::cout.rdbuf();
    std::ostringstream devnull;
    std::cout.rdbuf(devnull.rdbuf());
    ComputeCMAngleRange(ekin,
                        DETECTOR_THETA_CENTER_RAD,
                        DETECTOR_THETA_WINDOW,
                        theta_cm_min, theta_cm_max);
    std::cout.rdbuf(old);
}

// ====================================================================
// IntegrateXS: integrate dσ/dΩ over the detector θ_CM window
//
//   sigma_det = 2*pi * ∫ dσ/dΩ(θ_CM, ekin) * sin(θ_CM) dθ_CM
//
// Uses trapezoidal integration over n_steps uniform θ_CM points.
// Input:  theta_cm_min/max [degrees], ekin [MeV]
// Output: integrated cross section [mb]  (dσ/dΩ in mb/sr × sr)
// ====================================================================
Double_t IntegrateXS(Double_t theta_cm_min_deg,
                     Double_t theta_cm_max_deg,
                     Double_t ekin,
                     Bool_t   inelastic,
                     Int_t    n_steps = 10000)
{
    Double_t dtheta = (theta_cm_max_deg - theta_cm_min_deg) / n_steps;
    Double_t sum    = 0.;

    for (Int_t i = 0; i <= n_steps; i++) {
        Double_t theta_cm_deg = theta_cm_min_deg + i * dtheta;
        Double_t theta_cm_rad = theta_cm_deg * TMath::DegToRad();

        // Convert θ_CM → t for Meyer evaluation
        Double_t t_GeV2 = ConvertThetaCMtoT(theta_cm_deg, ekin);

        Double_t xs;  // dσ/dΩ [mb/sr]
        if (inelastic)
            xs = MeyerXS_Inelastic(t_GeV2, ekin);
        else
            xs = MeyerXS_Elastic(t_GeV2, ekin);

        // Trapezoidal weight (half weight at endpoints)
        Double_t w = (i == 0 || i == n_steps) ? 0.5 : 1.0;

        // dΩ element: sin(θ_CM) dθ_CM (the 2π from φ added outside)
        sum += w * xs * TMath::Sin(theta_cm_rad) * dtheta * TMath::DegToRad();
    }

    // Multiply by 2π for the azimuthal integration
    return 2. * TMath::Pi() * sum;  // [mb]
}

// ====================================================================
// ComputeNormalization
// ====================================================================
void ComputeNormalization(Double_t ekin         = 200.0,  // MeV
                          Double_t polarization = 0.80,
                          Int_t    n_steps      = 10000)
{
    cout << "\n╔══════════════════════════════════════════════════════════╗" << endl;
    cout << "║     Inelastic / Elastic Normalization Factor            ║" << endl;
    cout << "║     Meyer energy-dependent sampling                     ║" << endl;
    cout << "╚══════════════════════════════════════════════════════════╝\n" << endl;

    // ----------------------------------------------------------------
    // Print current detector configuration
    // ----------------------------------------------------------------
    cout << "Detector configuration:" << endl;
    cout << Form("  θ_lab centre : %.2f deg", DETECTOR_THETA_CENTER) << endl;
    cout << Form("  θ_lab window : ±%.4f rad = ±%.4f deg",
                 DETECTOR_THETA_WINDOW,
                 DETECTOR_THETA_WINDOW * TMath::RadToDeg()) << endl;
    cout << Form("  φ window     : ±%.4f rad = ±%.4f deg",
                 DETECTOR_PHI_WINDOW,
                 DETECTOR_PHI_WINDOW * TMath::RadToDeg()) << endl;
    cout << Form("  Beam energy  : %.1f MeV", ekin) << endl;
    cout << Form("  Polarization : %.0f%%\n", polarization*100) << endl;

    // ----------------------------------------------------------------
    // Get θ_CM range for the detector acceptance
    // ----------------------------------------------------------------
    Double_t theta_cm_min, theta_cm_max;
    GetCMRange(ekin, theta_cm_min, theta_cm_max);

    cout << Form("θ_CM acceptance range: [%.4f, %.4f] deg", theta_cm_min, theta_cm_max) << endl;
    cout << Form("  (maps to θ_lab: [%.4f, %.4f] deg)",
                 DETECTOR_THETA_CENTER - DETECTOR_THETA_WINDOW * TMath::RadToDeg(),
                 DETECTOR_THETA_CENTER + DETECTOR_THETA_WINDOW * TMath::RadToDeg()) << endl;
    cout << endl;

    // ----------------------------------------------------------------
    // Integrate elastic and inelastic XS over detector acceptance
    // ----------------------------------------------------------------
    Double_t sigma_elas = IntegrateXS(theta_cm_min, theta_cm_max, ekin, false, n_steps);
    Double_t sigma_inel = IntegrateXS(theta_cm_min, theta_cm_max, ekin, true,  n_steps);

    Double_t ratio = sigma_inel / sigma_elas;

    // ----------------------------------------------------------------
    // Print results
    // ----------------------------------------------------------------
    cout << "══════════════════════════════════════════════════════════" << endl;
    cout << "  Integrated cross sections over detector acceptance" << endl;
    cout << "══════════════════════════════════════════════════════════" << endl;
    cout << Form("  σ_elastic   (det) = %.6f mb", sigma_elas) << endl;
    cout << Form("  σ_inelastic (det) = %.6f mb", sigma_inel) << endl;
    cout << Form("  Ratio  R = σ_inel / σ_elas = %.6f", ratio) << endl;
    cout << "══════════════════════════════════════════════════════════" << endl;

    cout << "\n  Interpretation:" << endl;
    cout << Form("  If you generated N events for elastic AND N events for inelastic,") << endl;
    cout << Form("  weight each inelastic event by R = %.6f", ratio) << endl;
    cout << Form("  i.e. inelastic counts as %.4f%% of elastic rate.", ratio*100.) << endl;

    // ----------------------------------------------------------------
    // Cross-check: values at detector centre angle
    // ----------------------------------------------------------------
    Double_t theta_cm_ctr = 0.5 * (theta_cm_min + theta_cm_max);
    Double_t t_ctr        = ConvertThetaCMtoT(theta_cm_ctr, ekin);
    Double_t xs_el_ctr    = MeyerXS_Elastic  (t_ctr, ekin);
    Double_t xs_in_ctr    = MeyerXS_Inelastic(t_ctr, ekin);

    cout << "\n  Cross-check — dσ/dΩ at detector centre (θ_CM = " 
         << theta_cm_ctr << " deg):" << endl;
    cout << Form("    Elastic  : %.4f mb/sr", xs_el_ctr) << endl;
    cout << Form("    Inelastic: %.4f mb/sr", xs_in_ctr) << endl;
    cout << Form("    Point ratio: %.6f  (vs integrated %.6f)",
                 xs_in_ctr / xs_el_ctr, ratio) << endl;
    cout << "══════════════════════════════════════════════════════════\n" << endl;

    // ----------------------------------------------------------------
    // Scan over the Meyer energy range for reference table
    // ----------------------------------------------------------------
    cout << "  Reference table across 160-200 MeV (same detector config):" << endl;
    cout << Form("  %-8s  %-14s  %-14s  %-10s", "E [MeV]", "σ_elas [mb]", "σ_inel [mb]", "R = inel/elas") << endl;
    cout << "  ──────────────────────────────────────────────────────" << endl;

    Double_t energies[5] = {160., 170., 180., 190., 200.};
    for (Int_t ie = 0; ie < 5; ie++) {
        Double_t e = energies[ie];
        Double_t lo, hi;
        GetCMRange(e, lo, hi);
        Double_t se = IntegrateXS(lo, hi, e, false, n_steps);
        Double_t si = IntegrateXS(lo, hi, e, true,  n_steps);
        TString marker = (TMath::Abs(e - ekin) < 0.5) ? " ← current" : "";
        cout << Form("  %-8.0f  %-14.6f  %-14.6f  %-10.6f%s",
                     e, se, si, si/se, marker.Data()) << endl;
    }
    cout << endl;

    // ----------------------------------------------------------------
    // Diagnostic plots — 3-panel canvas
    //
    // Panel 1: dσ/dΩ vs θ_CM over full 5-35° range (log scale)
    //          Vertical band marks the detector acceptance window.
    //          Shows the full shape of both cross sections in context.
    //
    // Panel 2: Point-by-point ratio dσ/dΩ_inel / dσ/dΩ_elas vs θ_CM
    //          Dashed vertical lines mark the detector window.
    //          Horizontal dashed line marks the integrated ratio R.
    //
    // Panel 3: Zoom into the detector window (linear scale)
    //          Shows the actual integrand: dσ/dΩ × sin(θ_CM)
    //          The shaded area IS the integral — visually confirms
    //          what is being summed.
    // ----------------------------------------------------------------
    const Double_t plot_lo = 5.0;   // full-range plot limits [deg]
    const Double_t plot_hi = 35.0;
    const Int_t    npts_full = 500;
    const Int_t    npts_zoom = 200;

    // Full-range curves
    TGraph* g_el_full = new TGraph(npts_full);
    TGraph* g_in_full = new TGraph(npts_full);
    TGraph* g_ratio   = new TGraph(npts_full);  // point-by-point ratio

    for (Int_t i = 0; i < npts_full; i++) {
        Double_t th = plot_lo + (plot_hi - plot_lo) * i / (npts_full - 1.);
        Double_t t  = ConvertThetaCMtoT(th, ekin);
        Double_t xe = MeyerXS_Elastic  (t, ekin);
        Double_t xi = MeyerXS_Inelastic(t, ekin);
        g_el_full->SetPoint(i, th, xe);
        g_in_full->SetPoint(i, th, xi);
        g_ratio  ->SetPoint(i, th, (xe > 0) ? xi/xe : 0.);
    }

    // Zoom curves — integrand: dσ/dΩ * sin(θ_CM) [mb/sr]
    TGraph* g_el_zoom = new TGraph(npts_zoom);
    TGraph* g_in_zoom = new TGraph(npts_zoom);

    for (Int_t i = 0; i < npts_zoom; i++) {
        Double_t th  = theta_cm_min + (theta_cm_max - theta_cm_min) * i / (npts_zoom-1.);
        Double_t t   = ConvertThetaCMtoT(th, ekin);
        Double_t sin_th = TMath::Sin(th * TMath::DegToRad());
        g_el_zoom->SetPoint(i, th, MeyerXS_Elastic  (t, ekin) * sin_th);
        g_in_zoom->SetPoint(i, th, MeyerXS_Inelastic(t, ekin) * sin_th);
    }

    // Style
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.055);

    TCanvas* c = new TCanvas("c_norm", "Cross Section Normalization", 1500, 550);
    c->SetFillColor(kWhite);
    c->Divide(3, 1);

    // ── Panel 1: Full range, log scale ──────────────────────────────
    c->cd(1);
    gPad->SetLogy();
    gPad->SetGrid();
    gPad->SetLeftMargin(0.15);

    g_el_full->SetLineColor(kAzure+1);   g_el_full->SetLineWidth(2);
    g_in_full->SetLineColor(kOrange+7);  g_in_full->SetLineWidth(2);
    g_in_full->SetLineStyle(2);

    // Find y range
    Double_t ymin_f = 1e9, ymax_f = 0.;
    for (Int_t i = 0; i < npts_full; i++) {
        Double_t x, y;
        g_el_full->GetPoint(i, x, y); if (y>0) { ymin_f=TMath::Min(ymin_f,y); ymax_f=TMath::Max(ymax_f,y); }
        g_in_full->GetPoint(i, x, y); if (y>0) { ymin_f=TMath::Min(ymin_f,y); ymax_f=TMath::Max(ymax_f,y); }
    }

    TH1F* fr1 = gPad->DrawFrame(plot_lo, ymin_f*0.5, plot_hi, ymax_f*2.);
    fr1->SetTitle(Form("d#sigma/d#Omega  (E=%.0f MeV)"
                       ";#theta_{CM} [deg];d#sigma/d#Omega [mb/sr]", ekin));
    fr1->GetYaxis()->SetTitleOffset(1.5);

    // Detector window band
    TBox* band1 = new TBox(theta_cm_min, ymin_f*0.5, theta_cm_max, ymax_f*2.);
    band1->SetFillColorAlpha(kYellow, 0.35);
    band1->SetLineColor(kYellow+2);
    band1->Draw("SAME");

    g_el_full->Draw("L SAME");
    g_in_full->Draw("L SAME");

    TLegend* leg1 = new TLegend(0.35, 0.72, 0.94, 0.88);
    leg1->SetBorderSize(0); leg1->SetFillStyle(0); leg1->SetTextSize(0.042);
    leg1->AddEntry(g_el_full, "Elastic (Meyer)",   "l");
    leg1->AddEntry(g_in_full, "Inelastic (Meyer)", "l");
    leg1->AddEntry(band1,     "Detector window",   "f");
    leg1->Draw();

    // ── Panel 2: Point-by-point ratio vs θ_CM ───────────────────────
    c->cd(2);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.15);

    g_ratio->SetLineColor(kRed+1);
    g_ratio->SetLineWidth(2);

    // Find ratio range in detector window ± margin
    Double_t rlo=1e9, rhi=0.;
    for (Int_t i = 0; i < npts_full; i++) {
        Double_t x, y; g_ratio->GetPoint(i, x, y);
        if (x >= plot_lo && x <= plot_hi && y > 0) {
            rlo = TMath::Min(rlo, y);
            rhi = TMath::Max(rhi, y);
        }
    }
    Double_t rpad = 0.2 * (rhi - rlo);

    TH1F* fr2 = gPad->DrawFrame(plot_lo, TMath::Max(0., rlo-rpad), plot_hi, rhi+rpad);
    fr2->SetTitle(Form("Ratio d#sigma_{inel}/d#sigma_{elas}"
                       ";#theta_{CM} [deg];#sigma_{inel}/#sigma_{elas}"));
    fr2->GetYaxis()->SetTitleOffset(1.5);

    // Detector window lines
    TLine* lwin_lo = new TLine(theta_cm_min, TMath::Max(0.,rlo-rpad), theta_cm_min, rhi+rpad);
    TLine* lwin_hi = new TLine(theta_cm_max, TMath::Max(0.,rlo-rpad), theta_cm_max, rhi+rpad);
    lwin_lo->SetLineColor(kYellow+2); lwin_lo->SetLineWidth(2); lwin_lo->SetLineStyle(2);
    lwin_hi->SetLineColor(kYellow+2); lwin_hi->SetLineWidth(2); lwin_hi->SetLineStyle(2);

    // Horizontal line at the integrated ratio R
    TLine* lratio = new TLine(plot_lo, ratio, plot_hi, ratio);
    lratio->SetLineColor(kGreen+2); lratio->SetLineWidth(2); lratio->SetLineStyle(2);

    lwin_lo->Draw(); lwin_hi->Draw(); lratio->Draw();
    g_ratio->Draw("L SAME");

    TLatex* lat2 = new TLatex();
    lat2->SetNDC(); lat2->SetTextSize(0.042);
    lat2->SetTextColor(kGreen+2);
    lat2->DrawLatex(0.18, 0.88, Form("R_{integrated} = %.5f", ratio));
    lat2->SetTextColor(kGray+2);
    lat2->DrawLatex(0.18, 0.80, "Dashed: detector window");

    // ── Panel 3: Zoom — integrand dσ/dΩ × sin(θ) ───────────────────
    c->cd(3);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.15);

    g_el_zoom->SetLineColor(kAzure+1);  g_el_zoom->SetLineWidth(2);
    g_in_zoom->SetLineColor(kOrange+7); g_in_zoom->SetLineWidth(2);
    g_in_zoom->SetLineStyle(2);

    // Scan point arrays for min/max — TGraph::GetMinimum() returns -1111 if unset
    Double_t zmin_el=1e9, zmax_el=0., zmin_in=1e9, zmax_in=0.;
    for (Int_t i = 0; i < npts_zoom; i++) {
        Double_t x, y;
        g_el_zoom->GetPoint(i, x, y);
        zmin_el = TMath::Min(zmin_el, y); zmax_el = TMath::Max(zmax_el, y);
        g_in_zoom->GetPoint(i, x, y);
        zmin_in = TMath::Min(zmin_in, y); zmax_in = TMath::Max(zmax_in, y);
    }
    Double_t zmin = 0.;                               // anchor fills at zero
    Double_t zmax = 1.3 * TMath::Max(zmax_el, zmax_in);

    TH1F* fr3 = gPad->DrawFrame(theta_cm_min, zmin, theta_cm_max, zmax);
    fr3->SetTitle("Integrand: d#sigma/d#Omega #times sin(#theta_{CM})"
                  ";#theta_{CM} [deg];d#sigma/d#Omega #times sin#theta [mb/sr]");
    fr3->GetYaxis()->SetTitleOffset(1.7);

    // Fill area under each integrand curve.
    // Close polygon at y=0 (frame baseline) so the shaded area is visible.
    TGraph* g_el_fill = new TGraph(npts_zoom + 2);
    for (Int_t i = 0; i < npts_zoom; i++) {
        Double_t x, y; g_el_zoom->GetPoint(i, x, y);
        g_el_fill->SetPoint(i, x, y);
    }
    g_el_fill->SetPoint(npts_zoom,     theta_cm_max, 0.);
    g_el_fill->SetPoint(npts_zoom + 1, theta_cm_min, 0.);
    g_el_fill->SetFillColorAlpha(kAzure+1, 0.25);
    g_el_fill->SetLineWidth(0);
    g_el_fill->Draw("F SAME");

    TGraph* g_in_fill = new TGraph(npts_zoom + 2);
    for (Int_t i = 0; i < npts_zoom; i++) {
        Double_t x, y; g_in_zoom->GetPoint(i, x, y);
        g_in_fill->SetPoint(i, x, y);
    }
    g_in_fill->SetPoint(npts_zoom,     theta_cm_max, 0.);
    g_in_fill->SetPoint(npts_zoom + 1, theta_cm_min, 0.);
    g_in_fill->SetFillColorAlpha(kOrange+7, 0.25);
    g_in_fill->SetLineWidth(0);
    g_in_fill->Draw("F SAME");

    // Draw lines on top of fills
    g_el_zoom->Draw("L SAME");
    g_in_zoom->Draw("L SAME");

    TLegend* leg3 = new TLegend(0.16, 0.72, 0.92, 0.88);
    leg3->SetBorderSize(0); leg3->SetFillStyle(0); leg3->SetTextSize(0.040);
    leg3->AddEntry(g_el_zoom, Form("Elastic   (#int = %.5f mb/2#pi)", sigma_elas / (2.*TMath::Pi())), "lf");
    leg3->AddEntry(g_in_zoom, Form("Inelastic (#int = %.5f mb/2#pi)", sigma_inel / (2.*TMath::Pi())), "lf");
    leg3->Draw();

    TLatex* lat3 = new TLatex();
    lat3->SetNDC(); lat3->SetTextSize(0.038); lat3->SetTextColor(kGray+2);
    lat3->DrawLatex(0.16, 0.18, "Shaded area #propto #sigma_{det}");

    c->Update();
    TString outname = Form("Normalization_%.0fMeV_%dp%d.pdf",
                           ekin,
                           (int)DETECTOR_THETA_CENTER,
                           (int)(DETECTOR_THETA_CENTER*10)%10);
    c->SaveAs(outname.Data());
    cout << "Plot saved: " << outname << endl;
}