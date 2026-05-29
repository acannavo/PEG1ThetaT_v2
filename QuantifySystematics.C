// QuantifySystematics.C
// ====================================================================
// Quantify the systematic difference between THETA_CM and Meyer
// sampling modes for both elastic and inelastic scattering.
//
// For each energy in {160, 170, 180, 190, 200} MeV, this macro:
//   1. Compares dσ/dt from both XS sources at the detector angle
//   2. Compares A_N from both sources
//   3. Computes the expected asymmetry A = P * A_N for each mode
//   4. Estimates the pull |A_theta - A_meyer| / sigma_stat
//      as a function of event count N
//
// Output:
//   - Console table: all numbers at a glance
//   - PDF: XS ratio, A_N comparison, and asymmetry pull vs N plots
//
// Usage:
//   root [0] .L QuantifySystematics.C
//   root [1] QuantifySystematics()
// ====================================================================

void QuantifySystematics()
{
    // Load inelastic data if needed
    if (!g_inelastic_xs || !g_inelastic_AN) LoadInelasticData();

    const Double_t polarization = 0.80;

    // Energies to scan
    const Int_t    nE    = 5;
    Double_t ekin_arr[5] = {160., 170., 180., 190., 200.};

    // ================================================================
    // Storage for results
    // ================================================================
    // Elastic
    Double_t xs_theta_el[5], xs_meyer_el[5], xs_ratio_el[5];
    Double_t an_theta_el[5], an_meyer_el[5], an_diff_el[5];
    Double_t asym_theta_el[5], asym_meyer_el[5];

    // Inelastic
    Double_t xs_theta_in[5], xs_meyer_in[5], xs_ratio_in[5];
    Double_t an_theta_in[5], an_meyer_in[5], an_diff_in[5];
    Double_t asym_theta_in[5], asym_meyer_in[5];

    // ================================================================
    // Compute for each energy
    // ================================================================
    for (Int_t ie = 0; ie < nE; ie++) {
        Double_t ekin = ekin_arr[ie];

        // Detector CM angle at this energy
        Double_t theta_cm_lo, theta_cm_hi;
        {
            // Suppress ComputeCMAngleRange printout
            std::streambuf* old = std::cout.rdbuf();
            std::ostringstream devnull;
            std::cout.rdbuf(devnull.rdbuf());
            ComputeCMAngleRange(ekin, DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW,
                                theta_cm_lo, theta_cm_hi);
            std::cout.rdbuf(old);
        }
        Double_t theta_det = 0.5 * (theta_cm_lo + theta_cm_hi);
        Double_t t_det     = ConvertThetaCMtoT(theta_det, ekin);

        // ── ELASTIC ──────────────────────────────────────────────
        // THETA_CM source: optical model formula (evaluated correctly in degrees)
        const char* formula = GetXSFormula(ekin);
        TF1 f_el(Form("f_el_%d", ie), formula, 1., 80.);
        xs_theta_el[ie] = f_el.Eval(theta_det);

        // Meyer source
        xs_meyer_el[ie] = MeyerXS_Elastic(t_det, ekin);
        xs_ratio_el[ie] = xs_meyer_el[ie] / xs_theta_el[ie];

        // A_N: THETA_CM uses GetElasticAnalyzingPower(ekin, theta_lab_deg)
        Double_t theta_lab_det = ConvertThetaCMtoLab(ekin, theta_det);
        an_theta_el[ie] = GetElasticAnalyzingPower(ekin, theta_lab_det);
        an_meyer_el[ie] = MeyerAP_Elastic(t_det, ekin);
        an_diff_el[ie]  = an_meyer_el[ie] - an_theta_el[ie];

        asym_theta_el[ie] = polarization * an_theta_el[ie];
        asym_meyer_el[ie] = polarization * an_meyer_el[ie];

        // ── INELASTIC ─────────────────────────────────────────────
        // THETA_CM source: digitized CSV (only at 200 MeV — for other
        // energies we still use it as an approximation since it's the
        // only non-Meyer source available)
        xs_theta_in[ie] = GetInelasticCrossSection(theta_det);
        xs_meyer_in[ie] = MeyerXS_Inelastic(t_det, ekin);
        xs_ratio_in[ie] = xs_meyer_in[ie] / xs_theta_in[ie];

        an_theta_in[ie] = GetInelasticAnalyzingPower(theta_det);
        an_meyer_in[ie] = MeyerAP_Inelastic(t_det, ekin);
        an_diff_in[ie]  = an_meyer_in[ie] - an_theta_in[ie];

        asym_theta_in[ie] = polarization * an_theta_in[ie];
        asym_meyer_in[ie] = polarization * an_meyer_in[ie];
    }

    // ================================================================
    // Console table
    // ================================================================
    cout << "\n";
    cout << "╔══════════════════════════════════════════════════════════════════════════════╗" << endl;
    cout << "║            SYSTEMATIC COMPARISON: THETA_CM vs MEYER SAMPLING               ║" << endl;
    cout << "╠══════════════════════════════════════════════════════════════════════════════╣" << endl;
    cout << "║  ELASTIC                                                                    ║" << endl;
    cout << "╠══════╦══════════╦══════════╦════════╦══════════╦══════════╦════════════════╣" << endl;
    cout << "║ E    ║ XS_THETA ║ XS_MEYER ║ Ratio  ║ AN_THETA ║ AN_MEYER ║    ΔA_N        ║" << endl;
    cout << "║ MeV  ║ mb/sr    ║ mb/sr    ║ M/T    ║          ║          ║ Meyer-Theta    ║" << endl;
    cout << "╠══════╬══════════╬══════════╬════════╬══════════╬══════════╬════════════════╣" << endl;
    for (Int_t ie = 0; ie < nE; ie++) {
        cout << Form("║ %4.0f ║ %8.4f ║ %8.4f ║ %6.4f ║ %8.4f ║ %8.4f ║ %+.4f        ║",
                     ekin_arr[ie],
                     xs_theta_el[ie], xs_meyer_el[ie], xs_ratio_el[ie],
                     an_theta_el[ie], an_meyer_el[ie], an_diff_el[ie]) << endl;
    }
    cout << "╠══════════════════════════════════════════════════════════════════════════════╣" << endl;
    cout << "║  INELASTIC  (note: CSV source is 200 MeV only — used at all energies)      ║" << endl;
    cout << "╠══════╦══════════╦══════════╦════════╦══════════╦══════════╦════════════════╣" << endl;
    cout << "║ E    ║ XS_CSV   ║ XS_MEYER ║ Ratio  ║ AN_CSV   ║ AN_MEYER ║    ΔA_N        ║" << endl;
    cout << "╠══════╬══════════╬══════════╬════════╬══════════╬══════════╬════════════════╣" << endl;
    for (Int_t ie = 0; ie < nE; ie++) {
        cout << Form("║ %4.0f ║ %8.4f ║ %8.4f ║ %6.4f ║ %8.4f ║ %8.4f ║ %+.4f        ║",
                     ekin_arr[ie],
                     xs_theta_in[ie], xs_meyer_in[ie], xs_ratio_in[ie],
                     an_theta_in[ie], an_meyer_in[ie], an_diff_in[ie]) << endl;
    }
    cout << "╚══════╩══════════╩══════════╩════════╩══════════╩══════════╩════════════════╝" << endl;

    // ================================================================
    // Pull as function of N for elastic at each energy
    // ================================================================
    // |A_theta - A_meyer| / sqrt(2) * sigma_stat
    // where sigma_stat = sqrt(4*R*L/(R+L)^3) ≈ 1/sqrt(N) for large N
    // A ≈ P*A_N, so delta_A = P * |AN_theta - AN_meyer|
    // Pull = delta_A / (1/sqrt(N)) = delta_A * sqrt(N)
    // N_half = (delta_A)^-2  (pull = 1 at this N)

    cout << "\n════════════════════════════════════════════════════════════════" << endl;
    cout << "  Events needed to see the systematic at 1σ / 2σ / 3σ" << endl;
    cout << "════════════════════════════════════════════════════════════════" << endl;
    cout << Form("  %-10s  %-12s  %-12s  %-12s", "Energy", "Elastic", "Inelastic", "Note") << endl;
    cout << "  (MeV)       (1σ/2σ/3σ)   (1σ/2σ/3σ)" << endl;
    cout << "────────────────────────────────────────────────────────────────" << endl;

    for (Int_t ie = 0; ie < nE; ie++) {
        // For asymmetry: sigma_stat(A) ≈ 2*sqrt(R*L/(R+L)^3)
        // For A~0.8*0.8=0.64, R/L~5:1, sigma_stat ≈ 1.0/sqrt(N) roughly
        // More precisely: sigma_A = sqrt((1-A^2)/N) for large N
        Double_t dA_el = TMath::Abs(asym_theta_el[ie] - asym_meyer_el[ie]);
        Double_t dA_in = TMath::Abs(asym_theta_in[ie] - asym_meyer_in[ie]);

        // N for pull = k: N = (k / dA)^2 * (1-A^2)
        // Use A = average of the two
        Double_t A_el = 0.5*(asym_theta_el[ie] + asym_meyer_el[ie]);
        Double_t A_in = 0.5*(asym_theta_in[ie] + asym_meyer_in[ie]);
        Double_t var_el = 1.0 - A_el*A_el;
        Double_t var_in = 1.0 - A_in*A_in;

        auto N_for_pull = [](Double_t dA, Double_t var, Double_t k) -> Long64_t {
            if (dA <= 0) return -1;
            return (Long64_t)((k / dA) * (k / dA) * var);
        };

        Long64_t el1 = N_for_pull(dA_el, var_el, 1.);
        Long64_t el2 = N_for_pull(dA_el, var_el, 2.);
        Long64_t el3 = N_for_pull(dA_el, var_el, 3.);
        Long64_t in1 = N_for_pull(dA_in, var_in, 1.);
        Long64_t in2 = N_for_pull(dA_in, var_in, 2.);
        Long64_t in3 = N_for_pull(dA_in, var_in, 3.);

        TString note = (ekin_arr[ie] == 200.) ? "CSV valid" : "CSV=200MeV approx";
        cout << Form("  %-10.0f  %6lld / %6lld / %6lld   %6lld / %6lld / %6lld   %s",
                     ekin_arr[ie], el1, el2, el3, in1, in2, in3, note.Data()) << endl;
    }
    cout << "════════════════════════════════════════════════════════════════\n" << endl;

    // ================================================================
    // Plots: 2×2 canvas
    //   [0,0] Elastic XS ratio vs energy
    //   [0,1] Elastic ΔA_N vs energy
    //   [1,0] Inelastic XS ratio vs energy
    //   [1,1] Inelastic ΔA_N vs energy
    // ================================================================
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c_sys", "Systematic: THETA_CM vs Meyer", 1400, 900);
    c->Divide(2, 2);

    TGraph* g_xs_ratio_el = new TGraph(nE, ekin_arr, xs_ratio_el);
    TGraph* g_an_diff_el  = new TGraph(nE, ekin_arr, an_diff_el);
    TGraph* g_xs_ratio_in = new TGraph(nE, ekin_arr, xs_ratio_in);
    TGraph* g_an_diff_in  = new TGraph(nE, ekin_arr, an_diff_in);

    auto StyleSys = [](TGraph* g, Color_t c, Style_t m) {
        g->SetLineColor(c); g->SetLineWidth(2);
        g->SetMarkerColor(c); g->SetMarkerStyle(m); g->SetMarkerSize(1.4);
    };
    StyleSys(g_xs_ratio_el, kBlue+1,   20);
    StyleSys(g_an_diff_el,  kBlue+1,   20);
    StyleSys(g_xs_ratio_in, kOrange+7, 21);
    StyleSys(g_an_diff_in,  kOrange+7, 21);

    auto DrawUnity = [](Double_t xlo, Double_t xhi) {
        TLine* l = new TLine(xlo, 1.0, xhi, 1.0);
        l->SetLineColor(kGreen+2); l->SetLineWidth(2); l->SetLineStyle(2);
        l->Draw();
    };
    auto DrawZero = [](Double_t xlo, Double_t xhi) {
        TLine* l = new TLine(xlo, 0.0, xhi, 0.0);
        l->SetLineColor(kGreen+2); l->SetLineWidth(2); l->SetLineStyle(2);
        l->Draw();
    };

    // Panel 0,0: Elastic XS ratio
    c->cd(1); gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    {
        Double_t ylo = 0.8, yhi = 1.4;
        TH1F* fr = gPad->DrawFrame(155, ylo, 205, yhi);
        fr->SetTitle("Elastic: XS ratio Meyer/OptModel;E_{kin} [MeV];XS ratio (Meyer / Opt.Model)");
        fr->GetYaxis()->SetTitleOffset(1.4);
        DrawUnity(155, 205);
        g_xs_ratio_el->Draw("LP SAME");
        TLatex* l = new TLatex(); l->SetNDC(); l->SetTextSize(0.038);
        for (Int_t ie = 0; ie < nE; ie++)
            l->DrawLatex(0.18 + ie*0.15, 0.20,
                         Form("%.3f", xs_ratio_el[ie]));
    }

    // Panel 0,1: Elastic ΔA_N
    c->cd(2); gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    {
        Double_t amax = 0;
        for (Int_t ie = 0; ie < nE; ie++) amax = TMath::Max(amax, TMath::Abs(an_diff_el[ie]));
        amax = TMath::Max(amax*1.5, 0.05);
        TH1F* fr = gPad->DrawFrame(155, -amax, 205, amax);
        fr->SetTitle("Elastic: #DeltaA_{N} = A_{N}(Meyer) #minus A_{N}(OptModel);E_{kin} [MeV];#Delta A_{N}");
        fr->GetYaxis()->SetTitleOffset(1.4);
        DrawZero(155, 205);
        g_an_diff_el->Draw("LP SAME");
        TLatex* l = new TLatex(); l->SetNDC(); l->SetTextSize(0.038);
        for (Int_t ie = 0; ie < nE; ie++)
            l->DrawLatex(0.18 + ie*0.15, 0.20,
                         Form("%+.3f", an_diff_el[ie]));
    }

    // Panel 1,0: Inelastic XS ratio
    c->cd(3); gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    {
        Double_t ylo = 0.8, yhi = 1.4;
        TH1F* fr = gPad->DrawFrame(155, ylo, 205, yhi);
        fr->SetTitle("Inelastic: XS ratio Meyer/CSV;E_{kin} [MeV];XS ratio (Meyer / CSV)");
        fr->GetYaxis()->SetTitleOffset(1.4);
        DrawUnity(155, 205);
        g_xs_ratio_in->Draw("LP SAME");
        TLatex* l = new TLatex(); l->SetNDC(); l->SetTextSize(0.038);
        for (Int_t ie = 0; ie < nE; ie++)
            l->DrawLatex(0.18 + ie*0.15, 0.20,
                         Form("%.3f", xs_ratio_in[ie]));
        // Note about CSV energy validity
        TLatex* note = new TLatex();
        note->SetNDC(); note->SetTextSize(0.032); note->SetTextColor(kGray+2);
        note->DrawLatex(0.15, 0.88, "CSV data valid only at 200 MeV");
    }

    // Panel 1,1: Inelastic ΔA_N
    c->cd(4); gPad->SetGrid(); gPad->SetLeftMargin(0.14);
    {
        Double_t amax = 0;
        for (Int_t ie = 0; ie < nE; ie++) amax = TMath::Max(amax, TMath::Abs(an_diff_in[ie]));
        amax = TMath::Max(amax*1.5, 0.05);
        TH1F* fr = gPad->DrawFrame(155, -amax, 205, amax);
        fr->SetTitle("Inelastic: #DeltaA_{N} = A_{N}(Meyer) #minus A_{N}(CSV);E_{kin} [MeV];#Delta A_{N}");
        fr->GetYaxis()->SetTitleOffset(1.4);
        DrawZero(155, 205);
        g_an_diff_in->Draw("LP SAME");
        TLatex* l = new TLatex(); l->SetNDC(); l->SetTextSize(0.038);
        for (Int_t ie = 0; ie < nE; ie++)
            l->DrawLatex(0.18 + ie*0.15, 0.20,
                         Form("%+.3f", an_diff_in[ie]));
    }

    c->Update();
    c->SaveAs("Systematics_ThetaCM_vs_Meyer.pdf");
    cout << "Saved: Systematics_ThetaCM_vs_Meyer.pdf\n" << endl;
}