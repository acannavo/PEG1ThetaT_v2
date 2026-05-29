// CompareTheta_vs_Meyer.C
// Compare THETA_CM_SAMPLING vs T_SAMPLING_MEYER for elastic and inelastic
//
// Usage — generate files first, then run comparison:
//
//   Step 1: generate all four files (run once, takes ~1-2 min at 500k events)
//     root [0] GenerateAllForComparison(180.0, 500000, 0.80, +1)
//
//   Step 2: run elastic comparison
//     root [0] CompareElastic("pC_Elas_180MeV_MT_P80_SpinUp_16p2_THETA.root",
//                             "pC_Elas_180MeV_MT_P80_SpinUp_16p2_MEYER.root",
//                             180.0, 0.80, +1)
//
//   Step 3: run inelastic comparison
//     root [0] CompareInelastic("pC_Inel443_180MeV_MT_P80_SpinUp_16p2_THETA.root",
//                               "pC_Inel443_180MeV_MT_P80_SpinUp_16p2_MEYER.root",
//                               180.0, 0.80, +1)
//
//   Or run everything at once:
//     root [0] CompareAll(180.0, 500000, 0.80, +1)

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include <sstream>

// ====================================================================
// Helper: style a histogram pair consistently
// ====================================================================
void StylePair(TH1D* h_theta, TH1D* h_meyer,
               const char* title,
               const char* xtitle, const char* ytitle = "Events")
{
    h_theta->SetTitle(Form("%s;%s;%s", title, xtitle, ytitle));
    h_theta->SetLineColor(kAzure+1);
    h_theta->SetLineWidth(2);
    h_theta->SetFillColorAlpha(kAzure+1, 0.15);
    h_theta->SetFillStyle(1001);

    h_meyer->SetLineColor(kOrange+7);
    h_meyer->SetLineWidth(2);
    h_meyer->SetLineStyle(2);
    h_meyer->SetFillStyle(0);
}

// ====================================================================
// Helper: draw ratio panel (meyer/theta) with unity line
// Draws into the currently active pad — caller must cd() first
// ====================================================================
void DrawRatioPanel(TH1D* h_num, TH1D* h_den,
                    const char* xtitle,
                    Double_t ylo = 0.7, Double_t yhi = 1.3)
{
    TH1D* h_ratio = (TH1D*)h_num->Clone(
        Form("h_ratio_%s", h_num->GetName()));
    h_ratio->SetDirectory(0);
    h_ratio->Divide(h_den);
    h_ratio->SetTitle(Form(";%s;Meyer / #theta_{CM}", xtitle));
    h_ratio->SetLineColor(kOrange+7);
    h_ratio->SetLineWidth(2);
    h_ratio->SetMarkerColor(kOrange+7);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerSize(0.7);
    h_ratio->GetYaxis()->SetRangeUser(ylo, yhi);
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.10);
    h_ratio->GetYaxis()->SetLabelSize(0.09);
    h_ratio->GetXaxis()->SetTitleSize(0.10);
    h_ratio->GetXaxis()->SetLabelSize(0.09);
    h_ratio->GetYaxis()->SetTitleOffset(0.45);
    h_ratio->Draw("P");

    TLine* unity = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0,
                             h_ratio->GetXaxis()->GetXmax(), 1.0);
    unity->SetLineColor(kAzure+1);
    unity->SetLineWidth(2);
    unity->SetLineStyle(2);
    unity->Draw();
}

// ====================================================================
// Helper: compute asymmetry and its statistical uncertainty
//   A = (N_R - N_L) / (N_R + N_L)
//   δA = 2 * sqrt(N_R*N_L) / (N_R+N_L)²   [binomial error]
// ====================================================================
void ComputeAsymmetry(TH1D* h_phi,
                      Double_t& asym, Double_t& asym_err,
                      Double_t& n_right, Double_t& n_left)
{
    // Standard left-right asymmetry for a polarimeter:
    // RIGHT = phi in [-90°, +90°]  (detector at phi = 0°)
    // LEFT  = phi in [-180°,-90°] U [+90°,+180°]  (detector at phi = 180°)
    //
    // This is more robust than counting in a narrow window because
    // the detector phi acceptance (±5 mrad) maps to only ~1 histogram bin
    // when the histogram covers ±180°. Integrating over the full hemisphere
    // captures all the asymmetry signal regardless of binning.
    //
    // For a distribution W(phi) = 1 + A*cos(phi), integrating over [-90,+90]
    // gives N_R = N_tot/2 * (1 + 2A/pi) and over the complementary half gives
    // N_L = N_tot/2 * (1 - 2A/pi), so A_LR = pi/2 * A_phys.
    // We report A_LR directly (not corrected to A_phys) for the comparison.

    n_right = h_phi->Integral(h_phi->FindBin(-90.0 + 1e-6),
                              h_phi->FindBin( 90.0 - 1e-6));
    n_left  = h_phi->Integral(1,
                              h_phi->FindBin(-90.0 - 1e-6));
    n_left += h_phi->Integral(h_phi->FindBin( 90.0 + 1e-6),
                              h_phi->GetNbinsX());

    Double_t total = n_right + n_left;
    if (total <= 0.) { asym = 0.; asym_err = 0.; return; }

    asym     = (n_right - n_left) / total;
    // Binomial error: δA = sqrt(4*N_R*N_L/N_tot^3)
    asym_err = 2.0 * TMath::Sqrt(n_right * n_left) / (total * TMath::Sqrt(total));
}

// ====================================================================
// Helper: build a split pad (main + ratio) within a sub-canvas cell.
// Returns {pad_main, pad_ratio}. Caller does cd(pad_main) to draw.
// ====================================================================
// NOTE: ROOT does not support nested TPad::Divide cleanly inside a
// pre-divided TCanvas. Instead we use absolute coordinates.
// Caller must SetFillStyle(0) on both pads to avoid white rectangles.

// ====================================================================
// Helper: open file + tree with clean error messages
// Returns nullptr on failure (caller must check)
// ====================================================================
TTree* OpenTree(const char* fname, const char* treename = "data")
{
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) {
        cerr << "ERROR: Cannot open file: " << fname << endl;
        return nullptr;
    }
    TTree* t = (TTree*)f->Get(treename);
    if (!t) {
        cerr << "ERROR: No tree '" << treename << "' in " << fname << endl;
        f->Close();
        return nullptr;
    }
    return t;   // File stays open (ROOT tree needs it)
}

// ====================================================================
// Helper: summary pave text — fills pad 9
// ====================================================================
void DrawSummaryPave(
    const char* reaction_label,
    Double_t    energy,
    Double_t    polarization,
    Int_t       spin_state,
    Long64_t    n_theta,
    Long64_t    n_meyer,
    Double_t    asym_theta,  Double_t asym_theta_err,
    Double_t    asym_meyer,  Double_t asym_meyer_err,
    Double_t    ks_prob_theta, Double_t ks_prob_phi,
    Double_t    ks_prob_p)
{
    gPad->SetFillColor(kWhite);
    gPad->SetFrameFillColor(kWhite);

    TPaveText* pt = new TPaveText(0.05, 0.02, 0.95, 0.98, "NDC");
    pt->SetFillColor(kWhite);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.045);

    pt->AddText(Form("=== %s  |  E_{kin} = %.0f MeV ===", reaction_label, energy));
    pt->AddText(" ");
    pt->AddText(Form("Polarization: P = %.0f%%   Spin: %s",
                     polarization*100,
                     spin_state > 0 ? "UP (+1)" : "DOWN (-1)"));
    pt->AddText(Form("Events  #theta_{CM}: %lld     Meyer: %lld", n_theta, n_meyer));
    pt->AddText(" ");
    pt->AddText("─── Asymmetry  A = (R-L)/(R+L) ───────────────");
    pt->AddText(Form("  #theta_{CM}: %.4f #pm %.4f", asym_theta, asym_theta_err));
    pt->AddText(Form("  Meyer:       %.4f #pm %.4f", asym_meyer, asym_meyer_err));
    Double_t pull = (asym_theta_err > 0 || asym_meyer_err > 0)
                    ? (asym_theta - asym_meyer) /
                      TMath::Sqrt(asym_theta_err*asym_theta_err +
                                  asym_meyer_err*asym_meyer_err)
                    : 0.;
    pt->AddText(Form("  Pull: %.2f #sigma", pull));
    pt->AddText(" ");
    pt->AddText("─── KS compatibility tests ────────────────────");
    pt->AddText(Form("  #theta_{lab} :  p = %.4f %s",
                     ks_prob_theta,
                     ks_prob_theta > 0.05 ? "  #color[416]{OK}" : "  #color[632]{LOW}"));
    pt->AddText(Form("  #phi_{lab}   :  p = %.4f %s",
                     ks_prob_phi,
                     ks_prob_phi > 0.05 ? "  #color[416]{OK}" : "  #color[632]{LOW}"));
    pt->AddText(Form("  |p|          :  p = %.4f %s",
                     ks_prob_p,
                     ks_prob_p > 0.05 ? "  #color[416]{OK}" : "  #color[632]{LOW}"));
    pt->AddText(" ");

    Bool_t pass = (TMath::Abs(pull) < 2.0 &&
                   ks_prob_theta > 0.05 &&
                   ks_prob_phi   > 0.05 &&
                   ks_prob_p     > 0.05);

    TText* verdict = pt->AddText(pass ? "VERDICT:  #color[416]{GOOD AGREEMENT}"
                                      : "VERDICT:  #color[632]{DISTRIBUTIONS DIFFER}");
    verdict->SetTextSize(0.055);

    pt->Draw();
}

// ====================================================================
// ELASTIC COMPARISON
// ====================================================================
void CompareElastic(const char* fname_theta,
                    const char* fname_meyer,
                    Double_t    energy      = 180.0,
                    Double_t    polarization = 0.80,
                    Int_t       spin_state  = +1)
{
    // Unique suffix per call so ROOT histogram names never collide
    // when this function is called more than once in the same session.
    static Int_t elas_call = 0;
    TString sfx = Form("_e%d", elas_call++);

    cout << "\n╔══════════════════════════════════════════════════════════╗" << endl;
    cout << "║  ELASTIC: #theta_{CM} vs Meyer Sampling Comparison        ║" << endl;
    cout << "╚══════════════════════════════════════════════════════════╝\n" << endl;

    TTree* tree_th = OpenTree(fname_theta);
    TTree* tree_me = OpenTree(fname_meyer);
    if (!tree_th || !tree_me) return;

    Long64_t n_theta = tree_th->GetEntries();
    Long64_t n_meyer = tree_me->GetEntries();
    cout << Form("  Theta file:  %s  (%lld events)", fname_theta, n_theta) << endl;
    cout << Form("  Meyer file:  %s  (%lld events)", fname_meyer, n_meyer) << endl;

    // ----------------------------------------------------------------
    // Global style
    // ----------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(0.055);
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetFrameLineWidth(1);

    // ----------------------------------------------------------------
    // Canvas: 4 columns × 3 rows
    //   Row 1: theta_lab, phi_lab, |p|_proton, p_T proton
    //   Row 2: ratio panels for the above four
    //   Row 3: theta_C, 2D theta-phi (theta), 2D theta-phi (Meyer), summary
    // ----------------------------------------------------------------
    TCanvas* c = new TCanvas(("c_elas"+sfx).Data(), "Elastic: #theta_{CM} vs Meyer", 2000, 1400);
    c->SetFillColor(kWhite);

    // We use a 4x3 divide with manual sub-pad splitting for ratio rows
    // Layout: 12 cells in 4 columns × 3 rows
    // Cells 1-4:  main distributions
    // Cells 5-8:  ratio panels  (shorter height)
    // Cells 9-12: theta_C, 2D_theta, 2D_meyer, summary

    // Manual pad layout for split main+ratio per column
    // Each column: main pad covers y=[0.43,1.0], ratio pad y=[0.30,0.43]
    // Row 3 covers y=[0.00,0.30]
    const Int_t  ncol      = 4;
    const Double_t x_lo[4] = {0.00, 0.25, 0.50, 0.75};
    const Double_t x_hi[4] = {0.25, 0.50, 0.75, 1.00};
    const Double_t y_main_lo = 0.43, y_main_hi = 1.00;
    const Double_t y_rat_lo  = 0.30, y_rat_hi  = 0.43;
    const Double_t y_bot_lo  = 0.00, y_bot_hi  = 0.30;

    TPad* pad_main[4];
    TPad* pad_rat[4];
    TPad* pad_bot[4];

    for (Int_t i = 0; i < ncol; i++) {
        pad_main[i] = new TPad(Form("pm%d",i), "", x_lo[i], y_main_lo, x_hi[i], y_main_hi);
        pad_main[i]->SetFillColor(kWhite);
        pad_main[i]->SetLeftMargin(0.16); pad_main[i]->SetRightMargin(0.03);
        pad_main[i]->SetTopMargin(0.10);  pad_main[i]->SetBottomMargin(0.01);
        pad_main[i]->Draw();

        pad_rat[i] = new TPad(Form("pr%d",i), "", x_lo[i], y_rat_lo, x_hi[i], y_rat_hi);
        pad_rat[i]->SetFillColor(kWhite);
        pad_rat[i]->SetLeftMargin(0.16); pad_rat[i]->SetRightMargin(0.03);
        pad_rat[i]->SetTopMargin(0.01);  pad_rat[i]->SetBottomMargin(0.35);
        pad_rat[i]->Draw();

        pad_bot[i] = new TPad(Form("pb%d",i), "", x_lo[i], y_bot_lo, x_hi[i], y_bot_hi);
        pad_bot[i]->SetFillColor(kWhite);
        pad_bot[i]->SetLeftMargin(0.16); pad_bot[i]->SetRightMargin(0.03);
        if (i < 2) { pad_bot[i]->SetTopMargin(0.08); pad_bot[i]->SetBottomMargin(0.18); }
        else        { pad_bot[i]->SetTopMargin(0.08); pad_bot[i]->SetBottomMargin(0.05); }
        pad_bot[i]->Draw();
    }

    // ================================================================
    // Detector window [degrees] — used for histogram ranges
    // ================================================================
    Double_t th_ctr = DETECTOR_THETA_CENTER;            // e.g. 16.2
    Double_t th_win = DETECTOR_THETA_WINDOW * TMath::RadToDeg();  // e.g. 0.286
    // Use the exact detector acceptance so the KS test compares distributions
    // WITHIN the detector window, not the sampling domain padding.
    // (THETA_CM pads the CM range by 10%, which shifts some events outside
    // the acceptance; those are rejected by the lab cuts and irrelevant here.)
    Double_t th_win_deg = DETECTOR_THETA_WINDOW * TMath::RadToDeg();
    Double_t th_lo  = th_ctr - th_win_deg * 1.05;
    Double_t th_hi  = th_ctr + th_win_deg * 1.05;

    // ================================================================
    // COL 0: Proton theta_lab
    // ================================================================
    TH1D* h_th_th = new TH1D(("h_th_th"+sfx).Data(), "", 80, th_lo, th_hi);
    TH1D* h_th_me = new TH1D(("h_th_me"+sfx).Data(), "", 80, th_lo, th_hi);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    tree_th->Draw(Form("Particles.Theta()*TMath::RadToDeg()>>%s", ("h_th_th"+sfx).Data()), "Particles.ID()==14", "goff");
    tree_me->Draw(Form("Particles.Theta()*TMath::RadToDeg()>>%s", ("h_th_me"+sfx).Data()), "Particles.ID()==14", "goff");

    // Normalize per event so different N doesn't distort comparison
    if (h_th_th->Integral() > 0) h_th_th->Scale(1.0/h_th_th->Integral());
    if (h_th_me->Integral() > 0) h_th_me->Scale(1.0/h_th_me->Integral());

    StylePair(h_th_th, h_th_me,
              Form("Proton #theta_{lab}  (E=%.0f MeV)",energy),
              "#theta_{lab} [deg]", "Norm. events");

    Double_t ks_theta = h_th_th->KolmogorovTest(h_th_me);

    pad_main[0]->cd();
    Double_t ymax0 = 1.3 * TMath::Max(h_th_th->GetMaximum(), h_th_me->GetMaximum());
    h_th_th->GetYaxis()->SetRangeUser(0, ymax0);
    h_th_th->Draw("HIST");
    h_th_me->Draw("HIST SAME");
    TLegend* leg0 = new TLegend(0.18, 0.76, 0.72, 0.90);
    leg0->SetBorderSize(0); leg0->SetFillStyle(0); leg0->SetTextSize(0.048);
    leg0->AddEntry(h_th_th, "#theta_{CM} sampling", "lf");
    leg0->AddEntry(h_th_me, "Meyer sampling",        "l");
    leg0->Draw();
    TLatex* lat = new TLatex();
    lat->SetNDC(); lat->SetTextSize(0.042); lat->SetTextColor(kGray+2);
    lat->DrawLatex(0.18, 0.70, Form("KS p = %.3f", ks_theta));

    pad_rat[0]->cd();
    DrawRatioPanel(h_th_me, h_th_th, "#theta_{lab} [deg]");

    // ================================================================
    // COL 1: Proton phi_lab
    // ================================================================
    TH1D* h_ph_th = new TH1D(("h_ph_th"+sfx).Data(), "", 72, -180, 180);
    TH1D* h_ph_me = new TH1D(("h_ph_me"+sfx).Data(), "", 72, -180, 180);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    tree_th->Draw(Form("Particles.Phi()*TMath::RadToDeg()>>%s", ("h_ph_th"+sfx).Data()), "Particles.ID()==14", "goff");
    tree_me->Draw(Form("Particles.Phi()*TMath::RadToDeg()>>%s", ("h_ph_me"+sfx).Data()), "Particles.ID()==14", "goff");

    // Clone BEFORE normalizing — ComputeAsymmetry needs raw event counts
    TH1D* h_ph_th_raw = (TH1D*)h_ph_th->Clone(("h_ph_th_raw"+sfx).Data());
    TH1D* h_ph_me_raw = (TH1D*)h_ph_me->Clone(("h_ph_me_raw"+sfx).Data());
    h_ph_th_raw->SetDirectory(0); h_ph_me_raw->SetDirectory(0);

    if (h_ph_th->Integral() > 0) h_ph_th->Scale(1.0/h_ph_th->Integral());
    if (h_ph_me->Integral() > 0) h_ph_me->Scale(1.0/h_ph_me->Integral());

    StylePair(h_ph_th, h_ph_me,
              "#phi_{lab} distribution  (#color[416]{polarisation signal})",
              "#phi_{lab} [deg]", "Norm. events");

    Double_t ks_phi = h_ph_th->KolmogorovTest(h_ph_me);

    // Asymmetries — use the raw (pre-normalization) clones
    // histogram lives in gDirectory so tree->Draw(">>name") can find it
    // Raw clone already has event counts — no re-scaling needed

    Double_t asym_th, asym_th_err, nr_th, nl_th;
    Double_t asym_me, asym_me_err, nr_me, nl_me;
    ComputeAsymmetry(h_ph_th_raw, asym_th, asym_th_err, nr_th, nl_th);
    ComputeAsymmetry(h_ph_me_raw, asym_me, asym_me_err, nr_me, nl_me);

    pad_main[1]->cd();
    Double_t ymax1 = 1.3 * TMath::Max(h_ph_th->GetMaximum(), h_ph_me->GetMaximum());
    h_ph_th->GetYaxis()->SetRangeUser(0, ymax1);
    h_ph_th->Draw("HIST");
    h_ph_me->Draw("HIST SAME");
    TLatex* lat1 = new TLatex();
    lat1->SetNDC(); lat1->SetTextSize(0.042); lat1->SetTextColor(kGray+2);
    lat1->DrawLatex(0.18, 0.85, Form("A(#theta_{CM}) = %.4f #pm %.4f", asym_th, asym_th_err));
    lat1->DrawLatex(0.18, 0.77, Form("A(Meyer)       = %.4f #pm %.4f", asym_me, asym_me_err));
    lat1->DrawLatex(0.18, 0.69, Form("KS p = %.3f", ks_phi));
    TLegend* leg1b = new TLegend(0.55, 0.76, 0.96, 0.90);
    leg1b->SetBorderSize(0); leg1b->SetFillStyle(0); leg1b->SetTextSize(0.044);
    leg1b->AddEntry(h_ph_th, "#theta_{CM}", "lf");
    leg1b->AddEntry(h_ph_me, "Meyer",        "l");
    leg1b->Draw();

    pad_rat[1]->cd();
    DrawRatioPanel(h_ph_me, h_ph_th, "#phi_{lab} [deg]");

    // ================================================================
    // COL 2: Proton |p|
    // ================================================================
    // Momentum range for elastic: p ~ sqrt(2*mp*T + T²) with T ~ ekin
    Double_t ekinGeV = energy / 1000.;
    Double_t p_beam  = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    Double_t p_lo    = 0.85 * p_beam;
    Double_t p_hi    = 1.05 * p_beam;

    TH1D* h_p_th = new TH1D(("h_p_th"+sfx).Data(), "", 80, p_lo, p_hi);
    TH1D* h_p_me = new TH1D(("h_p_me"+sfx).Data(), "", 80, p_lo, p_hi);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    tree_th->Draw(Form("Particles.P()>>%s", ("h_p_th"+sfx).Data()), "Particles.ID()==14", "goff");
    tree_me->Draw(Form("Particles.P()>>%s", ("h_p_me"+sfx).Data()), "Particles.ID()==14", "goff");

    if (h_p_th->Integral() > 0) h_p_th->Scale(1.0/h_p_th->Integral());
    if (h_p_me->Integral() > 0) h_p_me->Scale(1.0/h_p_me->Integral());

    StylePair(h_p_th, h_p_me,
              "Proton |p|  (#color[616]{kinematics check})",
              "|p| [GeV/c]", "Norm. events");
    Double_t ks_p = h_p_th->KolmogorovTest(h_p_me);

    pad_main[2]->cd();
    Double_t ymax2 = 1.3 * TMath::Max(h_p_th->GetMaximum(), h_p_me->GetMaximum());
    h_p_th->GetYaxis()->SetRangeUser(0, ymax2);
    h_p_th->Draw("HIST");
    h_p_me->Draw("HIST SAME");
    TLatex* lat2 = new TLatex();
    lat2->SetNDC(); lat2->SetTextSize(0.042); lat2->SetTextColor(kGray+2);
    lat2->DrawLatex(0.18, 0.85, Form("Mean(#theta_{CM}) = %.4f GeV/c", h_p_th->GetMean()));
    lat2->DrawLatex(0.18, 0.77, Form("Mean(Meyer)       = %.4f GeV/c", h_p_me->GetMean()));
    lat2->DrawLatex(0.18, 0.69, Form("KS p = %.3f", ks_p));

    pad_rat[2]->cd();
    DrawRatioPanel(h_p_me, h_p_th, "|p| [GeV/c]");

    // ================================================================
    // COL 3: Proton p_T
    // ================================================================
    TH1D* h_pt_th = new TH1D(("h_pt_th"+sfx).Data(), "", 80, 0, 0.30);
    TH1D* h_pt_me = new TH1D(("h_pt_me"+sfx).Data(), "", 80, 0, 0.30);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    tree_th->Draw(Form("Particles.Perp()>>%s", ("h_pt_th"+sfx).Data()), "Particles.ID()==14", "goff");
    tree_me->Draw(Form("Particles.Perp()>>%s", ("h_pt_me"+sfx).Data()), "Particles.ID()==14", "goff");

    if (h_pt_th->Integral() > 0) h_pt_th->Scale(1.0/h_pt_th->Integral());
    if (h_pt_me->Integral() > 0) h_pt_me->Scale(1.0/h_pt_me->Integral());

    StylePair(h_pt_th, h_pt_me, "Proton p_{T}", "p_{T} [GeV/c]", "Norm. events");
    Double_t ks_pt = h_pt_th->KolmogorovTest(h_pt_me);

    pad_main[3]->cd();
    Double_t ymax3 = 1.3 * TMath::Max(h_pt_th->GetMaximum(), h_pt_me->GetMaximum());
    h_pt_th->GetYaxis()->SetRangeUser(0, ymax3);
    h_pt_th->Draw("HIST");
    h_pt_me->Draw("HIST SAME");
    TLatex* lat3 = new TLatex();
    lat3->SetNDC(); lat3->SetTextSize(0.042); lat3->SetTextColor(kGray+2);
    lat3->DrawLatex(0.18, 0.85, Form("KS p = %.3f", ks_pt));

    pad_rat[3]->cd();
    DrawRatioPanel(h_pt_me, h_pt_th, "p_{T} [GeV/c]");

    // ================================================================
    // BOT COL 0: Carbon theta_lab
    // ================================================================
    TH1D* h_thC_th = new TH1D(("h_thC_th"+sfx).Data(), "", 100, 70, 90);
    TH1D* h_thC_me = new TH1D(("h_thC_me"+sfx).Data(), "", 100, 70, 90);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    tree_th->Draw(Form("Particles.Theta()*TMath::RadToDeg()>>%s", ("h_thC_th"+sfx).Data()), "Particles.ID()==614", "goff");
    tree_me->Draw(Form("Particles.Theta()*TMath::RadToDeg()>>%s", ("h_thC_me"+sfx).Data()), "Particles.ID()==614", "goff");

    if (h_thC_th->Integral() > 0) h_thC_th->Scale(1.0/h_thC_th->Integral());
    if (h_thC_me->Integral() > 0) h_thC_me->Scale(1.0/h_thC_me->Integral());

    StylePair(h_thC_th, h_thC_me,
              "^{12}C recoil #theta_{lab}",
              "#theta_{lab}^{C} [deg]", "Norm. events");

    pad_bot[0]->cd();
    Double_t ymaxC = 1.3 * TMath::Max(h_thC_th->GetMaximum(), h_thC_me->GetMaximum());
    h_thC_th->GetYaxis()->SetRangeUser(0, ymaxC);
    h_thC_th->Draw("HIST");
    h_thC_me->Draw("HIST SAME");
    Double_t ks_thC = h_thC_th->KolmogorovTest(h_thC_me);
    TLatex* latC = new TLatex();
    latC->SetNDC(); latC->SetTextSize(0.07); latC->SetTextColor(kGray+2);
    latC->DrawLatex(0.18, 0.87, Form("KS p = %.3f", ks_thC));

    // ================================================================
    // BOT COL 1: Asymmetry bar chart (R-L)/(R+L) for both modes
    // ================================================================
    pad_bot[1]->cd();
    gPad->SetGridy();

    // Expected asymmetry: A = P * A_N  (rough estimate for label)
    // A_N at this theta from Meyer spline — evaluated at t corresponding
    // to the detector centre angle
    // Use the CM angle mid-point corresponding to the detector centre.
    // ComputeCMAngleRange is available from Kinematics.C; we call it silently
    // by redirecting stdout — or simply use the known mid-point from Stage 1
    // printout (~17.78° for 16.2° lab at 180 MeV).  To be fully general we
    // compute it properly:
    Double_t theta_cm_lo_tmp, theta_cm_hi_tmp;
    {
        // Suppress the ComputeCMAngleRange printout for this auxiliary call
        std::streambuf* oldCout = std::cout.rdbuf();
        std::ostringstream devnull;
        std::cout.rdbuf(devnull.rdbuf());
        ComputeCMAngleRange(energy,
                            DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW,
                            theta_cm_lo_tmp, theta_cm_hi_tmp);
        std::cout.rdbuf(oldCout);
    }
    Double_t theta_cm_ctr = 0.5 * (theta_cm_lo_tmp + theta_cm_hi_tmp);
    Double_t t_ctr    = ConvertThetaCMtoT(theta_cm_ctr, energy);
    Double_t AN_meyer = MeyerAP_Elastic(t_ctr, energy);
    Double_t A_exp    = polarization * AN_meyer;

    // Draw as a grouped bar chart with error bars
    Double_t x_pos[2]   = {1.0, 2.0};
    Double_t asym_vals[2]   = {asym_th, asym_me};
    Double_t asym_errs[2]   = {asym_th_err, asym_me_err};
    Double_t x_errs[2]  = {0., 0.};

    TGraphErrors* g_asym = new TGraphErrors(2, x_pos, asym_vals, x_errs, asym_errs);
    g_asym->SetTitle("Azimuthal Asymmetry;Sampling mode;A = (R-L)/(R+L)");
    g_asym->SetMarkerStyle(20);
    g_asym->SetMarkerSize(2.0);
    g_asym->SetMarkerColor(kBlack);
    g_asym->SetLineColor(kBlack);
    g_asym->SetLineWidth(2);
    g_asym->GetXaxis()->SetBinLabel(1, "#theta_{CM}");
    g_asym->GetYaxis()->SetRangeUser(
        TMath::Min(asym_th, asym_me) - 5.*TMath::Max(asym_th_err, asym_me_err),
        TMath::Max(asym_th, asym_me) + 5.*TMath::Max(asym_th_err, asym_me_err));
    g_asym->GetXaxis()->SetLimits(0, 3);
    g_asym->GetXaxis()->SetNdivisions(3);
    g_asym->Draw("AP");

    // Labels on points
    TLatex* lA = new TLatex();
    lA->SetTextSize(0.07); lA->SetTextAlign(22);
    lA->DrawLatex(1.0, asym_th - 3.*asym_th_err, "#theta_{CM}");
    lA->DrawLatex(2.0, asym_me - 3.*asym_me_err, "Meyer");

    // Expected value line
    TLine* lexp = new TLine(0, A_exp, 3, A_exp);
    lexp->SetLineColor(kGreen+2); lexp->SetLineWidth(2); lexp->SetLineStyle(2);
    lexp->Draw();
    TLatex* lAexp = new TLatex();
    lAexp->SetNDC(); lAexp->SetTextSize(0.065); lAexp->SetTextColor(kGreen+2);
    lAexp->DrawLatex(0.18, 0.87, Form("Expected: P#timesA_{N} = %.3f", A_exp));

    // ================================================================
    // BOT COL 2: 2D theta-phi for THETA mode
    // ================================================================
    pad_bot[2]->cd();
    gPad->SetRightMargin(0.15);

    TH2D* h2_th = new TH2D(("h2_th"+sfx).Data(),
        "#theta_{CM} sampling: #phi vs #theta;#phi_{lab} [deg];#theta_{lab} [deg]",
        36, -180, 180, 40, th_lo, th_hi);
    tree_th->Draw(
        Form("Particles.Theta()*TMath::RadToDeg():Particles.Phi()*TMath::RadToDeg()>>%s", ("h2_th"+sfx).Data()),
        "Particles.ID()==14", "goff");
    h2_th->Draw("COLZ");

    // ================================================================
    // BOT COL 3: 2D theta-phi for MEYER mode
    // ================================================================
    pad_bot[3]->cd();
    gPad->SetRightMargin(0.15);

    TH2D* h2_me = new TH2D(("h2_me"+sfx).Data(),
        "Meyer sampling: #phi vs #theta;#phi_{lab} [deg];#theta_{lab} [deg]",
        36, -180, 180, 40, th_lo, th_hi);
    tree_me->Draw(
        Form("Particles.Theta()*TMath::RadToDeg():Particles.Phi()*TMath::RadToDeg()>>%s", ("h2_me"+sfx).Data()),
        "Particles.ID()==14", "goff");
    h2_me->Draw("COLZ");

    // ================================================================
    // Console summary
    // ================================================================
    cout << "\n══════════════════════════════════════════════════════════" << endl;
    cout << Form("  ELASTIC SUMMARY  |  E = %.0f MeV  |  P = %.0f%%  |  Spin %s",
                 energy, polarization*100, spin_state>0?"UP":"DOWN") << endl;
    cout << "══════════════════════════════════════════════════════════" << endl;
    cout << Form("  %-20s  %12lld events", "#theta_CM sampling:", n_theta) << endl;
    cout << Form("  %-20s  %12lld events", "Meyer sampling:",     n_meyer) << endl;
    cout << endl;
    cout << "  KS compatibility tests:" << endl;
    cout << Form("    theta_lab : p = %.4f %s", ks_theta, ks_theta>0.05?"✓":"✗") << endl;
    cout << Form("    phi_lab   : p = %.4f %s", ks_phi,   ks_phi  >0.05?"✓":"✗") << endl;
    cout << Form("    |p|       : p = %.4f %s", ks_p,     ks_p    >0.05?"✓":"✗") << endl;
    cout << endl;
    cout << "  Azimuthal asymmetry A = (R-L)/(R+L):" << endl;
    cout << Form("    #theta_CM : %.4f +/- %.4f  (R=%.0f, L=%.0f)",
                 asym_th, asym_th_err, nr_th, nl_th) << endl;
    cout << Form("    Meyer     : %.4f +/- %.4f  (R=%.0f, L=%.0f)",
                 asym_me, asym_me_err, nr_me, nl_me) << endl;
    Double_t pull_elas = (asym_th_err>0||asym_me_err>0)
        ? (asym_th-asym_me)/TMath::Sqrt(asym_th_err*asym_th_err+asym_me_err*asym_me_err)
        : 0.;
    cout << Form("    Pull      : %.2f sigma", pull_elas) << endl;
    cout << "══════════════════════════════════════════════════════════\n" << endl;

    c->Update();

    TString outname = Form("Compare_Elastic_%.0fMeV_P%.0f_%s.pdf",
                           energy, polarization*100,
                           spin_state>0?"SpinUp":"SpinDown");
    c->SaveAs(outname.Data());
    cout << "Saved: " << outname << endl;
}

// ====================================================================
// INELASTIC COMPARISON
// ====================================================================
void CompareInelastic(const char* fname_theta,
                      const char* fname_meyer,
                      Double_t    energy       = 180.0,
                      Double_t    polarization = 0.80,
                      Int_t       spin_state   = +1)
{
    static Int_t inel_call = 0;
    TString sfx = Form("_i%d", inel_call++);

    cout << "\n╔══════════════════════════════════════════════════════════╗" << endl;
    cout << "║  INELASTIC: #theta_{CM} vs Meyer Sampling Comparison     ║" << endl;
    cout << "╚══════════════════════════════════════════════════════════╝\n" << endl;

    TTree* tree_th = OpenTree(fname_theta);
    TTree* tree_me = OpenTree(fname_meyer);
    if (!tree_th || !tree_me) return;

    Long64_t n_theta = tree_th->GetEntries();
    Long64_t n_meyer = tree_me->GetEntries();
    cout << Form("  Theta file:  %s  (%lld events)", fname_theta, n_theta) << endl;
    cout << Form("  Meyer file:  %s  (%lld events)", fname_meyer, n_meyer) << endl;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(0.055);

    Double_t th_ctr = DETECTOR_THETA_CENTER;
    Double_t th_win = DETECTOR_THETA_WINDOW * TMath::RadToDeg();
    // Use the exact detector acceptance so the KS test compares distributions
    // WITHIN the detector window, not the sampling domain padding.
    // (THETA_CM pads the CM range by 10%, which shifts some events outside
    // the acceptance; those are rejected by the lab cuts and irrelevant here.)
    Double_t th_win_deg = DETECTOR_THETA_WINDOW * TMath::RadToDeg();
    Double_t th_lo  = th_ctr - th_win_deg * 1.05;
    Double_t th_hi  = th_ctr + th_win_deg * 1.05;

    TCanvas* c = new TCanvas(("c_inel"+sfx).Data(), "Inelastic: #theta_{CM} vs Meyer", 2000, 1400);
    c->SetFillColor(kWhite);

    const Int_t    ncol      = 4;
    const Double_t x_lo[4]  = {0.00, 0.25, 0.50, 0.75};
    const Double_t x_hi[4]  = {0.25, 0.50, 0.75, 1.00};
    const Double_t y_main_lo = 0.43, y_main_hi = 1.00;
    const Double_t y_rat_lo  = 0.30, y_rat_hi  = 0.43;
    const Double_t y_bot_lo  = 0.00, y_bot_hi  = 0.30;

    TPad* pad_main[4];
    TPad* pad_rat[4];
    TPad* pad_bot[4];

    for (Int_t i = 0; i < ncol; i++) {
        pad_main[i] = new TPad(Form("ipm%d",i), "", x_lo[i], y_main_lo, x_hi[i], y_main_hi);
        pad_main[i]->SetFillColor(kWhite);
        pad_main[i]->SetLeftMargin(0.16); pad_main[i]->SetRightMargin(0.03);
        pad_main[i]->SetTopMargin(0.10);  pad_main[i]->SetBottomMargin(0.01);
        pad_main[i]->Draw();

        pad_rat[i] = new TPad(Form("ipr%d",i), "", x_lo[i], y_rat_lo, x_hi[i], y_rat_hi);
        pad_rat[i]->SetFillColor(kWhite);
        pad_rat[i]->SetLeftMargin(0.16); pad_rat[i]->SetRightMargin(0.03);
        pad_rat[i]->SetTopMargin(0.01);  pad_rat[i]->SetBottomMargin(0.35);
        pad_rat[i]->Draw();

        pad_bot[i] = new TPad(Form("ipb%d",i), "", x_lo[i], y_bot_lo, x_hi[i], y_bot_hi);
        pad_bot[i]->SetFillColor(kWhite);
        pad_bot[i]->SetLeftMargin(0.16); pad_bot[i]->SetRightMargin(0.03);
        pad_bot[i]->SetTopMargin(0.08);  pad_bot[i]->SetBottomMargin(0.18);
        pad_bot[i]->Draw();
    }

    // ================================================================
    // COL 0: Proton theta_lab
    // ================================================================
    TH1D* h_th_th = new TH1D(("ih_th_th"+sfx).Data(), "", 80, th_lo, th_hi);
    TH1D* h_th_me = new TH1D(("ih_th_me"+sfx).Data(), "", 80, th_lo, th_hi);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    tree_th->Draw(Form("Particles.Theta()*TMath::RadToDeg()>>%s", ("ih_th_th"+sfx).Data()), "Particles.ID()==14", "goff");
    tree_me->Draw(Form("Particles.Theta()*TMath::RadToDeg()>>%s", ("ih_th_me"+sfx).Data()), "Particles.ID()==14", "goff");

    if (h_th_th->Integral() > 0) h_th_th->Scale(1.0/h_th_th->Integral());
    if (h_th_me->Integral() > 0) h_th_me->Scale(1.0/h_th_me->Integral());

    StylePair(h_th_th, h_th_me,
              Form("Proton #theta_{lab}  [inel, E=%.0f MeV]", energy),
              "#theta_{lab} [deg]", "Norm. events");
    Double_t ks_theta = h_th_th->KolmogorovTest(h_th_me);

    pad_main[0]->cd();
    h_th_th->GetYaxis()->SetRangeUser(0, 1.3*TMath::Max(h_th_th->GetMaximum(),h_th_me->GetMaximum()));
    h_th_th->Draw("HIST");
    h_th_me->Draw("HIST SAME");
    TLegend* leg0i = new TLegend(0.18, 0.76, 0.72, 0.90);
    leg0i->SetBorderSize(0); leg0i->SetFillStyle(0); leg0i->SetTextSize(0.048);
    leg0i->AddEntry(h_th_th, "#theta_{CM} sampling", "lf");
    leg0i->AddEntry(h_th_me, "Meyer sampling",        "l");
    leg0i->Draw();
    TLatex* lati = new TLatex();
    lati->SetNDC(); lati->SetTextSize(0.042); lati->SetTextColor(kGray+2);
    lati->DrawLatex(0.18, 0.70, Form("KS p = %.3f", ks_theta));

    pad_rat[0]->cd();
    DrawRatioPanel(h_th_me, h_th_th, "#theta_{lab} [deg]");

    // ================================================================
    // COL 1: Proton phi_lab + asymmetry
    // ================================================================
    TH1D* h_ph_th = new TH1D(("ih_ph_th"+sfx).Data(), "", 72, -180, 180);
    TH1D* h_ph_me = new TH1D(("ih_ph_me"+sfx).Data(), "", 72, -180, 180);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    tree_th->Draw(Form("Particles.Phi()*TMath::RadToDeg()>>%s", ("ih_ph_th"+sfx).Data()), "Particles.ID()==14", "goff");
    tree_me->Draw(Form("Particles.Phi()*TMath::RadToDeg()>>%s", ("ih_ph_me"+sfx).Data()), "Particles.ID()==14", "goff");

    TH1D* h_ph_th_raw = (TH1D*)h_ph_th->Clone(("ih_ph_th_raw"+sfx).Data());
    TH1D* h_ph_me_raw = (TH1D*)h_ph_me->Clone(("ih_ph_me_raw"+sfx).Data());
    h_ph_th_raw->SetDirectory(0); h_ph_me_raw->SetDirectory(0);

    Double_t asym_th, asym_th_err, nr_th, nl_th;
    Double_t asym_me, asym_me_err, nr_me, nl_me;
    ComputeAsymmetry(h_ph_th_raw, asym_th, asym_th_err, nr_th, nl_th);
    ComputeAsymmetry(h_ph_me_raw, asym_me, asym_me_err, nr_me, nl_me);

    if (h_ph_th->Integral() > 0) h_ph_th->Scale(1.0/h_ph_th->Integral());
    if (h_ph_me->Integral() > 0) h_ph_me->Scale(1.0/h_ph_me->Integral());

    StylePair(h_ph_th, h_ph_me,
              "#phi_{lab}  [inel, #color[416]{polarisation signal}]",
              "#phi_{lab} [deg]", "Norm. events");
    Double_t ks_phi = h_ph_th->KolmogorovTest(h_ph_me);

    pad_main[1]->cd();
    h_ph_th->GetYaxis()->SetRangeUser(0, 1.3*TMath::Max(h_ph_th->GetMaximum(),h_ph_me->GetMaximum()));
    h_ph_th->Draw("HIST");
    h_ph_me->Draw("HIST SAME");
    TLatex* lat1i = new TLatex();
    lat1i->SetNDC(); lat1i->SetTextSize(0.042); lat1i->SetTextColor(kGray+2);
    lat1i->DrawLatex(0.18, 0.85, Form("A(#theta_{CM}) = %.4f #pm %.4f", asym_th, asym_th_err));
    lat1i->DrawLatex(0.18, 0.77, Form("A(Meyer)       = %.4f #pm %.4f", asym_me, asym_me_err));
    lat1i->DrawLatex(0.18, 0.69, Form("KS p = %.3f", ks_phi));

    pad_rat[1]->cd();
    DrawRatioPanel(h_ph_me, h_ph_th, "#phi_{lab} [deg]");

    // ================================================================
    // COL 2: Proton |p|  (inelastic: slightly lower than elastic)
    // ================================================================
    Double_t ekinGeV = energy / 1000.;
    Double_t p_beam  = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    // Inelastic proton loses ~4.43 MeV → momentum slightly lower
    Double_t p_lo    = 0.80 * p_beam;
    Double_t p_hi    = 1.02 * p_beam;

    TH1D* h_p_th = new TH1D(("ih_p_th"+sfx).Data(), "", 80, p_lo, p_hi);
    TH1D* h_p_me = new TH1D(("ih_p_me"+sfx).Data(), "", 80, p_lo, p_hi);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    tree_th->Draw(Form("Particles.P()>>%s", ("ih_p_th"+sfx).Data()), "Particles.ID()==14", "goff");
    tree_me->Draw(Form("Particles.P()>>%s", ("ih_p_me"+sfx).Data()), "Particles.ID()==14", "goff");

    if (h_p_th->Integral() > 0) h_p_th->Scale(1.0/h_p_th->Integral());
    if (h_p_me->Integral() > 0) h_p_me->Scale(1.0/h_p_me->Integral());

    StylePair(h_p_th, h_p_me, "Proton |p|  [inel]", "|p| [GeV/c]", "Norm. events");
    Double_t ks_p = h_p_th->KolmogorovTest(h_p_me);

    pad_main[2]->cd();
    h_p_th->GetYaxis()->SetRangeUser(0, 1.3*TMath::Max(h_p_th->GetMaximum(),h_p_me->GetMaximum()));
    h_p_th->Draw("HIST");
    h_p_me->Draw("HIST SAME");
    TLatex* lat2i = new TLatex();
    lat2i->SetNDC(); lat2i->SetTextSize(0.042); lat2i->SetTextColor(kGray+2);
    lat2i->DrawLatex(0.18, 0.85, Form("KS p = %.3f", ks_p));

    pad_rat[2]->cd();
    DrawRatioPanel(h_p_me, h_p_th, "|p| [GeV/c]");

    // ================================================================
    // COL 3: Carbon internal energy (E_internal = mass - mC)
    // KEY diagnostic for inelastic: should peak at 4.43 MeV for both
    // ================================================================
    TH1D* h_Eint_th = new TH1D(("ih_Eint_th"+sfx).Data(), "", 80, 0, 10);
    TH1D* h_Eint_me = new TH1D(("ih_Eint_me"+sfx).Data(), "", 80, 0, 10);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    // E_internal = sqrt(E²-p²) - mC_ground = M_Carbon - mC
    // In PParticle: M() gives the invariant mass in GeV
    tree_th->Draw(Form("(Particles.M()-0.11177)*1000>>%s", ("ih_Eint_th"+sfx).Data()), "Particles.ID()==614", "goff");
    tree_me->Draw(Form("(Particles.M()-0.11177)*1000>>%s", ("ih_Eint_me"+sfx).Data()), "Particles.ID()==614", "goff");

    if (h_Eint_th->Integral() > 0) h_Eint_th->Scale(1.0/h_Eint_th->Integral());
    if (h_Eint_me->Integral() > 0) h_Eint_me->Scale(1.0/h_Eint_me->Integral());

    StylePair(h_Eint_th, h_Eint_me,
              "^{12}C internal energy  #color[616]{(check 4.43 MeV)}",
              "E_{int}^{C} [MeV]", "Norm. events");

    pad_main[3]->cd();
    h_Eint_th->GetYaxis()->SetRangeUser(0, 1.3*TMath::Max(h_Eint_th->GetMaximum(),h_Eint_me->GetMaximum()));
    h_Eint_th->Draw("HIST");
    h_Eint_me->Draw("HIST SAME");
    TLine* l443 = new TLine(4.43, 0, 4.43, h_Eint_th->GetMaximum()*1.2);
    l443->SetLineColor(kGreen+2); l443->SetLineWidth(2); l443->SetLineStyle(2);
    l443->Draw();
    TLatex* lat3i = new TLatex();
    lat3i->SetNDC(); lat3i->SetTextSize(0.042); lat3i->SetTextColor(kGreen+2);
    lat3i->DrawLatex(0.55, 0.87, "4.43 MeV");

    pad_rat[3]->cd();
    DrawRatioPanel(h_Eint_me, h_Eint_th, "E_{int}^{C} [MeV]");

    // ================================================================
    // BOT COL 0: Carbon theta_lab
    // ================================================================
    TH1D* h_thC_th = new TH1D(("ih_thC_th"+sfx).Data(), "", 100, 70, 90);
    TH1D* h_thC_me = new TH1D(("ih_thC_me"+sfx).Data(), "", 100, 70, 90);
    // histogram lives in gDirectory so tree->Draw(">>name") can find it

    tree_th->Draw(Form("Particles.Theta()*TMath::RadToDeg()>>%s", ("ih_thC_th"+sfx).Data()), "Particles.ID()==614", "goff");
    tree_me->Draw(Form("Particles.Theta()*TMath::RadToDeg()>>%s", ("ih_thC_me"+sfx).Data()), "Particles.ID()==614", "goff");

    if (h_thC_th->Integral() > 0) h_thC_th->Scale(1.0/h_thC_th->Integral());
    if (h_thC_me->Integral() > 0) h_thC_me->Scale(1.0/h_thC_me->Integral());

    StylePair(h_thC_th, h_thC_me,
              "^{12}C recoil #theta_{lab}  [inel]",
              "#theta_{lab}^{C} [deg]", "Norm. events");

    pad_bot[0]->cd();
    h_thC_th->GetYaxis()->SetRangeUser(0, 1.3*TMath::Max(h_thC_th->GetMaximum(),h_thC_me->GetMaximum()));
    h_thC_th->Draw("HIST");
    h_thC_me->Draw("HIST SAME");
    TLatex* latCi = new TLatex();
    latCi->SetNDC(); latCi->SetTextSize(0.07); latCi->SetTextColor(kGray+2);
    latCi->DrawLatex(0.18, 0.87, Form("KS p = %.3f", h_thC_th->KolmogorovTest(h_thC_me)));

    // ================================================================
    // BOT COL 1: Asymmetry bar chart
    // ================================================================
    pad_bot[1]->cd();
    gPad->SetGridy();

    // Compute CM angle midpoint for the expected asymmetry reference line
    Double_t theta_cm_lo_i, theta_cm_hi_i;
    {
        std::streambuf* oldCout = std::cout.rdbuf();
        std::ostringstream devnull;
        std::cout.rdbuf(devnull.rdbuf());
        ComputeCMAngleRange(energy,
                            DETECTOR_THETA_CENTER_RAD, DETECTOR_THETA_WINDOW,
                            theta_cm_lo_i, theta_cm_hi_i);
        std::cout.rdbuf(oldCout);
    }
    Double_t theta_cm_ctr_i = 0.5 * (theta_cm_lo_i + theta_cm_hi_i);
    Double_t t_ctr_i  = ConvertThetaCMtoT(theta_cm_ctr_i, energy);
    Double_t AN_meyer_i = MeyerAP_Inelastic(t_ctr_i, energy);
    Double_t A_exp_i  = polarization * AN_meyer_i;

    Double_t x_pos[2]  = {1.0, 2.0};
    Double_t asym_vals[2] = {asym_th, asym_me};
    Double_t asym_errs[2] = {asym_th_err, asym_me_err};
    Double_t x_errs[2] = {0., 0.};

    TGraphErrors* g_asym_i = new TGraphErrors(2, x_pos, asym_vals, x_errs, asym_errs);
    g_asym_i->SetTitle("Azimuthal Asymmetry  [inel];Sampling mode;A = (R-L)/(R+L)");
    g_asym_i->SetMarkerStyle(20);
    g_asym_i->SetMarkerSize(2.0);
    g_asym_i->SetMarkerColor(kBlack);
    g_asym_i->SetLineColor(kBlack);
    g_asym_i->SetLineWidth(2);
    g_asym_i->GetXaxis()->SetLimits(0, 3);
    g_asym_i->GetXaxis()->SetNdivisions(3);
    Double_t ylo_i = TMath::Min(asym_th, asym_me) - 5.*TMath::Max(asym_th_err, asym_me_err);
    Double_t yhi_i = TMath::Max(asym_th, asym_me) + 5.*TMath::Max(asym_th_err, asym_me_err);
    g_asym_i->GetYaxis()->SetRangeUser(ylo_i, yhi_i);
    g_asym_i->Draw("AP");

    TLatex* lAi = new TLatex();
    lAi->SetTextSize(0.07); lAi->SetTextAlign(22);
    lAi->DrawLatex(1.0, asym_th - 3.*asym_th_err, "#theta_{CM}");
    lAi->DrawLatex(2.0, asym_me - 3.*asym_me_err, "Meyer");

    TLine* lexp_i = new TLine(0, A_exp_i, 3, A_exp_i);
    lexp_i->SetLineColor(kGreen+2); lexp_i->SetLineWidth(2); lexp_i->SetLineStyle(2);
    lexp_i->Draw();
    TLatex* lAexpi = new TLatex();
    lAexpi->SetNDC(); lAexpi->SetTextSize(0.065); lAexpi->SetTextColor(kGreen+2);
    lAexpi->DrawLatex(0.18, 0.87, Form("Expected: P#timesA_{N} = %.3f", A_exp_i));

    // ================================================================
    // BOT COL 2: 2D theta-phi for THETA mode
    // ================================================================
    pad_bot[2]->cd();
    gPad->SetRightMargin(0.15);
    TH2D* h2_th_i = new TH2D(("ih2_th"+sfx).Data(),
        "#theta_{CM}: #phi vs #theta  [inel];#phi_{lab} [deg];#theta_{lab} [deg]",
        36, -180, 180, 40, th_lo, th_hi);
    tree_th->Draw(
        Form("Particles.Theta()*TMath::RadToDeg():Particles.Phi()*TMath::RadToDeg()>>%s", ("ih2_th"+sfx).Data()),
        "Particles.ID()==14", "goff");
    h2_th_i->Draw("COLZ");

    // ================================================================
    // BOT COL 3: 2D theta-phi for MEYER mode
    // ================================================================
    pad_bot[3]->cd();
    gPad->SetRightMargin(0.15);
    TH2D* h2_me_i = new TH2D(("ih2_me"+sfx).Data(),
        "Meyer: #phi vs #theta  [inel];#phi_{lab} [deg];#theta_{lab} [deg]",
        36, -180, 180, 40, th_lo, th_hi);
    tree_me->Draw(
        Form("Particles.Theta()*TMath::RadToDeg():Particles.Phi()*TMath::RadToDeg()>>%s", ("ih2_me"+sfx).Data()),
        "Particles.ID()==14", "goff");
    h2_me_i->Draw("COLZ");

    // ================================================================
    // Console summary
    // ================================================================
    cout << "\n══════════════════════════════════════════════════════════" << endl;
    cout << Form("  INELASTIC SUMMARY  |  E = %.0f MeV  |  P = %.0f%%  |  Spin %s",
                 energy, polarization*100, spin_state>0?"UP":"DOWN") << endl;
    cout << "══════════════════════════════════════════════════════════" << endl;
    cout << Form("  %-20s  %12lld events", "#theta_CM sampling:", n_theta) << endl;
    cout << Form("  %-20s  %12lld events", "Meyer sampling:",     n_meyer) << endl;
    cout << endl;
    cout << "  KS compatibility tests:" << endl;
    cout << Form("    theta_lab : p = %.4f %s", ks_theta, ks_theta>0.05?"✓":"✗") << endl;
    cout << Form("    phi_lab   : p = %.4f %s", ks_phi,   ks_phi  >0.05?"✓":"✗") << endl;
    cout << Form("    |p|       : p = %.4f %s", ks_p,     ks_p    >0.05?"✓":"✗") << endl;
    cout << endl;
    cout << "  Azimuthal asymmetry A = (R-L)/(R+L):" << endl;
    cout << Form("    #theta_CM : %.4f +/- %.4f  (R=%.0f, L=%.0f)",
                 asym_th, asym_th_err, nr_th, nl_th) << endl;
    cout << Form("    Meyer     : %.4f +/- %.4f  (R=%.0f, L=%.0f)",
                 asym_me, asym_me_err, nr_me, nl_me) << endl;
    Double_t pull_inel = (asym_th_err>0||asym_me_err>0)
        ? (asym_th-asym_me)/TMath::Sqrt(asym_th_err*asym_th_err+asym_me_err*asym_me_err)
        : 0.;
    cout << Form("    Pull      : %.2f sigma", pull_inel) << endl;
    cout << "══════════════════════════════════════════════════════════\n" << endl;

    c->Update();

    TString outname = Form("Compare_Inelastic_%.0fMeV_P%.0f_%s.pdf",
                           energy, polarization*100,
                           spin_state>0?"SpinUp":"SpinDown");
    c->SaveAs(outname.Data());
    cout << "Saved: " << outname << endl;
}

// ====================================================================
// GenerateAllForComparison
// ====================================================================
// Generates all four files needed by the two comparison functions.
//
// The EventGenerator now appends the mode suffix (_THETA / _MEYER)
// directly to the filename, so no post-generation rename is needed.
//
// Output filenames (example at 180 MeV, P=80%, SpinUp, 16.2 deg):
//   pC_Elas_180MeV_MT_P80_SpinUp_16p2_THETA.root
//   pC_Elas_180MeV_MT_P80_SpinUp_16p2_MEYER.root
//   pC_Inel443_180MeV_MT_P80_SpinUp_16p2_THETA.root
//   pC_Inel443_180MeV_MT_P80_SpinUp_16p2_MEYER.root
// ====================================================================
void GenerateAllForComparison(Double_t energy       = 180.0,
                              Int_t    n_events      = 500000,
                              Double_t polarization  = 0.80,
                              Int_t    spin_state    = +1,
                              Int_t    num_threads   = 0)
{
    cout << "\n╔══════════════════════════════════════════════════════════╗" << endl;
    cout << "║  Generating all four comparison files                   ║" << endl;
    cout << Form("║  E = %.0f MeV   N = %d   P = %.0f%%   Spin: %s",
                 energy, n_events, polarization*100,
                 spin_state>0?"UP  ":"DOWN") << "             ║" << endl;
    cout << "╚══════════════════════════════════════════════════════════╝\n" << endl;

    Int_t    pol_int   = (int)(polarization * 100);
    TString  spinLabel = (spin_state > 0) ? "SpinUp" : "SpinDown";
    Int_t    th_deg    = (int)DETECTOR_THETA_CENTER;
    Int_t    th_dec    = (int)(DETECTOR_THETA_CENTER * 10) % 10;

    // Base pattern (without mode suffix) — mode suffix added by generator
    TString base = Form("_%.0fMeV_MT_P%d_%s_%dp%d",
                        energy, pol_int, spinLabel.Data(), th_deg, th_dec);

    // 1. Elastic — THETA_CM  →  output: pC_Elas<base>_THETA.root/.txt
    cout << "\n[1/4] Elastic THETA_CM sampling..." << endl;
    SingleRunMultithreadPolarized(energy, n_events, polarization, spin_state,
                                  num_threads, THETA_CM_SAMPLING);

    // 2. Elastic — MEYER  →  output: pC_Elas<base>_MEYER.root/.txt
    cout << "\n[2/4] Elastic Meyer sampling..." << endl;
    SingleRunMultithreadPolarized(energy, n_events, polarization, spin_state,
                                  num_threads, T_SAMPLING_MEYER);

    // 3. Inelastic — THETA_CM  →  output: pC_Inel443<base>_THETA.root/.txt
    cout << "\n[3/4] Inelastic THETA_CM sampling..." << endl;
    SingleRunMultithreadInelasticPolarized(energy, n_events, polarization, spin_state,
                                           num_threads, THETA_CM_SAMPLING);

    // 4. Inelastic — MEYER  →  output: pC_Inel443<base>_MEYER.root/.txt
    cout << "\n[4/4] Inelastic Meyer sampling..." << endl;
    SingleRunMultithreadInelasticPolarized(energy, n_events, polarization, spin_state,
                                           num_threads, T_SAMPLING_MEYER);

    // Print ready-to-use comparison commands
    cout << "\n══════════════════════════════════════════════════════════" << endl;
    cout << "All files generated. Run comparisons with:" << endl;
    cout << endl;
    cout << Form("  CompareElastic(\n"
                 "    \"pC_Elas%s_THETA.root\",\n"
                 "    \"pC_Elas%s_MEYER.root\",\n"
                 "    %.0f, %.2f, %d)",
                 base.Data(), base.Data(), energy, polarization, spin_state) << endl;
    cout << endl;
    cout << Form("  CompareInelastic(\n"
                 "    \"pC_Inel443%s_THETA.root\",\n"
                 "    \"pC_Inel443%s_MEYER.root\",\n"
                 "    %.0f, %.2f, %d)",
                 base.Data(), base.Data(), energy, polarization, spin_state) << endl;
    cout << "══════════════════════════════════════════════════════════\n" << endl;
}

// ====================================================================
// CompareAll — convenience wrapper: generate + compare in one call
// ====================================================================
void CompareAll(Double_t energy       = 180.0,
                Int_t    n_events     = 500000,
                Double_t polarization = 0.80,
                Int_t    spin_state   = +1,
                Int_t    num_threads  = 0)
{
    GenerateAllForComparison(energy, n_events, polarization, spin_state, num_threads);

    Int_t    pol_int   = (int)(polarization * 100);
    TString  spinLabel = (spin_state > 0) ? "SpinUp" : "SpinDown";
    Int_t    th_deg    = (int)DETECTOR_THETA_CENTER;
    Int_t    th_dec    = (int)(DETECTOR_THETA_CENTER * 10) % 10;
    TString  base      = Form("_%.0fMeV_MT_P%d_%s_%dp%d",
                               energy, pol_int, spinLabel.Data(), th_deg, th_dec);

    // Filenames now written directly by generator with mode suffix
    TString f_elas_theta = Form("pC_Elas%s_THETA.root",    base.Data());
    TString f_elas_meyer = Form("pC_Elas%s_MEYER.root",    base.Data());
    TString f_inel_theta = Form("pC_Inel443%s_THETA.root", base.Data());
    TString f_inel_meyer = Form("pC_Inel443%s_MEYER.root", base.Data());

    CompareElastic  (f_elas_theta.Data(), f_elas_meyer.Data(), energy, polarization, spin_state);
    CompareInelastic(f_inel_theta.Data(), f_inel_meyer.Data(), energy, polarization, spin_state);
}