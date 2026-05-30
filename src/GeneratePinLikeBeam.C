// ====================================================================
// GeneratePinLikeBeam.C
// ====================================================================
// Generate pin-like beam events at a fixed detector angle.
//
// All N_total events are generated at exactly one (theta_lab, energy).
// The right/left split is determined by the analyzing power at that
// angle and energy, weighted by beam polarization and spin state:
//
//   epsilon  = polarization * A_N(theta_lab, energy) * spin_state
//   N_right  = N_total * (1 + epsilon) / 2   [phi = 0]
//   N_left   = N_total * (1 - epsilon) / 2   [phi = pi]
//
// A_N is taken from the Meyer energy-interpolated splines so the
// correct analyzing power is used at any beam energy in [160,200] MeV.
//
// Usage (single angle):
//   RunPinLikeBeam(200.0, 16.2, 100000, 0.80, +1)
//   RunPinLikeBeam(200.0, 16.2, 100000, 0.80, -1, false)  // inelastic
//
// Usage (grid scan):
//   vector<Double_t> scan = {14.0, 15.0, 16.0, 16.2, 17.0, 18.0};
//   RunPinLikeBeamScan(200.0, scan, 100000, 0.80)
//
// Usage (range scan):
//   RunPinLikeBeamScanRange(200.0, 14.0, 18.0, 0.5, 100000, 0.80)
//
// Verification (requires both spin files to exist):
//   CheckPinLikeBeam("pC_Elas_Pin_200MeV_P80_SpinUp_16p2.root",
//                     200.0, 16.2, 0.80, +1)
//   CheckPinLikeBeamAN("pC_Elas_Pin_200MeV_P80_SpinUp_16p2.root",
//                       "pC_Elas_Pin_200MeV_P80_SpinDown_16p2.root",
//                       200.0, 16.2, 0.80)
//
// Output file naming:
//   pC_Elas_Pin_<E>MeV_P<pol>_<Spin>_<theta>p<dec>.root
//   pC_Inel443_Pin_<E>MeV_P<pol>_<Spin>_<theta>p<dec>.root
//
// Requirements:
//   All PEG modules must be loaded (rootlogon.C)
// ====================================================================

#include <vector>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TCanvas.h"

#include "PParticle.h"

#include "PhysicalConstants.h"
#include "DetectorConfig.h"
#include "Kinematics.h"
#include "ElasticScattering.h"
#include "InelasticScattering.h"
#include "MeyerScattering.h"
#include "AnalyzingPowerUtils.h"

using namespace std;

// ====================================================================
// Internal helpers
// ====================================================================

// Build the output filename
TString PinBeamFilename(Double_t energy,
                        Double_t theta_lab_deg,
                        Double_t polarization,
                        Int_t    spin_state,
                        Bool_t   elastic)
{
    int     pol_int = (int)(polarization * 100 + 0.5);
    TString spin    = (spin_state > 0) ? "SpinUp" : "SpinDown";
    TString channel = elastic ? "pC_Elas" : "pC_Inel443";

    int th_deg = (int)theta_lab_deg;
    int th_dec = (int)(theta_lab_deg * 10 + 0.5) % 10;

    return Form("%s_Pin_%3.0fMeV_P%d_%s_%dp%d.root",
                channel.Data(), energy, pol_int,
                spin.Data(), th_deg, th_dec);
}

// ====================================================================
// ConvertThetaCMtoLab_Inelastic
// ====================================================================
// Same as ConvertThetaCMtoLab but uses the inelastic CM momentum
// (reduced by the 4.43 MeV excitation energy).
// This is needed so the bisection inversion lands at the correct
// theta_lab when we later boost with pcm_inelastic.
// ====================================================================
Double_t ConvertThetaCMtoLab_Inelastic(Double_t energy,
                                        Double_t theta_cm_deg)
{
    Double_t ekinGeV = energy / 1000.;
    Double_t mom     = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);

    // Use inelastic CM momentum — key difference from elastic version
    Double_t pcm = CalculateCMMomentumInelastic(energy, E_EXCITATION);

    Double_t theta_cm_rad = theta_cm_deg * TMath::DegToRad();
    TLorentzVector p_cm;
    p_cm.SetXYZM(pcm * TMath::Sin(theta_cm_rad),
                 0.,
                 pcm * TMath::Cos(theta_cm_rad),
                 mp);
    p_cm.Boost(iState.BoostVector());

    return p_cm.Theta() * TMath::RadToDeg();
}

// ====================================================================
// ThetaLabToCM
// ====================================================================
// Invert theta_lab -> theta_CM via bisection.
// Uses the channel-appropriate forward conversion so the result is
// consistent with the pcm used during event generation.
// ====================================================================
Double_t ThetaLabToCM(Double_t energy,
                       Double_t theta_lab_deg,
                       Bool_t   elastic = true)
{
    Double_t lo = 0., hi = 50., mid = 0.;
    for (int iter = 0; iter < 60; iter++) {
        mid = 0.5*(lo + hi);
        Double_t th_lab = elastic
            ? ConvertThetaCMtoLab(energy, mid)
            : ConvertThetaCMtoLab_Inelastic(energy, mid);
        if (th_lab < theta_lab_deg) lo = mid;
        else                        hi = mid;
    }
    return 0.5*(lo + hi);
}

// ====================================================================
// Core generator — one fixed angle, elastic or inelastic
// ====================================================================
void GeneratePinLikeBeamCore(Double_t energy,
                              Double_t theta_lab_deg,
                              Int_t    n_total,
                              Double_t polarization,
                              Int_t    spin_state,
                              Bool_t   elastic = true)
{
    // ================================================================
    // 1. Kinematics setup
    // ================================================================
    Double_t ekinGeV = energy / 1000.;
    Double_t mom     = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);

    Double_t s   = iState.Mag2();
    Double_t pcm = elastic
        ? TMath::Sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s))
        : CalculateCMMomentumInelastic(energy, E_EXCITATION);

    // Invert theta_lab -> theta_CM using the channel-correct conversion
    Double_t theta_cm_deg = ThetaLabToCM(energy, theta_lab_deg, elastic);
    Double_t theta_cm_rad = theta_cm_deg * TMath::DegToRad();
    Double_t sin_th       = TMath::Sin(theta_cm_rad);
    Double_t cos_th       = TMath::Cos(theta_cm_rad);

    // Verify the round-trip before proceeding
    Double_t theta_lab_check = elastic
        ? ConvertThetaCMtoLab(energy, theta_cm_deg)
        : ConvertThetaCMtoLab_Inelastic(energy, theta_cm_deg);

    if (TMath::Abs(theta_lab_check - theta_lab_deg) > 1e-4) {
        cerr << "WARNING: theta_lab round-trip error = "
             << theta_lab_check - theta_lab_deg << " deg" << endl;
    }

    // ================================================================
    // 2. Analyzing power at this (theta_lab, energy) via Meyer
    // ================================================================
    Double_t t_GeV2 = ConvertThetaCMtoT(theta_cm_deg, energy);
    Double_t AN     = elastic ? MeyerAP_Elastic  (t_GeV2, energy)
                              : MeyerAP_Inelastic(t_GeV2, energy);

    // ================================================================
    // 3. R/L split from asymmetry
    // ================================================================
    Double_t epsilon = polarization * AN * spin_state;
    epsilon = TMath::Max(-1. + 1e-9, TMath::Min(1. - 1e-9, epsilon));

    Int_t n_right = (Int_t)(n_total * (1. + epsilon) / 2. + 0.5);
    Int_t n_left  = n_total - n_right;

    // ================================================================
    // 4. Print summary
    // ================================================================
    cout << "\n========================================" << endl;
    cout << "  PIN-LIKE BEAM GENERATOR" << endl;
    cout << "========================================" << endl;
    cout << "  Channel:      " << (elastic ? "Elastic" : "Inelastic 4.43 MeV") << endl;
    cout << "  Energy:       " << energy        << " MeV"   << endl;
    cout << "  theta_lab:    " << theta_lab_deg << " deg"   << endl;
    cout << "  theta_CM:     " << theta_cm_deg  << " deg"   << endl;
    cout << "  theta_lab check: " << theta_lab_check << " deg" << endl;
    cout << "  t:            " << t_GeV2        << " GeV2"  << endl;
    cout << "  A_N (Meyer):  " << AN                        << endl;
    cout << "  Polarization: " << polarization * 100 << "%" << endl;
    cout << "  Spin:         " << (spin_state > 0 ? "+1 (UP)" : "-1 (DOWN)") << endl;
    cout << "  epsilon:      " << epsilon                   << endl;
    cout << "  N_total:      " << n_total                   << endl;
    cout << "  N_right (phi=0):   " << n_right              << endl;
    cout << "  N_left  (phi=pi):  " << n_left               << endl;
    cout << "========================================" << endl;

    // ================================================================
    // 5. Output file
    // ================================================================
    TString fname = PinBeamFilename(energy, theta_lab_deg,
                                    polarization, spin_state, elastic);
    cout << "  Output: " << fname << endl;
    cout << "========================================\n" << endl;

    TFile* fout = new TFile(fname.Data(), "RECREATE");
    TTree* tree = new TTree("data", elastic
        ? "pC->pC pin-like beam (elastic)"
        : "pC->pC* pin-like beam (inelastic 4.43 MeV)");

    TClonesArray* part = new TClonesArray("PParticle");
    tree->Branch("Particles", &part, 32000, 0);

    Int_t   npart  = 2;
    Float_t impact = 0.;
    Float_t phi0   = 0.;
    tree->Branch("Npart",  &npart);
    tree->Branch("Impact", &impact);
    tree->Branch("Phi",    &phi0);

    // ================================================================
    // 6. Generate events
    // ================================================================
    auto WriteEvent = [&](Double_t phi_cm_rad) {
        TLorentzVector p;
        p.SetXYZM(pcm * sin_th * TMath::Cos(phi_cm_rad),
                  pcm * sin_th * TMath::Sin(phi_cm_rad),
                  pcm * cos_th,
                  mp);
        p.Boost(iState.BoostVector());
        TLorentzVector C = iState - p;

        part->Clear();
        new ((*part)[0]) PParticle(14, p.X(), p.Y(), p.Z());
        new ((*part)[1]) PParticle(614, C.X(), C.Y(), C.Z(),
                                   elastic ? mC : mC + E_EXCITATION, 1);
        tree->Fill();
    };

    for (Int_t i = 0; i < n_right; i++) WriteEvent(0.);
    for (Int_t i = 0; i < n_left;  i++) WriteEvent(TMath::Pi());

    // ================================================================
    // 7. Save and close
    // ================================================================
    tree->Write();
    fout->Close();
    delete fout;

    cout << "Done. Written " << n_total << " events to " << fname << endl;
    cout << "  Right (phi=0):  " << n_right << " events" << endl;
    cout << "  Left  (phi=pi): " << n_left  << " events\n" << endl;
}


// ====================================================================
// CheckPinLikeBeam
// ====================================================================
// Level 1: count N_right / N_left from a single file and compare the
// measured asymmetry to the expected epsilon = P * A_N * spin_state.
//
// Level 2: confirm every proton is at exactly theta_lab_deg.
//
// Usage:
//   CheckPinLikeBeam("pC_Elas_Pin_200MeV_P80_SpinUp_16p2.root",
//                     200.0, 16.2, 0.80, +1)
//   CheckPinLikeBeam("pC_Inel443_Pin_200MeV_P80_SpinUp_16p2.root",
//                     200.0, 16.2, 0.80, +1, false)
// ====================================================================
void CheckPinLikeBeam(const char* filename,
                      Double_t    energy,
                      Double_t    theta_lab_deg,
                      Double_t    polarization,
                      Int_t       spin_state,
                      Bool_t      elastic = true)
{
    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        cerr << "ERROR: cannot open " << filename << endl;
        return;
    }
    TTree* tree = (TTree*)f->Get("data");
    if (!tree) {
        cerr << "ERROR: no 'data' tree in " << filename << endl;
        f->Close();
        return;
    }

    TClonesArray* part = new TClonesArray("PParticle");
    tree->SetBranchAddress("Particles", &part);

    Int_t    nentries = tree->GetEntries();
    Int_t    N_right  = 0, N_left = 0;
    Double_t phi_tol  = 0.01;

    TH1D* h_theta = new TH1D("h_theta_check",
                              Form("theta_lab (expect spike at %.2f deg)"
                                   ";#theta_{lab} [deg];Events",
                                   theta_lab_deg),
                              200, theta_lab_deg - 1.0,
                                   theta_lab_deg + 1.0);
    h_theta->SetDirectory(0);

    for (Int_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        PParticle* proton = (PParticle*)part->At(0);
        Double_t phi      = proton->Phi();
        Double_t theta    = proton->Theta() * TMath::RadToDeg();

        if (TMath::Abs(phi) < phi_tol)
            N_right++;
        else if (TMath::Abs(TMath::Abs(phi) - TMath::Pi()) < phi_tol)
            N_left++;

        h_theta->Fill(theta);
    }

    // ── Level 1: asymmetry check ──────────────────────────────────────
    Double_t asym_measured = (Double_t)(N_right - N_left)
                           / (Double_t)(N_right + N_left);

    Double_t theta_cm = ThetaLabToCM(energy, theta_lab_deg, elastic);
    Double_t t_GeV2   = ConvertThetaCMtoT(theta_cm, energy);
    Double_t AN       = elastic ? MeyerAP_Elastic  (t_GeV2, energy)
                                : MeyerAP_Inelastic(t_GeV2, energy);
    Double_t asym_expected = polarization * AN * spin_state;

    // ── Level 2: theta distribution ───────────────────────────────────
    Double_t theta_mean = h_theta->GetMean();
    Double_t theta_rms  = h_theta->GetRMS();

    // ── Print results ──────────────────────────────────────────────────
    cout << "\n========================================" << endl;
    cout << "  PIN-LIKE BEAM VERIFICATION" << endl;
    cout << "========================================" << endl;
    cout << "  File:     " << filename  << endl;
    cout << "  Channel:  " << (elastic ? "Elastic" : "Inelastic") << endl;
    cout << "  Energy:   " << energy << " MeV" << endl;
    cout << "  Spin:     " << (spin_state > 0 ? "+1 (UP)" : "-1 (DOWN)") << endl;
    cout << "  N_total:  " << nentries << endl;
    cout << "  N_right:  " << N_right  << endl;
    cout << "  N_left:   " << N_left   << endl;
    cout << "\n  --- Level 1: Asymmetry ---" << endl;
    cout << "  A_N (Meyer):       " << AN             << endl;
    cout << "  epsilon expected:  " << asym_expected  << endl;
    cout << "  epsilon measured:  " << asym_measured  << endl;
    cout << "  Difference:        " << asym_measured - asym_expected << endl;
    Bool_t pass1 = (TMath::Abs(asym_measured - asym_expected) < 1e-4);
    cout << "  Result:  " << (pass1 ? "PASS" : "FAIL") << endl;
    cout << "\n  --- Level 2: Theta distribution ---" << endl;
    cout << "  theta_lab requested: " << theta_lab_deg << " deg" << endl;
    cout << "  theta_lab mean:      " << theta_mean    << " deg" << endl;
    cout << "  theta_lab RMS:       " << theta_rms     << " deg" << endl;
    Bool_t pass2 = (TMath::Abs(theta_mean - theta_lab_deg) < 1e-3
                    && theta_rms < 1e-3);
    cout << "  Result:  " << (pass2 ? "PASS" : "FAIL") << endl;
    cout << "========================================\n" << endl;

    // Quick plot of theta distribution
    TString cname = Form("c_pin_check_%s_%s",
                        elastic ? "elas" : "inel",
                        spin_state > 0 ? "up" : "down");
    TCanvas* c = new TCanvas(cname.Data(),
                            Form("Pin-like beam theta check [%s, spin %s]",
                                elastic ? "elastic" : "inelastic",
                                spin_state > 0 ? "UP" : "DOWN"),
                            600, 400);
    h_theta->SetLineColor(kBlue);
    h_theta->SetLineWidth(2);
    h_theta->Draw();
    c->Update();

    f->Close();
    delete part;
}


// ====================================================================
// CheckPinLikeBeamAN
// ====================================================================
// Full pipeline closure test: reconstruct A_N from both spin files
// using CalculateAsymmetry / AsymmetryToAnalyzingPower and compare
// to the Meyer value used during generation.
//
// Usage:
//   CheckPinLikeBeamAN("pC_Elas_Pin_200MeV_P80_SpinUp_16p2.root",
//                       "pC_Elas_Pin_200MeV_P80_SpinDown_16p2.root",
//                       200.0, 16.2, 0.80)
//   CheckPinLikeBeamAN("pC_Inel443_Pin_200MeV_P80_SpinUp_16p2.root",
//                       "pC_Inel443_Pin_200MeV_P80_SpinDown_16p2.root",
//                       200.0, 16.2, 0.80, false)
// ====================================================================
void CheckPinLikeBeamAN(const char* file_up,
                         const char* file_down,
                         Double_t    energy,
                         Double_t    theta_lab_deg,
                         Double_t    polarization,
                         Bool_t      elastic = true)
{
    Double_t phi_window = 0.01;

    Int_t N_R_up,   N_L_up;
    Int_t N_R_down, N_L_down;

    CountEventsInFileByPhi(TString(file_up),   phi_window, true,
                            N_R_up,   N_L_up);
    CountEventsInFileByPhi(TString(file_down), phi_window, false,
                            N_R_down, N_L_down);

    ANResult asym = CalculateAsymmetry(theta_lab_deg, theta_lab_deg,
                                        N_R_up,   N_R_down,
                                        N_L_up,   N_L_down);

    Double_t AN_measured = AsymmetryToAnalyzingPower(asym.AN, polarization);
    Double_t AN_err      = AsymmetryErrorToAnalyzingPowerError(asym.AN_error,
                                                                polarization);

    Double_t theta_cm = ThetaLabToCM(energy, theta_lab_deg, elastic);
    Double_t t_GeV2   = ConvertThetaCMtoT(theta_cm, energy);
    Double_t AN_meyer = elastic ? MeyerAP_Elastic  (t_GeV2, energy)
                                : MeyerAP_Inelastic(t_GeV2, energy);

    Double_t pull = (AN_err > 0)
                  ? (AN_measured - AN_meyer) / AN_err
                  : 999.;

    cout << "\n========================================" << endl;
    cout << "  PIN-LIKE BEAM A_N CLOSURE TEST" << endl;
    cout << "========================================" << endl;
    cout << "  Channel:  " << (elastic ? "Elastic" : "Inelastic") << endl;
    cout << "  Energy:   " << energy        << " MeV" << endl;
    cout << "  theta_lab:" << theta_lab_deg << " deg" << endl;
    cout << "\n  Event counts:" << endl;
    cout << "    N_R_up   = " << N_R_up    << endl;
    cout << "    N_L_up   = " << N_L_up    << endl;
    cout << "    N_R_down = " << N_R_down  << endl;
    cout << "    N_L_down = " << N_L_down  << endl;
    cout << "\n  Asymmetry (raw):   " << asym.AN
         << " +/- " << asym.AN_error      << endl;
    cout << "  A_N measured:      " << AN_measured
         << " +/- " << AN_err             << endl;
    cout << "  A_N Meyer (input): " << AN_meyer    << endl;
    cout << "  Difference:        " << AN_measured - AN_meyer << endl;
    cout << "  Pull:              " << pull << " sigma" << endl;
    Bool_t pass = (TMath::Abs(pull) < 2.0);
    cout << "  Result: " << (pass ? "PASS" : "FAIL") << endl;
    cout << "========================================\n" << endl;
}


// ====================================================================
// RunPinLikeBeam - single angle, single spin
// ====================================================================
void RunPinLikeBeam(Double_t energy        = 200.0,
                    Double_t theta_lab_deg = 16.2,
                    Int_t    n_total       = 100000,
                    Double_t polarization  = 0.80,
                    Int_t    spin_state    = +1,
                    Bool_t   elastic       = true)
{
    GeneratePinLikeBeamCore(energy, theta_lab_deg, n_total,
                             polarization, spin_state, elastic);
}

// ====================================================================
// RunPinLikeBeamBothSpins - both spin states at one angle
// ====================================================================
void RunPinLikeBeamBothSpins(Double_t energy        = 200.0,
                              Double_t theta_lab_deg = 16.2,
                              Int_t    n_total       = 100000,
                              Double_t polarization  = 0.80,
                              Bool_t   elastic       = true)
{
    cout << "\n>>> Generating SPIN UP <<<" << endl;
    GeneratePinLikeBeamCore(energy, theta_lab_deg, n_total,
                             polarization, +1, elastic);

    cout << "\n>>> Generating SPIN DOWN <<<" << endl;
    GeneratePinLikeBeamCore(energy, theta_lab_deg, n_total,
                             polarization, -1, elastic);
}

// ====================================================================
// RunPinLikeBeamBothSpinsAndCheck
// ====================================================================
// Generate both spins then immediately run all verification checks.
//
// Usage:
//   RunPinLikeBeamBothSpinsAndCheck(200.0, 16.2, 100000, 0.80)
//   RunPinLikeBeamBothSpinsAndCheck(200.0, 16.2, 100000, 0.80, false)
// ====================================================================
void RunPinLikeBeamBothSpinsAndCheck(Double_t energy        = 200.0,
                                      Double_t theta_lab_deg = 16.2,
                                      Int_t    n_total       = 100000,
                                      Double_t polarization  = 0.80,
                                      Bool_t   elastic       = true)
{
    RunPinLikeBeamBothSpins(energy, theta_lab_deg, n_total,
                             polarization, elastic);

    TString fup   = PinBeamFilename(energy, theta_lab_deg,
                                    polarization, +1, elastic);
    TString fdown = PinBeamFilename(energy, theta_lab_deg,
                                    polarization, -1, elastic);

    cout << "\n>>> Running individual file checks <<<" << endl;
    CheckPinLikeBeam(fup.Data(),   energy, theta_lab_deg,
                     polarization, +1, elastic);
    CheckPinLikeBeam(fdown.Data(), energy, theta_lab_deg,
                     polarization, -1, elastic);

    cout << "\n>>> Running A_N closure test <<<" << endl;
    CheckPinLikeBeamAN(fup.Data(), fdown.Data(),
                        energy, theta_lab_deg, polarization, elastic);
}

// ====================================================================
// RunPinLikeBeamScan - grid scan over a list of angles
// ====================================================================
void RunPinLikeBeamScan(Double_t         energy,
                         vector<Double_t> angles,
                         Int_t            n_total      = 100000,
                         Double_t         polarization = 0.80,
                         Bool_t           elastic      = true)
{
    cout << "\n========================================" << endl;
    cout << "  PIN-LIKE BEAM GRID SCAN" << endl;
    cout << "========================================" << endl;
    cout << "  Energy:       " << energy        << " MeV"  << endl;
    cout << "  Angles:       " << angles.size() << " positions" << endl;
    cout << "  N per config: " << n_total       << " per spin state" << endl;
    cout << "  Channel:      " << (elastic ? "Elastic" : "Inelastic") << endl;
    cout << "  Total files:  " << angles.size() * 2 << endl;
    cout << "========================================\n" << endl;

    TStopwatch timer;
    timer.Start();

    for (size_t i = 0; i < angles.size(); i++) {
        cout << "\n[" << (i+1) << "/" << angles.size() << "] "
             << "theta_lab = " << angles[i] << " deg" << endl;
        RunPinLikeBeamBothSpins(energy, angles[i], n_total,
                                polarization, elastic);
    }

    timer.Stop();
    cout << "\n========================================" << endl;
    cout << "  SCAN COMPLETE" << endl;
    cout << "  Angles:        " << angles.size()      << endl;
    cout << "  Files written: " << angles.size() * 2  << endl;
    cout << "  Total time:    " << timer.RealTime() << " s" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// RunPinLikeBeamScanRange - range version of the scan
// ====================================================================
void RunPinLikeBeamScanRange(Double_t energy        = 200.0,
                              Double_t angle_min     = 14.0,
                              Double_t angle_max     = 18.0,
                              Double_t angle_step    = 0.5,
                              Int_t    n_total       = 100000,
                              Double_t polarization  = 0.80,
                              Bool_t   elastic       = true)
{
    vector<Double_t> angles;
    for (Double_t a = angle_min; a <= angle_max + 1e-9; a += angle_step)
        angles.push_back(a);

    RunPinLikeBeamScan(energy, angles, n_total, polarization, elastic);
}