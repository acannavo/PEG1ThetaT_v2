// EventGenerator.C
// Implementation of multithreaded event generation
//
// ====================================================================
// SAMPLING MODE GUIDE — READ BEFORE RUNNING
// ====================================================================
//
// Three sampling modes are available via the SamplingMode enum:
//
//   THETA_CM_SAMPLING (0) — Legacy optical model
//   ─────────────────────────────────────────────
//   Elastic:   samples θ_CM from fXS*Op optical model formulas.
//              Valid at discrete energies: 150,160,...,240 MeV.
//              A_N from Wissink parametric formula (energy-dependent).
//   Inelastic: samples θ_CM from digitized CSV data at 200 MeV ONLY.
//              A_N from digitized CSV at 200 MeV ONLY.
//   ⚠ WARNING: For inelastic, this mode uses 200 MeV data at ALL
//              energies. DO NOT use for inelastic away from 200 MeV.
//
//   T_SAMPLING (1) — Energy-independent t-splines
//   ─────────────────────────────────────────────
//   Elastic only. Uses XStSpline(t) and APSpline(t) — single splines
//   with no energy dependence. Useful as a reference but not
//   recommended for physics runs.
//
//   T_SAMPLING_MEYER (2) — Meyer energy-interpolated [DEFAULT]
//   ─────────────────────────────────────────────────────────
//   Elastic & Inelastic. Uses Meyer splines interpolated between
//   160 MeV and 200 MeV data for both XS and A_N.
//   Valid range: 160–200 MeV (clamps to boundary outside this range).
//   ✓ Recommended for all energies and both reactions.
//
// SYSTEMATIC UNCERTAINTIES (measured at detector θ_lab ≈ 16.2°):
//   Elastic  XS:  Meyer/OptModel = 1.11–1.26 across 160–200 MeV
//   Elastic  A_N: |ΔA_N| < 0.04 (largest at 170–190 MeV)
//   Inelastic XS: Meyer/CSV ≈ 1.00–1.07 (good agreement)
//   Inelastic A_N: |ΔA_N| up to 0.30 at 160 MeV for THETA_CM mode
//                  (because CSV A_N is fixed at 200 MeV)
//
// DEFAULT: T_SAMPLING_MEYER for both elastic and inelastic.
// ====================================================================

#include "EventGenerator.h"
#include "ThreadUtils.h"
#include "ElasticScattering.h"
#include "InelasticScattering.h"
#include "Kinematics.h"
#include "DetectorConfig.h"
#include "PhysicalConstants.h"

// ROOT headers
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TMath.h"

// Meyer scattering module (energy-dependent XS and AP splines)
#include "MeyerScattering.h"

// Pluto headers
#include "PParticle.h"

// C++ standard headers
#include <ROOT/TThreadExecutor.hxx>
#include <iostream>
#include <chrono>
#include <thread>
#include <ctime>

using namespace std;



// ====================================================================
// OUTPUT CONFIGURATION
// ====================================================================
// Set this to your desired output directory
// Leave empty ("") to use current working directory
// Examples:
//   const TString OUTPUT_BASE_PATH = "";  // Current directory
//   const TString OUTPUT_BASE_PATH = "/home/user/output/";  // Specific path
//   const TString OUTPUT_BASE_PATH = "./results/";  // Relative path
// 
// IMPORTANT: Path must end with "/" if non-empty
// ====================================================================
const TString OUTPUT_BASE_PATH = "";  // Default: current directory

// Helper function to construct full output path
TString MakeOutputPath(const char* filename) {
    if (OUTPUT_BASE_PATH.IsNull() || OUTPUT_BASE_PATH == "") {
        return TString(filename);
    }
    // Ensure path ends with /
    TString path = OUTPUT_BASE_PATH;
    if (!path.EndsWith("/")) {
        path += "/";
    }
    return path + TString(filename);
}



// ====================================================================
// ELASTIC SCATTERING
// ====================================================================



// ====================================================================
// GenerateThreadEventsPolarized - Worker function with polarization
// ====================================================================
ThreadData GenerateThreadEventsPolarized(int thread_id, Double_t ekin, Int_t events_per_thread,
                                          Double_t theta_lab_target, Double_t theta_lab_window,
                                          Double_t phi_lab_window, const char* xsFormula,
                                          Double_t polarization, Int_t spin_state)
{
	ThreadData data;
	data.count = 0;
	
	Double_t ekinGeV = ekin/1000.;
	Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
	TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
	
	Double_t s = iState.Mag2();
	Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));
	
	TF1 *ffps = new TF1(Form("ffps_t%d", thread_id), xsFormula, 
	                     0.*TMath::DegToRad(), 50.*TMath::DegToRad());
	
	Double_t theta_cm_min, theta_cm_max;
	ComputeCMAngleRange(ekin, theta_lab_target, theta_lab_window, 
	                    theta_cm_min, theta_cm_max);
	
	// BUG FIX: XSLog*MeV splines use theta_CM in DEGREES as their x-axis
	// (x-array runs 1,2,...,80 degrees). The TF1 formula "fXS200Op(x)" passes
	// x directly to the spline, so we must evaluate it with theta_CM in degrees.
	// The original code called ffps_narrow->Eval(theta_deg * DegToRad()), which
	// passed radians (~0.31) to a spline expecting degrees — giving XS = 3e6 mb/sr
	// instead of ~40 mb/sr. We evaluate the formula directly in degrees instead.
	//
	// The TF1 is still constructed (needed for GetAnalyzingPower calls elsewhere)
	// but the histogram is now filled by calling the formula function directly.
	TF1 *ffps_narrow = new TF1(Form("ffps_narrow_t%d", thread_id), xsFormula,
	                            theta_cm_min, theta_cm_max);  // axis now in degrees
	
	TH1D* sigmacm_narrow = new TH1D(Form("sigmacm_narrow_t%d", thread_id), 
	                                 Form("narrow range t%d", thread_id), 
	                                 10000, theta_cm_min, theta_cm_max);
	for (int ibin = 1; ibin <= 10000; ibin++) {
		Double_t theta_deg = theta_cm_min + (theta_cm_max - theta_cm_min) * (ibin-0.5) / 10000.;
		// Evaluate with theta_CM in DEGREES — matches XSLog*MeV spline x-axis
		sigmacm_narrow->SetBinContent(ibin, ffps_narrow->Eval(theta_deg));
	}
	
	const Double_t theta_lab_min = theta_lab_target - theta_lab_window;
	const Double_t theta_lab_max = theta_lab_target + theta_lab_window;
	const Double_t phi_lab_min_0   = -phi_lab_window;
	const Double_t phi_lab_max_0   = phi_lab_window;
	const Double_t phi_lab_min_180 = TMath::Pi() - phi_lab_window;
	const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;
	
	TRandom3 rng(thread_id + std::time(nullptr));
	
	int li = 0;
	while ( li < events_per_thread )
	{
		Double_t theta = sigmacm_narrow->GetRandom();
		Double_t stheta = TMath::Sin(TMath::DegToRad() * theta);
		Double_t ctheta = TMath::Cos(TMath::DegToRad() * theta);
		
		Double_t phi;
		if (polarization == 0.0) {
			phi = rng.Rndm() * 2. * TMath::Pi();
		} else {
			Double_t theta_lab_deg = ConvertThetaCMtoLab(ekin, theta);
			Double_t AN = GetElasticAnalyzingPower(ekin, theta_lab_deg);
			Double_t asymmetry = spin_state * polarization * AN;
			
			Double_t phi_test, weight;
			do {
				phi_test = rng.Rndm() * 2. * TMath::Pi();
				weight = 1.0 + asymmetry * TMath::Cos(phi_test);
			} while (rng.Rndm() > weight / (1.0 + TMath::Abs(asymmetry)));
			phi = phi_test;
		}
		
		TLorentzVector p;
		p.SetXYZM(pcm * stheta * cos(phi), pcm * stheta * sin(phi), pcm * ctheta, mp);
		p.Boost(iState.BoostVector());
		TLorentzVector C = iState - p;
		
		Bool_t theta_accepted = (p.Theta() >= theta_lab_min && p.Theta() <= theta_lab_max);
		Bool_t phi_accepted_0   = (p.Phi() >= phi_lab_min_0   && p.Phi() <= phi_lab_max_0);
		Bool_t phi_accepted_180 = (p.Phi() >= phi_lab_min_180 || p.Phi() <= phi_lab_max_180);
		
		if ( theta_accepted && (phi_accepted_0 || phi_accepted_180) )
		{
			li++;
			data.event_ids.push_back(li);
			data.px.push_back(p.X());
			data.py.push_back(p.Y());
			data.pz.push_back(p.Z());
			data.cx.push_back(C.X());
			data.cy.push_back(C.Y());
			data.cz.push_back(C.Z());
		}
	}
	
	data.count = li;
	delete ffps;
	delete ffps_narrow;
	delete sigmacm_narrow;
	
	return data;
}

// ====================================================================
// GenerateThreadEventsPolarized_TSampling - t-based sampling worker
// ====================================================================
ThreadData GenerateThreadEventsPolarized_TSampling(
    int thread_id, 
    Double_t ekin, 
    Int_t events_per_thread,
    Double_t t_min,           // t range in GeV² (negative values)
    Double_t t_max,
    Double_t theta_lab_target, 
    Double_t theta_lab_window,
    Double_t phi_lab_window, 
    Double_t polarization, 
    Int_t spin_state)
{
	ThreadData data;
	data.count = 0;
	
	// Initial state kinematics (same as theta-based)
	Double_t ekinGeV = ekin/1000.;
	Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
	TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
	
	Double_t s = iState.Mag2();
	Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));
	
	// ================================================================
	// KEY DIFFERENCE: Create t-sampling histogram instead of theta
	// ================================================================
	TH1D* h_t_sampling = new TH1D(Form("h_t_sampling_t%d", thread_id), 
	                               "t sampling", 
	                               10000, t_min, t_max);
	
	// Fill histogram with cross section from spline
	// IMPORTANT: Splines use |t| in (MeV/c)², not t in GeV²
	for (int ibin = 1; ibin <= 10000; ibin++) {
		Double_t t_GeV2 = h_t_sampling->GetBinCenter(ibin);
		Double_t t_MeV2 = TMath::Abs(t_GeV2) * 1.0e6;  // |t| in (MeV/c)²
		Double_t xs = XStSpline(t_MeV2);
		h_t_sampling->SetBinContent(ibin, xs);
	}
	
	// Lab acceptance cuts (same as theta-based)
	const Double_t theta_lab_min = theta_lab_target - theta_lab_window;
	const Double_t theta_lab_max = theta_lab_target + theta_lab_window;
	const Double_t phi_lab_min_0   = -phi_lab_window;
	const Double_t phi_lab_max_0   = phi_lab_window;
	const Double_t phi_lab_min_180 = TMath::Pi() - phi_lab_window;
	const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;
	
	// Thread-local random number generator
	TRandom3 rng(thread_id + std::time(nullptr));
	
	int li = 0;
	while (li < events_per_thread)
	{
		// ============================================================
		// STEP 1: Sample t from histogram
		// ============================================================
		Double_t t_GeV2 = h_t_sampling->GetRandom();
		
		// ============================================================
		// STEP 2: Convert t → theta_CM using actual beam energy
		// ============================================================
		Double_t theta = ConvertTtoThetaCM(t_GeV2, ekin);
		Double_t stheta = TMath::Sin(TMath::DegToRad() * theta);
		Double_t ctheta = TMath::Cos(TMath::DegToRad() * theta);
		
		// ============================================================
		// STEP 3: Sample phi with t-based analyzing power
		// ============================================================
		Double_t phi;
		if (polarization == 0.0) {
			// Unpolarized: uniform phi
			phi = rng.Rndm() * 2. * TMath::Pi();
		} else {
			// Polarized: use APSpline(|t|) instead of GetElasticAnalyzingPower
			// CRITICAL: Use absolute value of t and convert to (MeV/c)²
			Double_t t_MeV2 = TMath::Abs(t_GeV2) * 1.0e6;  // |t| in (MeV/c)²
			Double_t AN = APSpline(t_MeV2);
			Double_t asymmetry = spin_state * polarization * AN;
			
			// Acceptance-rejection sampling (same as theta-based)
			Double_t phi_test, weight;
			do {
				phi_test = rng.Rndm() * 2. * TMath::Pi();
				weight = 1.0 + asymmetry * TMath::Cos(phi_test);
			} while (rng.Rndm() > weight / (1.0 + TMath::Abs(asymmetry)));
			phi = phi_test;
		}
		
		// ============================================================
		// STEP 4: Boost to lab frame (identical to theta-based)
		// ============================================================
		TLorentzVector p;
		p.SetXYZM(pcm * stheta * cos(phi), pcm * stheta * sin(phi), pcm * ctheta, mp);
		p.Boost(iState.BoostVector());
		TLorentzVector C = iState - p;
		
		// ============================================================
		// STEP 5: Apply lab acceptance cuts (identical to theta-based)
		// ============================================================
		Bool_t theta_accepted = (p.Theta() >= theta_lab_min && p.Theta() <= theta_lab_max);
		Bool_t phi_accepted_0   = (p.Phi() >= phi_lab_min_0   && p.Phi() <= phi_lab_max_0);
		Bool_t phi_accepted_180 = (p.Phi() >= phi_lab_min_180 || p.Phi() <= phi_lab_max_180);
		
		if (theta_accepted && (phi_accepted_0 || phi_accepted_180))
		{
			li++;
			data.event_ids.push_back(li);
			data.px.push_back(p.X());
			data.py.push_back(p.Y());
			data.pz.push_back(p.Z());
			data.cx.push_back(C.X());
			data.cy.push_back(C.Y());
			data.cz.push_back(C.Z());
		}
	}
	
	data.count = li;
	delete h_t_sampling;
	
	return data;
}

// ====================================================================
// GenerateThreadEventsPolarized_MeyerSampling
// ====================================================================
// ELASTIC Meyer worker — energy-interpolated XS and AP.
// Thread-safe: all state is local; envelope is shared read-only.
// ====================================================================
ThreadData GenerateThreadEventsPolarized_MeyerSampling(
    int thread_id,
    Double_t ekin,
    Int_t events_per_thread,
    TH1D* envelope,
    Double_t theta_lab_target,
    Double_t theta_lab_window,
    Double_t phi_lab_window,
    Double_t polarization,
    Int_t spin_state)
{
    ThreadData data;
    data.count = 0;

    // ================================================================
    // Initial state kinematics (identical to TSampling version)
    // ================================================================
    Double_t ekinGeV = ekin / 1000.;
    Double_t mom     = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);

    Double_t s   = iState.Mag2();
    Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));

    // ================================================================
    // Lab acceptance windows (identical to all other workers)
    // ================================================================
    const Double_t theta_lab_min   = theta_lab_target - theta_lab_window;
    const Double_t theta_lab_max   = theta_lab_target + theta_lab_window;
    const Double_t phi_lab_min_0   = -phi_lab_window;
    const Double_t phi_lab_max_0   =  phi_lab_window;
    const Double_t phi_lab_min_180 =  TMath::Pi() - phi_lab_window;
    const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;

    // ================================================================
    // Thread-local RNG — unique seed per thread.
    // We do NOT call SampleT_Meyer_Elastic() because it owns a static
    // TRandom3 that is shared across all callers and is not thread-safe.
    // The accept/reject logic is reproduced inline here instead.
    // ================================================================
    TRandom3 rng(thread_id + std::time(nullptr));

    int li = 0;
    while (li < events_per_thread)
    {
        // ============================================================
        // STEP 1: Meyer two-level accept/reject on t
        // ============================================================

        // 1a. Sample candidate t from the envelope — thread-safe.
        // We cannot call envelope->GetRandom() because TH1::GetRandom()
        // uses gRandom (the global ROOT RNG) internally and is not
        // thread-safe when multiple threads share the same TH1D object.
        // Instead we reproduce the inverse-CDF sampling manually using
        // the pre-built fIntegral array (read-only after ComputeIntegral)
        // and the thread-local rng.
        const Double_t* integral = envelope->GetIntegral();  // read-only
        Int_t nbins = envelope->GetNbinsX();
        Double_t u = rng.Rndm();  // uniform [0,1)
        // Binary search on the CDF array (size nbins+2: [0]=0, [nbins+1]=1)
        Int_t klow = 0, khig = nbins, khalf;
        while (khig - klow > 1) {
            khalf = (klow + khig) / 2;
            if (u > integral[khalf]) klow = khalf;
            else                     khig = khalf;
        }
        // Interpolate within the bin (linear)
        Double_t dI = integral[klow+1] - integral[klow];
        Double_t t_GeV2;
        if (dI > 0.)
            t_GeV2 = envelope->GetXaxis()->GetBinLowEdge(klow+1)
                   + envelope->GetXaxis()->GetBinWidth(klow+1)
                     * (u - integral[klow]) / dI;
        else
            t_GeV2 = envelope->GetXaxis()->GetBinCenter(klow+1);

        // 1b. Retrieve envelope upper bound at this t
        Int_t    bin         = envelope->FindBin(t_GeV2);
        Double_t xs_envelope = envelope->GetBinContent(bin);

        // Guard: skip empty bins that can appear at histogram edges
        if (xs_envelope <= 0.) continue;

        // 1c. Evaluate energy-interpolated cross section at this t and ekin
        Double_t xs_actual = MeyerXS_Elastic(t_GeV2, ekin);

        // 1d. Accept/reject: P(accept) = xs_actual / xs_envelope
        if (rng.Rndm() > xs_actual / xs_envelope) continue;  // REJECTED

        // t_GeV2 is now an accepted physical momentum transfer.

        // ============================================================
        // STEP 2: Convert t → theta_CM using actual beam energy
        // ============================================================
        Double_t theta  = ConvertTtoThetaCM(t_GeV2, ekin);  // [degrees]
        Double_t stheta = TMath::Sin(TMath::DegToRad() * theta);
        Double_t ctheta = TMath::Cos(TMath::DegToRad() * theta);

        // ============================================================
        // STEP 3: Sample phi using Meyer energy-interpolated A_N
        // ============================================================
        Double_t phi;
        if (polarization == 0.0) {
            // Unpolarized: uniform azimuthal distribution
            phi = rng.Rndm() * 2. * TMath::Pi();
        } else {
            // Polarized: sample from 1 + P·A_N(t, ekin)·cos(phi)
            // MeyerAP_Elastic expects t in GeV² (negative); it converts internally
            Double_t AN        = MeyerAP_Elastic(t_GeV2, ekin);
            Double_t asymmetry = spin_state * polarization * AN;

            Double_t phi_test, weight;
            do {
                phi_test = rng.Rndm() * 2. * TMath::Pi();
                weight   = 1.0 + asymmetry * TMath::Cos(phi_test);
            } while (rng.Rndm() > weight / (1.0 + TMath::Abs(asymmetry)));
            phi = phi_test;
        }

        // ============================================================
        // STEP 4: Build 4-momenta and boost to lab frame
        // ============================================================
        TLorentzVector p;
        p.SetXYZM(pcm * stheta * TMath::Cos(phi),
                  pcm * stheta * TMath::Sin(phi),
                  pcm * ctheta,
                  mp);
        p.Boost(iState.BoostVector());
        TLorentzVector C = iState - p;

        // ============================================================
        // STEP 5: Apply lab acceptance cuts
        // ============================================================
        Bool_t theta_accepted   = (p.Theta() >= theta_lab_min &&
                                   p.Theta() <= theta_lab_max);
        Bool_t phi_accepted_0   = (p.Phi()   >= phi_lab_min_0 &&
                                   p.Phi()   <= phi_lab_max_0);
        Bool_t phi_accepted_180 = (p.Phi()   >= phi_lab_min_180 ||
                                   p.Phi()   <= phi_lab_max_180);

        if (!theta_accepted || !(phi_accepted_0 || phi_accepted_180)) continue;

        // ============================================================
        // STEP 6: Store accepted event
        // ============================================================
        li++;
        data.event_ids.push_back(li);
        data.px.push_back(p.X());
        data.py.push_back(p.Y());
        data.pz.push_back(p.Z());
        data.cx.push_back(C.X());
        data.cy.push_back(C.Y());
        data.cz.push_back(C.Z());
    }

    data.count = li;
    // envelope is owned by the orchestrator — do NOT delete here
    return data;
}

// ====================================================================
// GenerateThreadEventsInelastic_MeyerSampling
// ====================================================================
// INELASTIC Meyer worker — energy-interpolated XS and AP.
// Mirrors the elastic Meyer worker with three key differences:
//   - Uses MeyerXS_Inelastic / MeyerAP_Inelastic
//   - Uses reduced CM momentum for the 4.43 MeV excitation
//   - Carbon recoil stored with excited-state mass (mC + E_EXCITATION)
// ====================================================================
ThreadData GenerateThreadEventsInelastic_MeyerSampling(
    int thread_id,
    Double_t ekin,
    Int_t events_per_thread,
    TH1D* envelope,
    Double_t theta_lab_target,
    Double_t theta_lab_window,
    Double_t phi_lab_window,
    Double_t polarization,
    Int_t spin_state)
{
    ThreadData data;
    data.count = 0;

    // ================================================================
    // Initial state kinematics
    // ================================================================
    Double_t ekinGeV = ekin / 1000.;
    Double_t mom     = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);

    // INELASTIC: reduced CM momentum accounting for 4.43 MeV excitation
    Double_t pcm = CalculateCMMomentumInelastic(ekin, 0.00443);

    // ================================================================
    // Lab acceptance windows
    // ================================================================
    const Double_t theta_lab_min   = theta_lab_target - theta_lab_window;
    const Double_t theta_lab_max   = theta_lab_target + theta_lab_window;
    const Double_t phi_lab_min_0   = -phi_lab_window;
    const Double_t phi_lab_max_0   =  phi_lab_window;
    const Double_t phi_lab_min_180 =  TMath::Pi() - phi_lab_window;
    const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;

    // Thread-local RNG
    TRandom3 rng(thread_id + std::time(nullptr));

    int li = 0;
    while (li < events_per_thread)
    {
        // ============================================================
        // STEP 1: Meyer two-level accept/reject on t (inelastic XS)
        // ============================================================

        // Thread-safe inverse-CDF sampling (see elastic worker comment)
        const Double_t* integral_i = envelope->GetIntegral();
        Int_t nbins_i = envelope->GetNbinsX();
        Double_t u_i = rng.Rndm();
        Int_t klow_i = 0, khig_i = nbins_i, khalf_i;
        while (khig_i - klow_i > 1) {
            khalf_i = (klow_i + khig_i) / 2;
            if (u_i > integral_i[khalf_i]) klow_i = khalf_i;
            else                           khig_i = khalf_i;
        }
        Double_t dI_i = integral_i[klow_i+1] - integral_i[klow_i];
        Double_t t_GeV2;
        if (dI_i > 0.)
            t_GeV2 = envelope->GetXaxis()->GetBinLowEdge(klow_i+1)
                   + envelope->GetXaxis()->GetBinWidth(klow_i+1)
                     * (u_i - integral_i[klow_i]) / dI_i;
        else
            t_GeV2 = envelope->GetXaxis()->GetBinCenter(klow_i+1);

        Int_t    bin         = envelope->FindBin(t_GeV2);
        Double_t xs_envelope = envelope->GetBinContent(bin);

        if (xs_envelope <= 0.) continue;

        Double_t xs_actual = MeyerXS_Inelastic(t_GeV2, ekin);

        if (rng.Rndm() > xs_actual / xs_envelope) continue;  // REJECTED

        // ============================================================
        // STEP 2: Convert t → theta_CM
        // ============================================================
        Double_t theta  = ConvertTtoThetaCM(t_GeV2, ekin);  // [degrees]
        Double_t stheta = TMath::Sin(TMath::DegToRad() * theta);
        Double_t ctheta = TMath::Cos(TMath::DegToRad() * theta);

        // ============================================================
        // STEP 3: Sample phi using inelastic Meyer A_N
        // ============================================================
        Double_t phi;
        if (polarization == 0.0) {
            phi = rng.Rndm() * 2. * TMath::Pi();
        } else {
            Double_t AN        = MeyerAP_Inelastic(t_GeV2, ekin);
            Double_t asymmetry = spin_state * polarization * AN;

            Double_t phi_test, weight;
            do {
                phi_test = rng.Rndm() * 2. * TMath::Pi();
                weight   = 1.0 + asymmetry * TMath::Cos(phi_test);
            } while (rng.Rndm() > weight / (1.0 + TMath::Abs(asymmetry)));
            phi = phi_test;
        }

        // ============================================================
        // STEP 4: Build 4-momenta with inelastic pcm and boost to lab
        // ============================================================
        TLorentzVector p;
        p.SetXYZM(pcm * stheta * TMath::Cos(phi),
                  pcm * stheta * TMath::Sin(phi),
                  pcm * ctheta,
                  mp);
        p.Boost(iState.BoostVector());
        TLorentzVector C = iState - p;

        // ============================================================
        // STEP 5: Apply lab acceptance cuts
        // ============================================================
        Bool_t theta_accepted   = (p.Theta() >= theta_lab_min &&
                                   p.Theta() <= theta_lab_max);
        Bool_t phi_accepted_0   = (p.Phi()   >= phi_lab_min_0 &&
                                   p.Phi()   <= phi_lab_max_0);
        Bool_t phi_accepted_180 = (p.Phi()   >= phi_lab_min_180 ||
                                   p.Phi()   <= phi_lab_max_180);

        if (!theta_accepted || !(phi_accepted_0 || phi_accepted_180)) continue;

        // ============================================================
        // STEP 6: Store accepted event
        // Carbon recoil carries excited-state mass (status flag = 1)
        // — same convention as GenerateThreadEventsInelasticPolarized
        // ============================================================
        li++;
        data.event_ids.push_back(li);
        data.px.push_back(p.X());
        data.py.push_back(p.Y());
        data.pz.push_back(p.Z());
        data.cx.push_back(C.X());
        data.cy.push_back(C.Y());
        data.cz.push_back(C.Z());
    }

    data.count = li;
    return data;
}

// ====================================================================
// SingleRunMultithreadPolarized - Polarized parallel event generation
// ====================================================================
void SingleRunMultithreadPolarized(Double_t energy, Int_t number_total, 
                                    Double_t polarization, Int_t spin_state,
                                    Int_t num_threads,
                                    SamplingMode mode)
{
	std::cout << "\n=== Multithreaded POLARIZED Event Generation ===" << std::endl;
	std::cout << "Energy: " << energy << " MeV" << std::endl;
	std::cout << "Total events: " << number_total << std::endl;
	std::cout << "Polarization: " << polarization * 100 << "%" << std::endl;
	std::cout << "Spin state: " << (spin_state > 0 ? "UP" : "DOWN") << std::endl;
	
	// Display sampling mode
	const char* mode_str;
	if (mode == T_SAMPLING) {
		mode_str = "t-based (energy-independent)";
	} else if (mode == T_SAMPLING_MEYER) {
		mode_str = "t-based Meyer (energy-dependent)";
	} else {
		mode_str = "θ_CM-based (legacy)";
	}
	std::cout << "Sampling mode: " << mode_str << std::endl;

	// ── Energy range validation ────────────────────────────────
	if (mode == T_SAMPLING_MEYER && (energy < 160.0 || energy > 200.0)) {
		std::cout << "WARNING: T_SAMPLING_MEYER is validated for 160-200 MeV." << std::endl;
		std::cout << "         Energy " << energy << " MeV is outside this range." << std::endl;
		std::cout << "         Splines will clamp to the nearest boundary." << std::endl;
	}
	if (mode == THETA_CM_SAMPLING) {
		// Check if this energy has an optical model formula
		const char* formula_check = GetXSFormula(energy);
		if (std::string(formula_check) == "cos(x)") {
			std::cout << "WARNING: THETA_CM_SAMPLING has no optical model data" << std::endl;
			std::cout << "         at " << energy << " MeV. Using cos(x) fallback." << std::endl;
			std::cout << "         Consider using T_SAMPLING_MEYER instead." << std::endl;
		}
	}
	
	auto start = std::chrono::high_resolution_clock::now();
	
	if (num_threads == 0) {
		num_threads = std::thread::hardware_concurrency();
	}
	std::cout << "Number of threads: " << num_threads << std::endl;
	
	Int_t events_per_thread = number_total / num_threads;
	
	const Double_t theta_lab_target = DETECTOR_THETA_CENTER_RAD;
	const Double_t theta_lab_window = DETECTOR_THETA_WINDOW;
	const Double_t phi_lab_window   = DETECTOR_PHI_WINDOW;
	
	ROOT::TThreadExecutor executor(num_threads);
	
	// ================================================================
	// BRANCH: Choose worker function based on sampling mode
	// ================================================================
	std::vector<ThreadData> results;
	
	if (mode == T_SAMPLING_MEYER) {
		// ============================================================
		// MEYER SAMPLING MODE — energy-interpolated XS and A_N
		// ============================================================
		// Build the envelope once here (outside threads).
		// All threads share it as read-only; the orchestrator owns it
		// and deletes it only after executor.Map() returns.
		// ============================================================
		std::cout << "\nCalculating unified t range for Meyer envelope..." << std::endl;

		Double_t t_min, t_max;
		ComputeUnifiedTRange(energy, energy,
		                     theta_lab_target, theta_lab_window,
		                     t_min, t_max);

		std::cout << "  t range: [" << t_min << ", " << t_max << "] GeV²" << std::endl;
		std::cout << "Building Meyer elastic envelope (50k bins)..." << std::endl;

		TH1D* meyer_envelope = CreateMeyerEnvelope_Elastic(t_min, t_max);
		meyer_envelope->SetDirectory(0);
		// THREAD SAFETY: pre-build fIntegral before threads launch.
		// TH1::GetRandom() lazily initialises fIntegral and uses gRandom
		// internally — both are not thread-safe under TThreadExecutor.
		// After ComputeIntegral() the array is read-only: safe to share.
		meyer_envelope->ComputeIntegral();

		std::cout << "  Envelope max: " << meyer_envelope->GetMaximum()
		          << " mb/GeV²" << std::endl;

		auto thread_work = [energy, events_per_thread,
		                    meyer_envelope,
		                    theta_lab_target, theta_lab_window, phi_lab_window,
		                    polarization, spin_state]
		                   (int tid) {
			return GenerateThreadEventsPolarized_MeyerSampling(
				tid, energy, events_per_thread,
				meyer_envelope,
				theta_lab_target, theta_lab_window, phi_lab_window,
				polarization, spin_state);
		};

		std::cout << "Starting multithreaded Meyer elastic generation..." << std::endl;
		results = executor.Map(thread_work, ROOT::TSeq<int>(0, num_threads));

		// All threads have finished — safe to release the envelope
		delete meyer_envelope;

	} else if (mode == T_SAMPLING) {
		// ============================================================
		// ORIGINAL t-BASED SAMPLING MODE (energy-independent splines)
		// ============================================================
		std::cout << "\nCalculating unified t range..." << std::endl;

		Double_t t_min, t_max;
		ComputeUnifiedTRange(energy, energy,
		                     theta_lab_target, theta_lab_window,
		                     t_min, t_max);

		auto thread_work = [energy, events_per_thread,
		                    t_min, t_max,
		                    theta_lab_target, theta_lab_window, phi_lab_window,
		                    polarization, spin_state]
		                   (int tid) {
			return GenerateThreadEventsPolarized_TSampling(
				tid, energy, events_per_thread,
				t_min, t_max,
				theta_lab_target, theta_lab_window, phi_lab_window,
				polarization, spin_state);
		};

		std::cout << "Starting multithreaded t-based generation..." << std::endl;
		results = executor.Map(thread_work, ROOT::TSeq<int>(0, num_threads));

	} else {
		// ============================================================
		// THETA-BASED SAMPLING MODE (original)
		// ============================================================
		const char* xsFormula = GetXSFormula(energy);
		
		// Define theta-based worker lambda
		auto thread_work = [energy, events_per_thread, 
		                    theta_lab_target, theta_lab_window, phi_lab_window,
		                    xsFormula, polarization, spin_state]
		                   (int tid) {
			return GenerateThreadEventsPolarized(
				tid, energy, events_per_thread,
				theta_lab_target, theta_lab_window, phi_lab_window,
				xsFormula, polarization, spin_state);
		};
		
		std::cout << "Starting multithreaded theta-based generation..." << std::endl;
		results = executor.Map(thread_work, ROOT::TSeq<int>(0, num_threads));
	}
	
	// ================================================================
	// COMMON: Write output files (same for both modes)
	// ================================================================
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	
	std::cout << "\n=== Writing to ROOT file ===" << std::endl;
	
	// Prepare output files with standardized names
	int polar_int = (int)(polarization * 100);
	TString spinLabel  = (spin_state > 0) ? "SpinUp" : "SpinDown";
	TString modeLabel  = (mode == T_SAMPLING_MEYER) ? "MEYER"
	                   : (mode == T_SAMPLING)        ? "T"
	                   :                               "THETA";
	TString basename = Form("pC_Elas_%3.0fMeV_MT_P%d_%s_%dp%d_%s",
	                    energy, polar_int, spinLabel.Data(),
	                    (int)DETECTOR_THETA_CENTER, (int)(DETECTOR_THETA_CENTER*10)%10,
	                    modeLabel.Data());
	TString nameo = MakeOutputPath(basename + ".txt");
	TString namer = MakeOutputPath(basename + ".root");

	std::cout << "Output files:" << std::endl;
	std::cout << "  ROOT: " << namer << std::endl;
	std::cout << "  TXT:  " << nameo << std::endl;
	std::cout << "========================================\n" << std::endl;
	    
	// Create ROOT file and tree - STANDARD PLUTO FORMAT
	TFile f(namer.Data(), "RECREATE");
	TTree* tree = new TTree("data", "pC->pC elastic");
	TClonesArray* part = new TClonesArray("PParticle");
	tree->Branch("Particles", &part, 32000, 0);
	
	// Standard PLUTO branches ONLY
	Int_t npart = 2;
	Float_t impact = 0.;
	Float_t phi0 = 0.;
	tree->Branch("Npart", &npart);
	tree->Branch("Impact", &impact);
	tree->Branch("Phi", &phi0);
	
	// Open ASCII output file
	FILE* fout = fopen(nameo.Data(), "w");
	
	// Write ASCII header
	Double_t ekinGeV = energy/1000.;
	Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
	int codr = 301;  // 301 for elastic
	int nbp = 2;
	double wei = 1.000E+00;
	int fform = 10000000 + codr;
	fprintf(fout, " %10d%10.2E%10.3f%10.3f%10.3f%10.3f%7d\n", 
	        fform, 7.36E-8, mom, 0.f, 0.f, 0.f, number_total);
	fprintf(fout, " REAC,CROSS(mb),B. MOM,  A1,    A2,    A3, # EVENTS\n");
	fprintf(fout, "  0.00000 0.00000 110.000 450.000  33.000  40.000 0.90000 2.30000\n");
	fprintf(fout, "  AFD7GH,ATS7GH,ZFD7GH,ZTS7GH,RCD7GH,ZCD7H,ASL7GH,ASU7GH\n");
	
	// Merge results from all threads
	std::cout << "Merging results from " << num_threads << " threads..." << std::endl;
	
	Int_t total_events = 0;
	for (const auto& thread_data : results) {
		for (size_t i = 0; i < thread_data.event_ids.size(); i++) {
			total_events++;
			
			// Fill tree
			part->Clear();
			new ((*part)[0]) PParticle(14, thread_data.px[i], thread_data.py[i], thread_data.pz[i]);
			new ((*part)[1]) PParticle(614, thread_data.cx[i], thread_data.cy[i], thread_data.cz[i], mC + E_EXCITATION, 1);
			tree->Fill();
			
			// Write to ASCII file
			fprintf(fout, " %10d%10d%10d%10.4f%10.3E%10.3E\n", 
			        total_events, codr, nbp, mom, wei, wei);
			fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 
			        1, thread_data.px[i], thread_data.py[i], thread_data.pz[i], 14);
			fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 
			        2, thread_data.cx[i], thread_data.cy[i], thread_data.cz[i], 614);
		}
	}
	
	std::cout << "Total events written: " << total_events << std::endl;
	
	// Close files
	fclose(fout);
	tree->Write();
	f.Close();
	
	std::cout << "\nGeneration time: " << elapsed.count() << " seconds" << std::endl;
	std::cout << "Events per second: " << number_total / (double)elapsed.count() << std::endl;
	std::cout << "========================================\n" << std::endl;
}





// ====================================================================
// INELASTIC SCATTERING
// ====================================================================

// ====================================================================
// GenerateThreadEventsInelasticPolarized - Worker function
// ====================================================================
ThreadData GenerateThreadEventsInelasticPolarized(int thread_id, Double_t ekin, 
                                                    Int_t events_per_thread,
                                                    Double_t theta_lab_target, 
                                                    Double_t theta_lab_window,
                                                    Double_t phi_lab_window,
                                                    Double_t polarization, 
                                                    Int_t spin_state)
{
    ThreadData data;
    data.count = 0;
    
    // Initial state kinematics
    Double_t ekinGeV = ekin/1000.;
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    
    Double_t s = iState.Mag2();
    
    // INELASTIC: Use reduced CM momentum
    Double_t pcm = CalculateCMMomentumInelastic(ekin, 0.00443);
    
    // Calculate CM angle range (silently - no printout)
    Double_t theta_cm_min, theta_cm_max;
    ComputeCMAngleRange(ekin, theta_lab_target, theta_lab_window, 
                        theta_cm_min, theta_cm_max);
    
    // INELASTIC: Create sampling histogram from digitized data
    // FIX: Give each thread a unique histogram name!
    TH1D* sigmacm_narrow = new TH1D(Form("h_inelastic_sampling_t%d", thread_id), 
                                     "Inelastic cross section sampling", 
                                     10000, theta_cm_min, theta_cm_max);
    
    // Fill histogram with interpolated cross section values
    for (Int_t ibin = 1; ibin <= 10000; ibin++) {
        Double_t theta_cm = sigmacm_narrow->GetBinCenter(ibin);
        Double_t xs = GetInelasticCrossSection(theta_cm);
        sigmacm_narrow->SetBinContent(ibin, xs);
    }
    
    // FIX: Set directory to 0 to prevent ROOT from managing it
    sigmacm_narrow->SetDirectory(0);
    
    // Detector acceptance in LAB
    const Double_t theta_lab_min = theta_lab_target - theta_lab_window;
    const Double_t theta_lab_max = theta_lab_target + theta_lab_window;
    const Double_t phi_lab_min_0   = -phi_lab_window;
    const Double_t phi_lab_max_0   = phi_lab_window;
    const Double_t phi_lab_min_180 = TMath::Pi() - phi_lab_window;
    const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;
    
    // Thread-local random number generator
    TRandom3 rng(thread_id + std::time(nullptr));
    
    int li = 0;
    while (li < events_per_thread)
    {
        // 1. Sample theta from inelastic cross section
        Double_t theta = sigmacm_narrow->GetRandom();
        Double_t stheta = TMath::Sin(TMath::DegToRad() * theta);
        Double_t ctheta = TMath::Cos(TMath::DegToRad() * theta);
        
        // 2. Sample phi from polarized distribution
        Double_t phi;
        if (polarization == 0.0) {
            // Unpolarized: uniform phi
            phi = rng.Rndm() * 2. * TMath::Pi();
        } else {
            // Polarized: sample from 1 + P·A_N·cos(phi)
            // INELASTIC: Use inelastic analyzing power!
            Double_t AN = GetInelasticAnalyzingPower(theta);
            Double_t asymmetry = spin_state * polarization * AN;
            
            // Acceptance-rejection sampling
            Double_t phi_test, weight;
            do {
                phi_test = rng.Rndm() * 2. * TMath::Pi();
                weight = 1.0 + asymmetry * TMath::Cos(phi_test);
            } while (rng.Rndm() > weight / (1.0 + TMath::Abs(asymmetry)));
            phi = phi_test;
        }
        
        // 3. Create 4-momentum in CM frame with INELASTIC pcm
        TLorentzVector p;
        p.SetXYZM(pcm * stheta * TMath::Cos(phi), 
                  pcm * stheta * TMath::Sin(phi), 
                  pcm * ctheta, 
                  mp);
        
        // 4. Boost to LAB frame
        p.Boost(iState.BoostVector());
        
        // 5. Calculate Carbon recoil
        TLorentzVector C = iState - p;
        
        // 6. Apply LAB acceptance cuts
        Bool_t theta_accepted = (p.Theta() >= theta_lab_min && p.Theta() <= theta_lab_max);
        Bool_t phi_accepted_0 = (p.Phi() >= phi_lab_min_0 && p.Phi() <= phi_lab_max_0);
        Bool_t phi_accepted_180 = (p.Phi() >= phi_lab_min_180 || p.Phi() <= phi_lab_max_180);
        
        if (!theta_accepted || !(phi_accepted_0 || phi_accepted_180)) {
            continue;
        }
        
        // 7. Store event
        li++;
        data.event_ids.push_back(li);
        data.px.push_back(p.X());
        data.py.push_back(p.Y());
        data.pz.push_back(p.Z());
        data.cx.push_back(C.X());
        data.cy.push_back(C.Y());
        data.cz.push_back(C.Z());
    }
    
    data.count = li;
    
    // Cleanup
    delete sigmacm_narrow;
    
    return data;
}

// ====================================================================
// SingleRunMultithreadInelasticPolarized - Main multithreaded function
// ====================================================================
void SingleRunMultithreadInelasticPolarized(Double_t energy, Int_t number_total,
                                             Double_t polarization, Int_t spin_state,
                                             Int_t num_threads,
                                             SamplingMode mode)
{
    std::cout << "\n=== Multithreaded POLARIZED INELASTIC Event Generation ===" << std::endl;
    std::cout << "Energy: " << energy << " MeV" << std::endl;
    std::cout << "Total events: " << number_total << std::endl;
    std::cout << "Excitation: 4.43 MeV (2+ state)" << std::endl;
    std::cout << "Polarization: " << polarization * 100 << "%" << std::endl;
    std::cout << "Spin state: " << (spin_state > 0 ? "UP" : "DOWN") << std::endl;
    std::cout << "Sampling mode: "
              << (mode == T_SAMPLING_MEYER ? "t-based Meyer (energy-dependent)"
                                           : "θ_CM-based (legacy)")
              << std::endl;

    // ── Energy range validation ────────────────────────────────
    if (mode == T_SAMPLING_MEYER && (energy < 160.0 || energy > 200.0)) {
        std::cout << "WARNING: T_SAMPLING_MEYER is validated for 160-200 MeV." << std::endl;
        std::cout << "         Energy " << energy << " MeV is outside this range." << std::endl;
        std::cout << "         Splines will clamp to the nearest boundary." << std::endl;
    }
    if (mode == THETA_CM_SAMPLING && energy != 200.0) {
        std::cout << "WARNING: THETA_CM_SAMPLING for INELASTIC uses CSV data" << std::endl;
        std::cout << "         measured at 200 MeV ONLY. At " << energy << " MeV," << std::endl;
        std::cout << "         the A_N will be wrong by up to 0.30." << std::endl;
        std::cout << "         Use T_SAMPLING_MEYER for correct energy dependence." << std::endl;
    }

    auto start = std::chrono::high_resolution_clock::now();

    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
    }
    std::cout << "Number of threads: " << num_threads << std::endl;

    Int_t events_per_thread = number_total / num_threads;

    const Double_t theta_lab_target = DETECTOR_THETA_CENTER_RAD;
    const Double_t theta_lab_window = DETECTOR_THETA_WINDOW;
    const Double_t phi_lab_window   = DETECTOR_PHI_WINDOW;

    ROOT::TThreadExecutor executor(num_threads);
    std::vector<ThreadData> results;

    if (mode == T_SAMPLING_MEYER) {
        // ============================================================
        // MEYER INELASTIC SAMPLING MODE
        // Build inelastic envelope once; share read-only across threads
        // ============================================================
        std::cout << "\nCalculating unified t range for inelastic Meyer envelope..." << std::endl;

        Double_t t_min, t_max;
        ComputeUnifiedTRange(energy, energy,
                             theta_lab_target, theta_lab_window,
                             t_min, t_max);

        std::cout << "  t range: [" << t_min << ", " << t_max << "] GeV²" << std::endl;
        std::cout << "Building Meyer inelastic envelope (50k bins)..." << std::endl;

        TH1D* meyer_envelope = CreateMeyerEnvelope_Inelastic(t_min, t_max);
        meyer_envelope->SetDirectory(0);
        // THREAD SAFETY: pre-build fIntegral before threads launch.
        meyer_envelope->ComputeIntegral();

        std::cout << "  Envelope max: " << meyer_envelope->GetMaximum()
                  << " mb/GeV²" << std::endl;

        auto thread_work = [energy, events_per_thread,
                            meyer_envelope,
                            theta_lab_target, theta_lab_window, phi_lab_window,
                            polarization, spin_state]
                           (int tid) {
            return GenerateThreadEventsInelastic_MeyerSampling(
                tid, energy, events_per_thread,
                meyer_envelope,
                theta_lab_target, theta_lab_window, phi_lab_window,
                polarization, spin_state);
        };

        std::cout << "Starting multithreaded Meyer inelastic generation..." << std::endl;
        results = executor.Map(thread_work, ROOT::TSeq<int>(0, num_threads));

        delete meyer_envelope;

    } else {
        // ============================================================
        // LEGACY THETA-BASED INELASTIC SAMPLING MODE
        // ============================================================
        // Load inelastic data if not already loaded
        if (!g_inelastic_xs || !g_inelastic_AN) {
            LoadInelasticData();
        }

        auto thread_work = [energy, events_per_thread, theta_lab_target,
                             theta_lab_window, phi_lab_window,
                             polarization, spin_state]
                            (int tid) {
            return GenerateThreadEventsInelasticPolarized(tid, energy, events_per_thread,
                                                           theta_lab_target, theta_lab_window,
                                                           phi_lab_window,
                                                           polarization, spin_state);
        };

        std::cout << "Starting multithreaded generation..." << std::endl;
        results = executor.Map(thread_work, ROOT::TSeq<int>(0, num_threads));
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    
    std::cout << "\n=== Writing to ROOT file ===" << std::endl;
    
    // Prepare output files with standardized names
    int polar_int = (int)(polarization * 100);
    TString spinLabel  = (spin_state > 0) ? "SpinUp" : "SpinDown";
    TString modeLabel  = (mode == T_SAMPLING_MEYER) ? "MEYER" : "THETA";
    TString basename = Form("pC_Inel443_%3.0fMeV_MT_P%d_%s_%dp%d_%s",
                        energy, polar_int, spinLabel.Data(),
                        (int)DETECTOR_THETA_CENTER, (int)(DETECTOR_THETA_CENTER*10)%10,
                        modeLabel.Data());
    TString nameo = MakeOutputPath(basename + ".txt");
    TString namer = MakeOutputPath(basename + ".root");

    std::cout << "Output files:" << std::endl;
    std::cout << "  ROOT: " << namer << std::endl;
    std::cout << "  TXT:  " << nameo << std::endl;
    std::cout << "========================================\n" << std::endl;
        
    // Create ROOT file and tree - STANDARD PLUTO FORMAT
    TFile f(namer.Data(), "RECREATE");
    TTree* tree = new TTree("data", "pC->pC* inelastic (4.43 MeV)");
    TClonesArray* part = new TClonesArray("PParticle");
    tree->Branch("Particles", &part, 32000, 0);
    
    // Standard PLUTO branches ONLY
    Int_t npart = 2;
    Float_t impact = 0.;
    Float_t phi0 = 0.;
    tree->Branch("Npart", &npart);
    tree->Branch("Impact", &impact);
    tree->Branch("Phi", &phi0);
    
    // Open ASCII output file
    FILE* fout = fopen(nameo.Data(), "w");
    
    // Write ASCII header
    Double_t ekinGeV = energy/1000.;
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    int codr = 302;  // 302 for inelastic
    int nbp = 2;
    double wei = 1.000E+00;
    int fform = 10000000 + codr;
    fprintf(fout, " %10d%10.2E%10.3f%10.3f%10.3f%10.3f%7d\n", 
            fform, 7.36E-8, mom, 0.f, 0.f, 0.f, number_total);
    fprintf(fout, " REAC,CROSS(mb),B. MOM,  A1,    A2,    A3, # EVENTS\n");
    fprintf(fout, "  0.00000 0.00000 110.000 450.000  33.000  40.000 0.90000 2.30000\n");
    fprintf(fout, "  AFD7GH,ATS7GH,ZFD7GH,ZTS7GH,RCD7GH,ZCD7H,ASL7GH,ASU7GH\n");
    
    // Merge results from all threads
    std::cout << "Merging results from " << num_threads << " threads..." << std::endl;
    
    Int_t total_events = 0;
    for (const auto& thread_data : results) {
        for (size_t i = 0; i < thread_data.event_ids.size(); i++) {
            total_events++;
            
            // Fill tree
            part->Clear();
            new ((*part)[0]) PParticle(14, thread_data.px[i], thread_data.py[i], thread_data.pz[i]);
            new ((*part)[1]) PParticle(614, thread_data.cx[i], thread_data.cy[i], thread_data.cz[i], mC + E_EXCITATION, 1);
            tree->Fill();
            
            // Write to ASCII file
            fprintf(fout, " %10d%10d%10d%10.4f%10.3E%10.3E\n", 
                    total_events, codr, nbp, mom, wei, wei);
            fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 
                    1, thread_data.px[i], thread_data.py[i], thread_data.pz[i], 14);
            fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 
                    2, thread_data.cx[i], thread_data.cy[i], thread_data.cz[i], 614);
        }
    }
    
    std::cout << "Total events written: " << total_events << std::endl;
    
    // Close files
    fclose(fout);
    tree->Write();
    f.Close();
    
    std::cout << "\nGeneration time: " << elapsed.count() << " seconds" << std::endl;
    std::cout << "Events per second: " << number_total / (double)elapsed.count() << std::endl;
    std::cout << "========================================\n" << std::endl;
}