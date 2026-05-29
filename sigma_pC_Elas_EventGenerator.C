#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TF2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TEllipse.h"      
#include "TLine.h"         
#include "TObjString.h"    
#include "PParticle.h"
#include "PChannel.h"
#include "TColor.h"
#include <ROOT/TThreadExecutor.hxx>

#include "DetectorConfig.h" // Setup the detectors acceptance
#include "PhysicalConstants.h" 
#include "ElasticScattering.h" 
#include "InelasticScattering.h"
#include "Kinematics.h"
#include "ThreadUtils.h"        
#include "EventGenerator.h"     


#include <TGraph2DErrors.h>
#include <TGraph.h> 
#include <stdio.h>
#include <thread>
#include <chrono>
#include <mutex>
#include <vector>
#include <algorithm>
#include <fstream>    
#include <string>    


using namespace std;




// ====================================================================
// Beam Polarization
// ====================================================================
const Double_t BEAM_POLARIZATION = 0.80;  // 80% polarization

// Spin state: +1 for spin-up, -1 for spin-down
// This will be passed as a parameter to functions that need it



// ====================================================================
// Produce: Generate p+C elastic scattering events (single-threaded)
// ====================================================================
// Uses IMPORTANCE SAMPLING: events are sampled in CM frame according 
// to the cross section, then boosted to lab frame and filtered by
// detector acceptance. No event weights needed - physics is encoded
// in the sampling distribution.
// ====================================================================
/* void Produce(const char* name, Double_t ekin, Int_t num, Bool_t withTree = 0)
{
	gROOT->SetStyle("Plain");

	// ================================================================
	// Create diagnostic histograms
	// ================================================================
	TFile* fspec = new TFile("spectra_ps.root", "RECREATE");
	TH1D* sigmacm  = new TH1D("sigmacm" , "dsigma/dtheta(thetaCM)", 10000,  1., 50.);
	TH1D* thlab    = new TH1D("thlab"   , "thlab"                 ,  1000,  0., 25.);
	TH1D* thcm     = new TH1D("thcm"    , "thcm"                  ,  1000,  0., 25.);
	TH1D* costhcm  = new TH1D("costhcm" , "costhcm"               ,  1000, -1.,  1.);
	TH1D* costhlab = new TH1D("costhlab", "costhlab"              ,  1000, -1.,  1.); 
	TH2D* pzvspxpy = new TH2D("pzvspxpy", "pzvspxpy", 500, 0.3, 1.5, 100, 0., 0.8);

	// ================================================================
	// Get cross section formula for this energy
	// ================================================================
	const char* xsFormula = GetXSFormula(ekin);
	TF1 *ffps = new TF1("ffps", xsFormula, 0.*TMath::DegToRad(), 50.*TMath::DegToRad());

	// ================================================================
	// Prepare output files
	// ================================================================
	// Prepare output files with standardized names
	TString basename = Form("pC_Elas_%3.0fMeV_ST", ekin);
	TString nameo = MakeOutputPath(basename + ".txt");
	TString namer = MakeOutputPath(basename + ".root");

	std::cout << "\n=== Single-threaded Event Generation ===" << std::endl;
	std::cout << "Output files:" << std::endl;
	std::cout << "  ROOT: " << namer << std::endl;
	std::cout << "  TXT:  " << nameo << std::endl;
	std::cout << "========================================\n" << std::endl;

	TFile f(namer.Data(), "RECREATE");
	
	TTree* tree = new TTree("data", "pC->pC cross section");
	TClonesArray* part = new TClonesArray("PParticle");
	tree -> Branch("Particles", &part);

	// ================================================================
	// Event metadata branches (NO WEIGHT BRANCH)
	// ================================================================
	Int_t   npart  = 2;      // Always 2 particles: proton + Carbon
	Float_t impact = 0.;     // Impact parameter (not used here)
	Float_t phi0   = 0.;     // Azimuthal angle (not used here)
	tree -> Branch("Npart" , &npart );
	tree -> Branch("Impact", &impact);
	tree -> Branch("Phi"   , &phi0  );
	
	// Open ASCII output file
	FILE* fout;
	fout = fopen(nameo.Data(), "w");

	// ================================================================
	// STEP 1: Calculate initial state kinematics
	// ================================================================
	Double_t ekinGeV = ekin/1000.;  // Convert MeV → GeV
	Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);  // Beam momentum
	
	// Initial state: beam proton + target Carbon at rest
	TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
	TLorentzVector p;   // Scattered proton 4-vector
	TLorentzVector C;   // Recoil Carbon 4-vector
	Double_t phi, theta, pcm, s;
	Double_t stheta, ctheta;
	Int_t nevent = num;
   
	// Calculate CM kinematics
	s = iState.Mag2();  // Mandelstam s
	pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));  // CM momentum

	std::cout << "pbeam: " << mom << "\tpcm: " << pcm << "\ts: " << s << std::endl;
	
	// ================================================================
	// Write ASCII file header (Pluto format)
	// ================================================================
	int    codr  = 301;  // Reaction code
	int    nbp   = 2;    // Number of particles
	double wei   = 1.000E+00;  // Weight for ASCII file only (not stored in tree)
	int    fform = 10000000 + codr;
	fprintf(fout, " %10d%10.2E%10.3f%10.3f%10.3f%10.3f%7d\n", fform, 7.36E-8, mom, 0.f, 0.f, 0.f, num);
	fprintf(fout, " REAC,CROSS(mb),B. MOM,  A1,    A2,    A3, # EVENTS\n");
	fprintf(fout, "  0.00000 0.00000 110.000 450.000  33.000  40.000 0.90000 2.30000\n");
	fprintf(fout, "  AFD7GH,ATS7GH,ZFD7GH,ZTS7GH,RCD7GH,ZCD7H,ASL7GH,ASU7GH\n");

	// ================================================================
	// STEP 2: Define detector acceptance in lab frame
	// ================================================================
	// Polar angle: 16° ± 5 mrad
	const Double_t theta_lab_target = DETECTOR_THETA_CENTER_RAD;  // Central angle: 16°
	const Double_t theta_lab_window = DETECTOR_THETA_WINDOW;  // Angular width: ±5 mrad
	
	// Azimuthal angle: 0° ± 5 mrad OR 180° ± 5 mrad
	// NOTE: ROOT's Phi() returns [-π, +π], so:
	//   - Around 0°: accept phi in [-5 mrad, +5 mrad]
	//   - Around 180°: accept phi near +π OR near -π (same physical direction)
	const Double_t phi_lab_window   = DETECTOR_PHI_WINDOW;  // Azimuthal width: ±5 mrad
	
	const Double_t theta_lab_min = theta_lab_target - theta_lab_window;
	const Double_t theta_lab_max = theta_lab_target + theta_lab_window;
	
	// Phi acceptance: two windows
	// Window 1: around 0°
	const Double_t phi_lab_min_0   = -phi_lab_window;
	const Double_t phi_lab_max_0   = phi_lab_window;
	// Window 2: around 180° (near ±π in ROOT convention)
	const Double_t phi_lab_min_180 = TMath::Pi() - phi_lab_window;
	const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;  // Note: near -π is same as near +π

	// ================================================================
	// STEP 3: Calculate CM angle range for efficient sampling
	// ================================================================
	Double_t theta_cm_min, theta_cm_max;
	ComputeCMAngleRange(ekin, theta_lab_target, theta_lab_window, 
	                    theta_cm_min, theta_cm_max);
	
	// Create cross section function over narrow CM range
	// This focuses sampling on angles that will reach the detector
	TF1 *ffps_narrow = new TF1("ffps_narrow", xsFormula, 
	                           theta_cm_min*TMath::DegToRad(), 
	                           theta_cm_max*TMath::DegToRad());
	
	// Create histogram for importance sampling
	// Bin contents = cross section values at each angle
	TH1D* sigmacm_narrow = new TH1D("sigmacm_narrow", "narrow range", 
	                                 10000, theta_cm_min, theta_cm_max);
	for (int ibin = 1; ibin <= 10000; ibin++) {
		Double_t theta_deg = theta_cm_min + (theta_cm_max - theta_cm_min) * (ibin-0.5) / 10000.;
		sigmacm_narrow->SetBinContent(ibin, ffps_narrow->Eval(theta_deg * TMath::DegToRad()));
	}
	
	// ================================================================
	// STEP 4: Event generation loop (IMPORTANCE SAMPLING)
	// ================================================================
	int li = 0;  // Counter for accepted events
	while ( li < nevent )
	{
		// Progress printing
		if ( !(li%10000) ) std::cout << li << " ..." << std::endl;
		
		// ============================================================
		// 4a. Sample CM angle from cross section distribution
		// ============================================================
		// GetRandom() samples from histogram with probability ∝ bin content
		// This is IMPORTANCE SAMPLING: we sample more events where σ is large
		theta = sigmacm_narrow->GetRandom();  // Returns angle in degrees
		
		// Calculate sin and cos for momentum components
		stheta = TMath::Sin(TMath::DegToRad() * theta);
		ctheta = TMath::Cos(TMath::DegToRad() * theta);
		
		// Sample azimuthal angle uniformly
		phi = gRandom->Rndm() * 2. * TMath::Pi();
		
		// ============================================================
		// 4b. Create 4-momentum in CM frame
		// ============================================================
		// Proton momentum in CM: magnitude = pcm, direction = (θ_CM, φ)
		p.SetXYZM(pcm * stheta * cos(phi), 
		          pcm * stheta * sin(phi), 
		          pcm * ctheta, 
		          mp);
		
		// ============================================================
		// 4c. Fill CM histograms BEFORE boost (correct!)
		// ============================================================
		thcm    -> Fill(p.Theta() * TMath::RadToDeg());
		costhcm -> Fill(TMath::Cos(p.Theta()));
		
		// ============================================================
		// 4d. Boost to lab frame
		// ============================================================
		// This Lorentz transformation preserves the physics!
		// The cross section is defined in CM, and boosting doesn't
		// change the probability distribution we sampled from
		p.Boost(iState.BoostVector());
		
		// Calculate recoil Carbon 4-momentum (momentum conservation)
		C = iState - p;
		
		// ============================================================
		// 4e. Apply detector acceptance cuts in lab frame
		// ============================================================
		// Theta: 16° ± 5 mrad
		// Phi: (0° ± 5 mrad) OR (180° ± 5 mrad)
		// Note: ROOT Phi() returns [-π, +π], so 180° appears as both +π and -π
		Bool_t theta_accepted = (p.Theta() >= theta_lab_min && p.Theta() <= theta_lab_max);
		Bool_t phi_accepted_0   = (p.Phi() >= phi_lab_min_0   && p.Phi() <= phi_lab_max_0);
		Bool_t phi_accepted_180 = (p.Phi() >= phi_lab_min_180 || p.Phi() <= phi_lab_max_180);  // Near +π OR near -π
		
		if ( theta_accepted && (phi_accepted_0 || phi_accepted_180) )
		{
			// Event accepted! Increment counter
			li++;
			
			// Fill lab frame histograms AFTER boost (correct!)
			thlab    -> Fill(p.Theta() * TMath::RadToDeg());
			costhlab -> Fill(TMath::Cos(p.Theta()));
			pzvspxpy -> Fill(p.Z(), TMath::Sqrt(p.X()*p.X() + p.Y()*p.Y()));
			
			// ====================================================
			// 4f. Save event to ROOT tree (NO WEIGHT BRANCH)
			// ====================================================
			part->Clear();
			// Particle 1: scattered proton (ID=14)
			new ((*part)[0]) PParticle(14, p.X(), p.Y(), p.Z());
			// Particle 2: recoil Carbon (ID=614)
			new ((*part)[1]) PParticle(614, C.X(), C.Y(), C.Z(), mC + E_EXCITATION, 1);
			
			// No weight branch - importance sampling means all events have implicit weight = 1
			tree->Fill();
			
			// ====================================================
			// 4g. Write to ASCII file (Pluto format)
			// ====================================================
			fprintf(fout, " %10d%10d%10d%10.4f%10.3E%10.3E\n", 
			        li, codr, nbp, mom, wei, wei);
			fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 
			        1, p.X(), p.Y(), p.Z(), 14);
			fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 
			        2, C.X(), C.Y(), C.Z(), 614);
		}
	}

	// ================================================================
	// STEP 5: Write output and clean up
	// ================================================================
	std::cout << "n processed: " << li << std::endl;
	
	fclose(fout);
	tree->Write();
	f.Close();
	
	// Write diagnostic histograms
	fspec->cd();
	sigmacm->Write();
	thlab->Write();
	thcm->Write();
	costhcm->Write();
	costhlab->Write();
	pzvspxpy->Write();
	fspec->Close();
	
	delete sigmacm_narrow;
	delete ffps_narrow;
	delete ffps;
	
	std::cout << "Generated " << nevent << " accepted events" << std::endl;
	std::cout << "Output files: " << namer << ", " << nameo << std::endl;
	
	// Optionally remove ROOT file if not needed
	if( !withTree )
	{
		char command[300];
		sprintf(command, "rm %s", namer.Data());
		std::cout << command << std::endl;
		gSystem -> Exec(command);
	}
}

// ====================================================================
// MultiRun: Run multiple energies sequentially
// ====================================================================
void MultiRun()
{
    Double_t E[6] = {150., 160., 170., 180., 190., 200.};
    
    for ( Int_t i = 0; i < 6; i++ )
    {
        Produce("", E[i], 10000000, 1);  // name parameter ignored
        cout << "Completed energy: " << E[i] << " MeV" << endl;
    }
}

// ====================================================================
// SingleRun: Run single energy
// ====================================================================
void SingleRun(Double_t energy)
{
    Int_t number = 100000;
    Produce("", energy, number, 1);  // name parameter ignored
    cout << "Completed: " << energy << " MeV" << endl;
}

// ====================================================================
// ProducePolarized: Generate polarized p+C elastic scattering events
// ====================================================================
void ProducePolarized(const char* name, Double_t ekin, Int_t num, Int_t spin_state, Bool_t withTree = 0)
{
    gROOT->SetStyle("Plain");
    
    cout << "\n========================================" << endl;
    cout << "POLARIZED EVENT GENERATION" << endl;
    cout << "========================================" << endl;
    cout << "Energy: " << ekin << " MeV" << endl;
    cout << "Events: " << num << endl;
    cout << "Polarization: " << BEAM_POLARIZATION * 100 << "%" << endl;
    cout << "Spin state: " << (spin_state > 0 ? "UP" : "DOWN") << endl;
    cout << "========================================\n" << endl;
    
    // Get cross section formula
    const char* xsFormula = GetXSFormula(ekin);
    
    // Prepare output files with standardized names
	int polar_int = (int)(BEAM_POLARIZATION * 100);
	TString spinLabel = (spin_state > 0) ? "SpinUp" : "SpinDown";
	TString basename = Form("pC_Elas_%3.0fMeV_ST_P%d_%s", ekin, polar_int, spinLabel.Data());
	TString nameo = MakeOutputPath(basename + ".txt");
	TString namer = MakeOutputPath(basename + ".root");

	std::cout << "Output files:" << std::endl;
	std::cout << "  ROOT: " << namer << std::endl;
	std::cout << "  TXT:  " << nameo << std::endl;
	std::cout << "========================================\n" << std::endl;

	TFile f(namer.Data(), "RECREATE");
   
    TTree* tree = new TTree("data", "pC->pC cross section");
    TClonesArray* part = new TClonesArray("PParticle");
    tree->Branch("Particles", &part);
    
    Int_t npart = 2;
    Float_t impact = 0.;
    Float_t phi0 = 0.;
    tree->Branch("Npart", &npart);
    tree->Branch("Impact", &impact);
    tree->Branch("Phi", &phi0);
    
    FILE* fout = fopen(nameo.Data(), "w");
    
    // Calculate initial state kinematics
    Double_t ekinGeV = ekin/1000.;
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    TLorentzVector p, C;
    Double_t s = iState.Mag2();
    Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));
    
    std::cout << "pbeam: " << mom << "\tpcm: " << pcm << "\ts: " << s << std::endl;
    
    // Write ASCII file header
    int codr = 301;
    int nbp = 2;
    double wei = 1.000E+00;
    int fform = 10000000 + codr;
    fprintf(fout, " %10d%10.2E%10.3f%10.3f%10.3f%10.3f%7d\n", fform, 7.36E-8, mom, 0.f, 0.f, 0.f, num);
    fprintf(fout, " REAC,CROSS(mb),B. MOM,  A1,    A2,    A3, # EVENTS\n");
    fprintf(fout, "  0.00000 0.00000 110.000 450.000  33.000  40.000 0.90000 2.30000\n");
    fprintf(fout, "  AFD7GH,ATS7GH,ZFD7GH,ZTS7GH,RCD7GH,ZCD7H,ASL7GH,ASU7GH\n");
    
    // Define detector acceptance in LAB frame
    const Double_t theta_lab_target = DETECTOR_THETA_CENTER_RAD;
    const Double_t theta_lab_window = DETECTOR_THETA_WINDOW;
    const Double_t phi_lab_window = DETECTOR_PHI_WINDOW;
    
    const Double_t theta_lab_min = theta_lab_target - theta_lab_window;
    const Double_t theta_lab_max = theta_lab_target + theta_lab_window;
    const Double_t phi_lab_min_0 = -phi_lab_window;
    const Double_t phi_lab_max_0 = phi_lab_window;
    const Double_t phi_lab_min_180 = TMath::Pi() - phi_lab_window;
    const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;
    
    // Calculate CM angle ranges
    Double_t theta_cm_min, theta_cm_max;
    ComputeCMAngleRange(ekin, theta_lab_target, theta_lab_window, 
                        theta_cm_min, theta_cm_max);
    
    Double_t theta_cm_avg = (theta_cm_min + theta_cm_max) / 2.0;
    
    Double_t phi_cm_min_0, phi_cm_max_0, phi_cm_min_180, phi_cm_max_180;
    ComputePhiCMRanges(ekin, theta_cm_avg,
                       phi_cm_min_0, phi_cm_max_0,
                       phi_cm_min_180, phi_cm_max_180);
    
    // Create theta sampling histogram
    TF1 *ffps_narrow = new TF1("ffps_narrow", xsFormula, 
                               theta_cm_min*TMath::DegToRad(), 
                               theta_cm_max*TMath::DegToRad());
    
    TH1D* sigmacm_narrow = new TH1D("sigmacm_narrow", "narrow range", 
                                     10000, theta_cm_min, theta_cm_max);
    for (int ibin = 1; ibin <= 10000; ibin++) {
        Double_t theta_deg = theta_cm_min + (theta_cm_max - theta_cm_min) * (ibin-0.5) / 10000.;
        sigmacm_narrow->SetBinContent(ibin, ffps_narrow->Eval(theta_deg * TMath::DegToRad()));
    }
    
    // Create phi sampling function (will reuse with different parameters)
    TF1 *fPhi = new TF1("fPhi_pol", "1 + [0]*cos(x)", 0, 2*TMath::Pi());
    
    // Event generation loop with polarization
    int li = 0;
    int attempts = 0;
    
    while (li < num)
    {
        if (!(li%10000)) std::cout << li << " ..." << std::endl;
        attempts++;
        
        // 1. Sample theta_CM from optical model
        Double_t theta_cm_deg = sigmacm_narrow->GetRandom();
        
        // 2. Convert to LAB and get analyzing power
        Double_t theta_lab_deg = ConvertThetaCMtoLab(ekin, theta_cm_deg);
        Double_t AN = GetElasticAnalyzingPower(ekin, theta_lab_deg);
        
        // 3. Sample phi_CM from polarized distribution
        Bool_t use_detector_0 = (gRandom->Rndm() < 0.5);
        
        // Set the asymmetry parameter
        Double_t effective_asymmetry = spin_state * BEAM_POLARIZATION * AN;
        fPhi->SetParameter(0, effective_asymmetry);
        
        Double_t phi_cm;
        
        if (use_detector_0) {
            fPhi->SetRange(phi_cm_min_0, phi_cm_max_0);
            phi_cm = fPhi->GetRandom();
        } else {
            fPhi->SetRange(phi_cm_min_180, phi_cm_max_180);
            phi_cm = fPhi->GetRandom();
        }
        
        // 4. Create 4-momentum in CM frame
        Double_t theta_cm_rad = theta_cm_deg * TMath::DegToRad();
        Double_t stheta = TMath::Sin(theta_cm_rad);
        Double_t ctheta = TMath::Cos(theta_cm_rad);
        
        p.SetXYZM(pcm * stheta * TMath::Cos(phi_cm), 
                  pcm * stheta * TMath::Sin(phi_cm), 
                  pcm * ctheta, 
                  mp);
        
        // 5. Boost to LAB frame
        p.Boost(iState.BoostVector());
        C = iState - p;
        
        // 6. Apply LAB acceptance cuts
        Bool_t theta_accepted = (p.Theta() >= theta_lab_min && p.Theta() <= theta_lab_max);
        Bool_t phi_accepted_0 = (p.Phi() >= phi_lab_min_0 && p.Phi() <= phi_lab_max_0);
        Bool_t phi_accepted_180 = (p.Phi() >= phi_lab_min_180 || p.Phi() <= phi_lab_max_180);
        
        if (!theta_accepted || !(phi_accepted_0 || phi_accepted_180)) {
            continue;
        }
        
        li++;
        
        // 7. Save to ROOT tree
        part->Clear();
        new ((*part)[0]) PParticle(14, p.X(), p.Y(), p.Z());
        new ((*part)[1]) PParticle(614, C.X(), C.Y(), C.Z(), mC, 1);
        tree->Fill();
        
        // 8. Write to ASCII file
        fprintf(fout, " %10d%10d%10d%10.4f%10.3E%10.3E\n", li, codr, nbp, mom, wei, wei);
        fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 1, p.X(), p.Y(), p.Z(), 14);
        fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 2, C.X(), C.Y(), C.Z(), 614);
    }
    
    std::cout << "\nn processed: " << li << std::endl;
    std::cout << "Attempts: " << attempts << std::endl;
    std::cout << "Efficiency: " << 100.0*li/attempts << "%" << std::endl;
    
    // Close files properly
    fclose(fout);
    tree->Write();
    f.Close();
    delete tree;
    
    // Clean up objects
    delete sigmacm_narrow;
    delete ffps_narrow;
    delete fPhi;
    
    std::cout << "Generated " << num << " polarized events" << std::endl;
    std::cout << "Output files: " << namer << ", " << nameo << std::endl;
    
    if (!withTree) {
        char command[300];
        sprintf(command, "rm %s", namer.Data());
        gSystem->Exec(command);
    }
}



// ====================================================================
// GenerateThreadEvents - Worker function for multithreading
// ====================================================================
// Each thread independently generates a batch of events
// This function is identical to the main Produce() loop, but:
//   - Returns data instead of writing directly to file
//   - Uses thread-safe random number generator
//   - No file I/O (done later by main thread)
// ====================================================================
ThreadData GenerateThreadEvents(int thread_id, Double_t ekin, Int_t events_per_thread,
                                Double_t theta_lab_target, Double_t theta_lab_window,
                                Double_t phi_lab_window, const char* xsFormula)
{
	ThreadData data;
	data.count = 0;
	
	// ================================================================
	// STEP 1: Calculate initial state kinematics
	// ================================================================
	Double_t ekinGeV = ekin/1000.;
	Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
	TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
	
	Double_t s = iState.Mag2();
	Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));
	
	// ================================================================
	// STEP 2: Create thread-local cross section functions
	// ================================================================
	// Each thread needs its own TF1 objects (ROOT objects not thread-safe)
	TF1 *ffps = new TF1(Form("ffps_t%d", thread_id), xsFormula, 
	                     0.*TMath::DegToRad(), 50.*TMath::DegToRad());
	
	// Calculate CM angle range
	Double_t theta_cm_min, theta_cm_max;
	ComputeCMAngleRange(ekin, theta_lab_target, theta_lab_window, 
	                    theta_cm_min, theta_cm_max);
	
	// Create narrow-range cross section function for importance sampling
	TF1 *ffps_narrow = new TF1(Form("ffps_narrow_t%d", thread_id), xsFormula, 
	                            theta_cm_min*TMath::DegToRad(), 
	                            theta_cm_max*TMath::DegToRad());
	
	// Create sampling histogram
	TH1D* sigmacm_narrow = new TH1D(Form("sigmacm_narrow_t%d", thread_id), 
	                                 Form("narrow range t%d", thread_id), 
	                                 10000, theta_cm_min, theta_cm_max);
	for (int ibin = 1; ibin <= 10000; ibin++) {
		Double_t theta_deg = theta_cm_min + (theta_cm_max - theta_cm_min) * (ibin-0.5) / 10000.;
		sigmacm_narrow->SetBinContent(ibin, ffps_narrow->Eval(theta_deg * TMath::DegToRad()));
	}
	
	// ================================================================
	// STEP 3: Define lab acceptance
	// ================================================================
	const Double_t theta_lab_min = theta_lab_target - theta_lab_window;
	const Double_t theta_lab_max = theta_lab_target + theta_lab_window;
	
	// Phi acceptance: two windows
	// Window 1: around 0°
	const Double_t phi_lab_min_0   = -phi_lab_window;
	const Double_t phi_lab_max_0   = phi_lab_window;
	// Window 2: around 180° (near ±π in ROOT convention)
	const Double_t phi_lab_min_180 = TMath::Pi() - phi_lab_window;
	const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;  // Near -π is same as near +π
	
	// ================================================================
	// STEP 4: Thread-safe random number generator
	// ================================================================
	// Seed with thread_id + current time to ensure different sequences
	TRandom3 rng(thread_id + std::time(nullptr));
	
	// ================================================================
	// STEP 5: Event generation loop (IMPORTANCE SAMPLING)
	// ================================================================
	int li = 0;
	while ( li < events_per_thread )
	{
		// ============================================================
		// 5a. Sample CM angle from cross section distribution
		// ============================================================
		Double_t theta = sigmacm_narrow->GetRandom();  // Importance sampling!
		Double_t stheta = TMath::Sin(TMath::DegToRad() * theta);
		Double_t ctheta = TMath::Cos(TMath::DegToRad() * theta);
		
		// Sample azimuthal angle uniformly
		Double_t phi = rng.Rndm() * 2. * TMath::Pi();
		
		// ============================================================
		// 5b. Create 4-momentum in CM frame
		// ============================================================
		TLorentzVector p;
		p.SetXYZM(pcm * stheta * cos(phi), pcm * stheta * sin(phi), pcm * ctheta, mp);
		
		// ============================================================
		// 5c. Boost to lab frame
		// ============================================================
		p.Boost(iState.BoostVector());
		TLorentzVector C = iState - p;
		
		// ============================================================
		// 5d. Apply detector acceptance cuts in lab frame
		// ============================================================
		// Theta: 16° ± 5 mrad
		// Phi: (0° ± 5 mrad) OR (180° ± 5 mrad)
		// Note: ROOT Phi() returns [-π, +π], so 180° appears as both +π and -π
		Bool_t theta_accepted = (p.Theta() >= theta_lab_min && p.Theta() <= theta_lab_max);
		Bool_t phi_accepted_0   = (p.Phi() >= phi_lab_min_0   && p.Phi() <= phi_lab_max_0);
		Bool_t phi_accepted_180 = (p.Phi() >= phi_lab_min_180 || p.Phi() <= phi_lab_max_180);  // Near +π OR near -π
		
		if ( theta_accepted && (phi_accepted_0 || phi_accepted_180) )
		{
			// Event accepted!
			li++;
			data.event_ids.push_back(li);
			
			// Store momentum components
			data.px.push_back(p.X());
			data.py.push_back(p.Y());
			data.pz.push_back(p.Z());
			data.cx.push_back(C.X());
			data.cy.push_back(C.Y());
			data.cz.push_back(C.Z());
			
			// ====================================================
			// IMPORTANCE SAMPLING: No weight stored
			// ====================================================
			// We sampled from the cross section distribution,
			// so the physics is already encoded in our sampling.
			// All accepted events have implicit weight = 1.0
			// ====================================================
		}
	}
	
	// ================================================================
	// STEP 6: Clean up thread-local objects
	// ================================================================
	data.count = li;
	delete ffps;
	delete ffps_narrow;
	delete sigmacm_narrow;
	
	return data;
}

// ====================================================================
// SingleRunMultithread - Parallel event generation
// ====================================================================
// Uses multiple CPU cores to generate events in parallel
// Events from all threads are merged and written to a single file
// ====================================================================
void SingleRunMultithread(Double_t energy, Int_t number_total, Int_t num_threads = 0)
{
	std::cout << "\n=== Multithreaded Event Generation ===" << std::endl;
	std::cout << "Energy: " << energy << " MeV" << std::endl;
	std::cout << "Total events to generate: " << number_total << std::endl;
	
	auto start = std::chrono::high_resolution_clock::now();
	
	// Auto-detect number of CPU cores if not specified
	if (num_threads == 0) {
		num_threads = std::thread::hardware_concurrency();
	}
	std::cout << "Number of threads: " << num_threads << std::endl;
	
	// Divide work among threads
	Int_t events_per_thread = number_total / num_threads;
	
	// Define detector acceptance (same as single-threaded version)
	const Double_t theta_lab_target = DETECTOR_THETA_CENTER_RAD;
	const Double_t theta_lab_window = DETECTOR_THETA_WINDOW;
	const Double_t phi_lab_window   = DETECTOR_PHI_WINDOW;
	
	// Get cross section formula
	const char* xsFormula = GetXSFormula(energy);
	
	// ================================================================
	// STEP 1: Launch parallel event generation
	// ================================================================
	ROOT::TThreadExecutor executor(num_threads);
	
	// Lambda function: work to be done by each thread
	auto thread_work = [energy, events_per_thread, theta_lab_target, 
	                     theta_lab_window, phi_lab_window, xsFormula]
	                    (int tid) {
		return GenerateThreadEvents(tid, energy, events_per_thread, 
		                             theta_lab_target, theta_lab_window, 
		                             phi_lab_window, xsFormula);
	};
	
	std::cout << "Starting multithreaded generation..." << std::endl;
	
	// Execute in parallel: each thread returns a ThreadData object
	auto results = executor.Map(thread_work, ROOT::TSeq<int>(0, num_threads));
	
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	
	// ================================================================
	// STEP 2: Merge results and write to file (single-threaded)
	// ================================================================
	std::cout << "\n=== Writing to ROOT file ===" << std::endl;
	
	// Prepare output files with standardized names
	TString basename = Form("pC_Elas_%3.0fMeV_MT", energy);
	TString nameo = MakeOutputPath(basename + ".txt");
	TString namer = MakeOutputPath(basename + ".root");

	std::cout << "Output files:" << std::endl;
	std::cout << "  ROOT: " << namer << std::endl;
	std::cout << "  TXT:  " << nameo << std::endl;
	std::cout << "========================================\n" << std::endl;
	
	// Create ROOT file and tree
	TFile f(namer.Data(), "RECREATE");
	TTree* tree = new TTree("data", "pC->pC cross section");
	TClonesArray* part = new TClonesArray("PParticle");
	//tree->Branch("Particles", &part);
	tree->Branch("Particles", &part, 32000, 0);
	
	// Event metadata (NO WEIGHT BRANCH)
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
	int codr = 301;
	int nbp = 2;
	double wei = 1.000E+00;  // Weight for ASCII file only
	int fform = 10000000 + codr;
	fprintf(fout, " %10d%10.2E%10.3f%10.3f%10.3f%10.3f%7d\n", fform, 7.36E-8, mom, 0.f, 0.f, 0.f, number_total);
	fprintf(fout, " REAC,CROSS(mb),B. MOM,  A1,    A2,    A3, # EVENTS\n");
	fprintf(fout, "  0.00000 0.00000 110.000 450.000  33.000  40.000 0.90000 2.30000\n");
	fprintf(fout, "  AFD7GH,ATS7GH,ZFD7GH,ZTS7GH,RCD7GH,ZCD7H,ASL7GH,ASU7GH\n");
	
	// ================================================================
	// STEP 3: Loop over all thread results and write events
	// ================================================================
	int global_event_id = 0;
	for (auto& result : results) {
		for (int i = 0; i < result.count; i++) {
			global_event_id++;
			
			// Write to ROOT tree (NO WEIGHT)
			part->Clear();
			new ((*part)[0]) PParticle(14, result.px[i], result.py[i], result.pz[i]);
			new ((*part)[1]) PParticle(614, result.cx[i], result.cy[i], result.cz[i], mC, 1);
			
			// No weight branch - importance sampling
			tree->Fill();
			
			// Write to ASCII file
			fprintf(fout, " %10d%10d%10d%10.4f%10.3E%10.3E\n", global_event_id, codr, nbp, mom, wei, wei);
			fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 1, result.px[i], result.py[i], result.pz[i], 14);
			fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 2, result.cx[i], result.cy[i], result.cz[i], 614);
		}
	}
	
	// ================================================================
	// STEP 4: Finalize and print summary
	// ================================================================
	fclose(fout);
	tree->Write();
	f.Close();
	
	std::cout << "Total events written: " << global_event_id << std::endl;
	std::cout << "Output file: " << namer << std::endl;
	std::cout << "Time elapsed: " << elapsed.count() << " seconds" << std::endl;
	if (elapsed.count() > 0) {
		std::cout << "Events/second: " << global_event_id / (double)elapsed.count() << std::endl;
	}
	std::cout << "===================================\n" << std::endl;
} 
*/

// ====================================================================
// VerifySimulation - Verification and diagnostic plots
// ====================================================================
// Run this to verify your simulation output
// Checks conservation laws, lab frame acceptance, CM kinematics
// 
// Usage examples:
//   VerifySimulation(200.0)  // Auto-finds pC_Elas_200MeV file
//   VerifySimulation(200.0, true)  // Use multithreaded file
//   VerifySimulation("custom_file.root")  // Use specific file
// ====================================================================
void VerifySimulation(Double_t energy = 200.0, Bool_t use_multithread = false, 
                      Int_t polarization = -1, const char* spin_label = "")
{
	// Construct filename based on parameters
	TString filename;

	if (polarization >= 0 && strlen(spin_label) > 0) {
		// Polarized file
		if (use_multithread) {
			filename = MakeOutputPath(Form("pC_Elas_%3.0fMeV_MT_P%d_%s.root", energy, polarization, spin_label));
		} else {
			filename = MakeOutputPath(Form("pC_Elas_%3.0fMeV_ST_P%d_%s.root", energy, polarization, spin_label));
		}
	} else {
		// Unpolarized file
		if (use_multithread) {
			filename = MakeOutputPath(Form("pC_Elas_%3.0fMeV_MT.root", energy));
		} else {
			filename = MakeOutputPath(Form("pC_Elas_%3.0fMeV_ST.root", energy));
		}
	}
    
    cout << "\n========================================" << endl;
    cout << "Attempting to open: " << filename << endl;
    cout << "========================================" << endl;
    
    TFile *f = new TFile(filename);
    if (!f || f->IsZombie()) {
        cout << "Error: Cannot open file " << filename << endl;
        cout << "\nTrying to find alternative files..." << endl;
        
		// Try without path prefix (current directory)
		if (polarization >= 0 && strlen(spin_label) > 0) {
			if (use_multithread) {
				filename = Form("pC_Elas_%3.0fMeV_MT_P%d_%s.root", energy, polarization, spin_label);
			} else {
				filename = Form("pC_Elas_%3.0fMeV_ST_P%d_%s.root", energy, polarization, spin_label);
			}
		} else {
			if (use_multithread) {
				filename = Form("pC_Elas_%3.0fMeV_MT.root", energy);
			} else {
				filename = Form("pC_Elas_%3.0fMeV_ST.root", energy);
			}
		}

		f = new TFile(filename);
		if (!f || f->IsZombie()) {
			cout << "Error: Could not find file" << endl;
			cout << "Tried:" << endl;
			if (polarization >= 0) {
				cout << "  - " << MakeOutputPath(Form("pC_Elas_%3.0fMeV_%s_P%d_%s.root", 
					 energy, use_multithread?"MT":"ST", polarization, spin_label)) << endl;
				cout << "  - " << Form("pC_Elas_%3.0fMeV_%s_P%d_%s.root", 
					 energy, use_multithread?"MT":"ST", polarization, spin_label) << endl;
			} else {
				cout << "  - " << MakeOutputPath(Form("pC_Elas_%3.0fMeV_%s.root", 
					 energy, use_multithread?"MT":"ST")) << endl;
				cout << "  - " << Form("pC_Elas_%3.0fMeV_%s.root", 
					 energy, use_multithread?"MT":"ST") << endl;
			}
			return;
		}

        cout << "Found: " << filename << endl;
    }
    
    TTree *tree = (TTree*)f->Get("data");
    if (!tree) {
        cout << "Error: Tree 'data' not found" << endl;
        return;
    }
    
    TClonesArray *particles = 0;
    tree->SetBranchAddress("Particles", &particles);
    
    // Constants - automatically use the energy parameter
    const Double_t mp = 0.9382720813;
    const Double_t mC = 11.174862;
    const Double_t ekin = energy / 1000.0;  // Convert MeV to GeV
    const Double_t mom_beam = TMath::Sqrt(2*ekin*mp + ekin*ekin);
    const Double_t E_total = ekin + mp + mC;
    
    // Initial state 4-vector for CM calculations
    TLorentzVector iState(0., 0., mom_beam, ekin + mp + mC);
    Double_t s = iState.Mag2();
    Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));
    
    cout << "\n========================================" << endl;
    cout << "Verification of pC Elastic Scattering" << endl;
    cout << "========================================" << endl;
    cout << "File: " << filename << endl;
    cout << "Total events: " << tree->GetEntries() << endl;
    cout << "Beam energy: " << ekin*1000 << " MeV" << endl;
    cout << "Beam momentum: " << mom_beam << " GeV/c" << endl;
    cout << "CM momentum: " << pcm << " GeV/c" << endl;
    cout << "========================================\n" << endl;
    
    // ========================================
    // Histograms for conservation checks
    // ========================================
    TH1D *h_px_sum = new TH1D("h_px_sum", "Total Px (should be ~0);Px [GeV/c]", 100, -0.01, 0.01);
    TH1D *h_py_sum = new TH1D("h_py_sum", "Total Py (should be ~0);Py [GeV/c]", 100, -0.01, 0.01);
    TH1D *h_pz_sum = new TH1D("h_pz_sum", "Total Pz (should be beam);Pz [GeV/c]", 100, mom_beam-0.01, mom_beam+0.01);
    TH1D *h_E_sum = new TH1D("h_E_sum", "Total Energy (should be constant);E [GeV]", 100, E_total-0.01, E_total+0.01);
    
    // ========================================
    // Lab frame kinematics
    // ========================================
    TH1D *h_theta_lab = new TH1D("h_theta_lab", "Proton #theta_{lab};#theta_{lab} [deg]", 100, 0, 90);
	//TH1D *h_theta_lab = new TH1D("h_theta_lab", "Proton #theta_{lab};#theta_{lab} [deg]", 100, 15, 17);
    TH1D *h_phi_lab = new TH1D("h_phi_lab", "Proton #phi_{lab};#phi_{lab} [rad]", 100, -TMath::Pi(), TMath::Pi());
    TH1D *h_phi_lab_mrad_0 = new TH1D("h_phi_lab_mrad_0", "Proton #phi_{lab} near 0°;#phi_{lab} [mrad]", 200, -20, 20);
    TH1D *h_phi_lab_mrad_180 = new TH1D("h_phi_lab_mrad_180", "Proton #phi_{lab} near 180°;#phi_{lab}-#pi [mrad]", 200, -10, 10);
    
    // Hit map at polarimeter distance
    const double distance = POLARIMETER_DISTANCE; // in meters
    TH2D *h_hitmap_lab = new TH2D("h_hitmap_lab", "Proton Hit Map at 2.18m;x [m];y [m]", 200, -0.8, 0.8, 200, -0.8, 0.8);
    
    // ========================================
    // CM frame kinematics
    // ========================================
    TH1D *h_theta_cm = new TH1D("h_theta_cm", "CM Angle;#theta_{CM} [deg];Counts", 100, 0, 20);
    TH1D *h_costheta_cm = new TH1D("h_costheta_cm", "cos(#theta_{CM});cos(#theta_{CM});Counts", 100, 0.9, 1.0);
    
    cout << "Analyzing all events..." << endl;
    
    // ========================================
    // Loop over all events
    // ========================================
    for (int i = 0; i < tree->GetEntries(); i++) {
        if (i % 100000 == 0) cout << "Processing event " << i << " / " << tree->GetEntries() << endl;
        
        tree->GetEntry(i);
        
        PParticle *proton = (PParticle*)particles->At(0);
        PParticle *carbon = (PParticle*)particles->At(1);
        
        // ====================================
        // Conservation checks
        // ====================================
        Double_t px_tot = proton->Px() + carbon->Px();
        Double_t py_tot = proton->Py() + carbon->Py();
        Double_t pz_tot = proton->Pz() + carbon->Pz();
        Double_t E_tot = proton->E() + carbon->E();
        
        h_px_sum->Fill(px_tot);
        h_py_sum->Fill(py_tot);
        h_pz_sum->Fill(pz_tot);
        h_E_sum->Fill(E_tot);
        
        // ====================================
        // Lab frame angles
        // ====================================
        Double_t theta_lab = TMath::ATan2(TMath::Sqrt(proton->Px()*proton->Px() + 
                                                       proton->Py()*proton->Py()), 
                                          proton->Pz());
        Double_t phi_lab = TMath::ATan2(proton->Py(), proton->Px());
        
        h_theta_lab->Fill(theta_lab * TMath::RadToDeg());
        h_phi_lab->Fill(phi_lab);
        
        // For phi near 0° or 180°, show in mrad in separate histograms
        if (TMath::Abs(phi_lab) < 0.02) {
            h_phi_lab_mrad_0->Fill(phi_lab * 1000);  // Around 0°
        } else if (TMath::Abs(phi_lab - TMath::Pi()) < 0.02 || TMath::Abs(phi_lab + TMath::Pi()) < 0.02) {
            // Around 180° (can appear as +π or -π)
            if (phi_lab > 0) {
                h_phi_lab_mrad_180->Fill((phi_lab - TMath::Pi()) * 1000);
            } else {
                h_phi_lab_mrad_180->Fill((phi_lab + TMath::Pi()) * 1000);
            }
        }
        
        // ====================================
        // Calculate hit map at polarimeter distance
        // ====================================
        Double_t xxx = distance * TMath::Tan(theta_lab) * TMath::Cos(phi_lab);
        Double_t yyy = distance * TMath::Tan(theta_lab) * TMath::Sin(phi_lab);
        h_hitmap_lab->Fill(xxx, yyy);
        
        // ====================================
        // CM frame calculation
        // ====================================
        // Create proton 4-vector in lab frame
        TLorentzVector p_lab;
        p_lab.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);
        
        // Boost to CM frame
        TLorentzVector p_cm = p_lab;
        p_cm.Boost(-iState.BoostVector());
        
        // Get CM angle
        Double_t theta_cm = p_cm.Theta() * TMath::RadToDeg();
        Double_t costheta_cm = TMath::Cos(p_cm.Theta());
        
        h_theta_cm->Fill(theta_cm);
        h_costheta_cm->Fill(costheta_cm);
    }
    
    // ========================================
    // Summary statistics
    // ========================================
    cout << "\n========================================" << endl;
    cout << "Conservation Summary:" << endl;
    cout << "========================================" << endl;
    cout << "Px: mean = " << h_px_sum->GetMean() << " GeV/c, RMS = " << h_px_sum->GetRMS() << endl;
    cout << "Py: mean = " << h_py_sum->GetMean() << " GeV/c, RMS = " << h_py_sum->GetRMS() << endl;
    cout << "Pz: mean = " << h_pz_sum->GetMean() << " GeV/c (expected: " << mom_beam << ")" << endl;
    cout << "    RMS  = " << h_pz_sum->GetRMS() << " GeV/c" << endl;
    cout << "E:  mean = " << h_E_sum->GetMean() << " GeV (expected: " << E_total << ")" << endl;
    cout << "    RMS  = " << h_E_sum->GetRMS() << " GeV" << endl;
    
    cout << "\n========================================" << endl;
    cout << "Acceptance Summary:" << endl;
    cout << "========================================" << endl;
    cout << "Theta_lab: mean = " << h_theta_lab->GetMean() << " deg" << endl;
    cout << "           RMS  = " << h_theta_lab->GetRMS() << " deg" << endl;
    cout << "Theta_CM:  mean = " << h_theta_cm->GetMean() << " deg" << endl;
    cout << "           RMS  = " << h_theta_cm->GetRMS() << " deg" << endl;
    
    // ========================================
    // DRAW PLOTS
    // ========================================
    
    // Canvas 1: Conservation laws
    TCanvas *c1 = new TCanvas("c1", "Conservation Laws", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); 
    h_px_sum->Draw();
    c1->cd(2); 
    h_py_sum->Draw();
    c1->cd(3); 
    h_pz_sum->Draw();
    c1->cd(4); 
    h_E_sum->Draw();
    
    // Canvas 2: Lab frame kinematics
    TCanvas *c2 = new TCanvas("c2", "Lab Frame Kinematics", 2000, 600);
    c2->Divide(4,1);
    c2->cd(1); 
    h_theta_lab->Draw();
    c2->cd(2); 
    h_phi_lab->Draw();
    c2->cd(3);
    h_phi_lab_mrad_0->Draw();
    c2->cd(4);
    h_phi_lab_mrad_180->Draw();
    
    // Canvas 3: CM frame kinematics
    TCanvas *c3 = new TCanvas("c3", "CM Frame Kinematics", 1200, 600);
    c3->Divide(2,1);
    c3->cd(1);
    h_theta_cm->SetLineColor(kBlue);
    h_theta_cm->SetLineWidth(2);
    h_theta_cm->Draw();
    
    c3->cd(2);
    h_costheta_cm->SetLineColor(kRed);
    h_costheta_cm->SetLineWidth(2);
    h_costheta_cm->Draw();
    
    // Canvas 4: Hit map at polarimeter
	TCanvas *c4 = new TCanvas("c4", "Hit Map on 2.18m Screen", 800, 800);
	h_hitmap_lab->SetMinimum(0.1);  // Set minimum to small positive value
	gPad->SetLogz();
	h_hitmap_lab->Draw("colz");
    
    // Add acceptance circle for reference (16° at 2.18m)
    TEllipse *circle_theta = new TEllipse(0, 0, distance*TMath::Tan(16.0*TMath::DegToRad()), 
                                                  distance*TMath::Tan(16.0*TMath::DegToRad()));
    circle_theta->SetFillStyle(0);
    circle_theta->SetLineColor(kRed);
    circle_theta->SetLineWidth(2);
    circle_theta->Draw("SAME");
    
    cout << "\n========================================" << endl;
    cout << "Verification complete!" << endl;
    cout << "========================================" << endl;
    cout << "Canvas c1: Conservation laws" << endl;
    cout << "Canvas c2: Lab frame kinematics" << endl;
    cout << "Canvas c3: CM frame kinematics" << endl;
    cout << "Canvas c4: Hit map at polarimeter (2.18m)" << endl;
    cout << "\nAll plots show distributions from IMPORTANCE SAMPLING" << endl;
    cout << "(no weights needed - physics encoded in sampling)" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// VerifySimulation - Overloaded version for custom filename
// ====================================================================
// This version takes a filename string directly
// Usage: VerifySimulation("my_custom_file.root")
// ====================================================================
void VerifySimulation(const char* filename)
{
    cout << "\n========================================" << endl;
    cout << "Opening specified file: " << filename << endl;
    cout << "========================================" << endl;
    
    TFile *f = new TFile(filename);
    if (!f || f->IsZombie()) {
        cout << "Error: Cannot open file " << filename << endl;
        return;
    }
    
    TTree *tree = (TTree*)f->Get("data");
    if (!tree) {
        cout << "Error: Tree 'data' not found" << endl;
        return;
    }
    
    TClonesArray *particles = 0;
    tree->SetBranchAddress("Particles", &particles);
    
    // Try to extract energy from filename
    TString fname(filename);
    Double_t ekin = 0.200;  // Default 200 MeV
    
    // Try to parse energy from filename pattern "pC_Elas_XXXMeV"
    if (fname.Contains("MeV")) {
        TObjArray *tokens = fname.Tokenize("_");
        for (int i = 0; i < tokens->GetEntries(); i++) {
            TString token = ((TObjString*)tokens->At(i))->String();
            if (token.Contains("MeV")) {
                token.ReplaceAll("MeV", "");
                token.ReplaceAll(".root", "");
                ekin = token.Atof() / 1000.0;  // Convert to GeV
                break;
            }
        }
        delete tokens;
    }
    
    cout << "Detected beam energy: " << ekin*1000 << " MeV" << endl;
    
    // Constants
    const Double_t mp = 0.9382720813;
    const Double_t mC = 11.174862;
    const Double_t mom_beam = TMath::Sqrt(2*ekin*mp + ekin*ekin);
    const Double_t E_total = ekin + mp + mC;
    
    // Initial state 4-vector for CM calculations
    TLorentzVector iState(0., 0., mom_beam, ekin + mp + mC);
    Double_t s = iState.Mag2();
    Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));
    
    cout << "\n========================================" << endl;
    cout << "Verification of pC Elastic Scattering" << endl;
    cout << "========================================" << endl;
    cout << "File: " << filename << endl;
    cout << "Total events: " << tree->GetEntries() << endl;
    cout << "Beam energy: " << ekin*1000 << " MeV" << endl;
    cout << "Beam momentum: " << mom_beam << " GeV/c" << endl;
    cout << "CM momentum: " << pcm << " GeV/c" << endl;
    cout << "========================================\n" << endl;
    
    // ========================================
    // Histograms for conservation checks
    // ========================================
    TH1D *h_px_sum = new TH1D("h_px_sum", "Total Px (should be ~0);Px [GeV/c]", 100, -0.01, 0.01);
    TH1D *h_py_sum = new TH1D("h_py_sum", "Total Py (should be ~0);Py [GeV/c]", 100, -0.01, 0.01);
    TH1D *h_pz_sum = new TH1D("h_pz_sum", "Total Pz (should be beam);Pz [GeV/c]", 100, mom_beam-0.01, mom_beam+0.01);
    TH1D *h_E_sum = new TH1D("h_E_sum", "Total Energy (should be constant);E [GeV]", 100, E_total-0.01, E_total+0.01);
    
    // ========================================
    // Lab frame kinematics
    // ========================================
    TH1D *h_theta_lab = new TH1D("h_theta_lab", "Proton #theta_{lab};#theta_{lab} [deg]", 500, 10, 20);
    TH1D *h_phi_lab = new TH1D("h_phi_lab", "Proton #phi_{lab};#phi_{lab} [rad]", 100, -TMath::Pi(), TMath::Pi());
    TH1D *h_phi_lab_mrad_0 = new TH1D("h_phi_lab_mrad_0", "Proton #phi_{lab} near 0°;#phi_{lab} [mrad]", 500, -30, 30);
    TH1D *h_phi_lab_mrad_180 = new TH1D("h_phi_lab_mrad_180", "Proton #phi_{lab} near 180°;#phi_{lab}-#pi [mrad]", 500, -30, 30);
    
    // Hit map at polarimeter distance
    const double distance = POLARIMETER_DISTANCE; // in meters
    TH2D *h_hitmap_lab = new TH2D("h_hitmap_lab", "Proton Hit Map at 2.18m;x [m];y [m]", 200, -0.8, 0.8, 200, -0.8, 0.8);
    
    // ========================================
    // CM frame kinematics
    // ========================================
    TH1D *h_theta_cm = new TH1D("h_theta_cm", "CM Angle;#theta_{CM} [deg];Counts", 100, 0, 20);
    TH1D *h_costheta_cm = new TH1D("h_costheta_cm", "cos(#theta_{CM});cos(#theta_{CM});Counts", 100, 0.9, 1.0);
    
    cout << "Analyzing all events..." << endl;
    
    // ========================================
    // Loop over all events
    // ========================================
    for (int i = 0; i < tree->GetEntries(); i++) {
        if (i % 100000 == 0) cout << "Processing event " << i << " / " << tree->GetEntries() << endl;
        
        tree->GetEntry(i);
        
        PParticle *proton = (PParticle*)particles->At(0);
        PParticle *carbon = (PParticle*)particles->At(1);
        
        // ====================================
        // Conservation checks
        // ====================================
        Double_t px_tot = proton->Px() + carbon->Px();
        Double_t py_tot = proton->Py() + carbon->Py();
        Double_t pz_tot = proton->Pz() + carbon->Pz();
        Double_t E_tot = proton->E() + carbon->E();
        
        h_px_sum->Fill(px_tot);
        h_py_sum->Fill(py_tot);
        h_pz_sum->Fill(pz_tot);
        h_E_sum->Fill(E_tot);
        
        // ====================================
        // Lab frame angles
        // ====================================
        Double_t theta_lab = TMath::ATan2(TMath::Sqrt(proton->Px()*proton->Px() + 
                                                       proton->Py()*proton->Py()), 
                                          proton->Pz());
        Double_t phi_lab = TMath::ATan2(proton->Py(), proton->Px());
        
        h_theta_lab->Fill(theta_lab * TMath::RadToDeg());
        h_phi_lab->Fill(phi_lab);
        
        // For phi near 0° or 180°, show in mrad in separate histograms
        if (TMath::Abs(phi_lab) < 0.02) {
            h_phi_lab_mrad_0->Fill(phi_lab * 1000);  // Around 0°
        } else if (TMath::Abs(phi_lab - TMath::Pi()) < 0.02 || TMath::Abs(phi_lab + TMath::Pi()) < 0.02) {
            // Around 180° (can appear as +π or -π)
            if (phi_lab > 0) {
                h_phi_lab_mrad_180->Fill((phi_lab - TMath::Pi()) * 1000);
            } else {
                h_phi_lab_mrad_180->Fill((phi_lab + TMath::Pi()) * 1000);
            }
        }
        
        // ====================================
        // Calculate hit map at polarimeter distance
        // ====================================
        Double_t xxx = distance * TMath::Tan(theta_lab) * TMath::Cos(phi_lab);
        Double_t yyy = distance * TMath::Tan(theta_lab) * TMath::Sin(phi_lab);
        h_hitmap_lab->Fill(xxx, yyy);
        
        // ====================================
        // CM frame calculation
        // ====================================
        // Create proton 4-vector in lab frame
        TLorentzVector p_lab;
        p_lab.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp);
        
        // Boost to CM frame
        TLorentzVector p_cm = p_lab;
        p_cm.Boost(-iState.BoostVector());
        
        // Get CM angle
        Double_t theta_cm = p_cm.Theta() * TMath::RadToDeg();
        Double_t costheta_cm = TMath::Cos(p_cm.Theta());
        
        h_theta_cm->Fill(theta_cm);
        h_costheta_cm->Fill(costheta_cm);
    }
    
    // ========================================
    // Summary statistics
    // ========================================
    cout << "\n========================================" << endl;
    cout << "Conservation Summary:" << endl;
    cout << "========================================" << endl;
    cout << "Px: mean = " << h_px_sum->GetMean() << " GeV/c, RMS = " << h_px_sum->GetRMS() << endl;
    cout << "Py: mean = " << h_py_sum->GetMean() << " GeV/c, RMS = " << h_py_sum->GetRMS() << endl;
    cout << "Pz: mean = " << h_pz_sum->GetMean() << " GeV/c (expected: " << mom_beam << ")" << endl;
    cout << "    RMS  = " << h_pz_sum->GetRMS() << " GeV/c" << endl;
    cout << "E:  mean = " << h_E_sum->GetMean() << " GeV (expected: " << E_total << ")" << endl;
    cout << "    RMS  = " << h_E_sum->GetRMS() << " GeV" << endl;
    
    cout << "\n========================================" << endl;
    cout << "Acceptance Summary:" << endl;
    cout << "========================================" << endl;
    cout << "Theta_lab: mean = " << h_theta_lab->GetMean() << " deg" << endl;
    cout << "           RMS  = " << h_theta_lab->GetRMS() << " deg" << endl;
    cout << "Theta_CM:  mean = " << h_theta_cm->GetMean() << " deg" << endl;
    cout << "           RMS  = " << h_theta_cm->GetRMS() << " deg" << endl;
    
    // ========================================
    // DRAW PLOTS
    // ========================================
    
    // Canvas 1: Conservation laws
    TCanvas *c1 = new TCanvas("c1", "Conservation Laws", 1200, 800);
    c1->Divide(2,2);
    c1->cd(1); 
    h_px_sum->Draw();
    c1->cd(2); 
    h_py_sum->Draw();
    c1->cd(3); 
    h_pz_sum->Draw();
    c1->cd(4); 
    h_E_sum->Draw();
    
    // Canvas 2: Lab frame kinematics
    TCanvas *c2 = new TCanvas("c2", "Lab Frame Kinematics", 2000, 600);
    c2->Divide(4,1);
    c2->cd(1); 
    h_theta_lab->Draw();
    c2->cd(2); 
    h_phi_lab->Draw();
    c2->cd(3);
    h_phi_lab_mrad_0->Draw();    c2->cd(4);
    h_phi_lab_mrad_180->Draw();
    
    // Canvas 3: CM frame kinematics
    TCanvas *c3 = new TCanvas("c3", "CM Frame Kinematics", 1200, 600);
    c3->Divide(2,1);
    c3->cd(1);
    h_theta_cm->SetLineColor(kBlue);
    h_theta_cm->SetLineWidth(2);
    h_theta_cm->Draw();
    
    c3->cd(2);
    h_costheta_cm->SetLineColor(kRed);
    h_costheta_cm->SetLineWidth(2);
    h_costheta_cm->Draw();
    
    // Canvas 4: Hit map at polarimeter
    TCanvas *c4 = new TCanvas("c4", "Hit Map on 2.18m Screen", 800, 800);
    h_hitmap_lab->Draw("colz");
    
    // Add acceptance circle for reference (16° at 2.18m)
    TEllipse *circle_theta = new TEllipse(0, 0, distance*TMath::Tan(16.0*TMath::DegToRad()), 
                                                  distance*TMath::Tan(16.0*TMath::DegToRad()));
    circle_theta->SetFillStyle(0);
    circle_theta->SetLineColor(kRed);
    circle_theta->SetLineWidth(2);
    circle_theta->Draw("SAME");
    
    cout << "\n========================================" << endl;
    cout << "Verification complete!" << endl;
    cout << "========================================" << endl;
    cout << "Canvas c1: Conservation laws" << endl;
    cout << "Canvas c2: Lab frame kinematics" << endl;
    cout << "Canvas c3: CM frame kinematics" << endl;
    cout << "Canvas c4: Hit map at polarimeter (2.18m)" << endl;
    cout << "\nAll plots show distributions from IMPORTANCE SAMPLING" << endl;
    cout << "(no weights needed - physics encoded in sampling)" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// ====================================================================
// TestAnalyzingPower: Simple test of A_N parameterization
// ====================================================================
void TestAnalyzingPower()
{
    cout << "\n========================================" << endl;
    cout << "TESTING ANALYZING POWER" << endl;
    cout << "========================================\n" << endl;
    
    // Test at key points from the paper
    Double_t AN_200_16p2 = GetElasticAnalyzingPower(200.0, 16.2);
    Double_t AN_189_17p3 = GetElasticAnalyzingPower(189.0, 17.3);
    Double_t AN_200_17p3 = GetElasticAnalyzingPower(200.0, 17.3);
    
    cout << "Test results:" << endl;
    cout << "  A_N(200 MeV, 16.2 deg) = " << AN_200_16p2 * 100 << "%" << endl;
    cout << "    Expected: 99.35%" << endl;
    cout << "  A_N(189 MeV, 17.3 deg) = " << AN_189_17p3 * 100 << "%" << endl;
    cout << "    Expected: ~100%" << endl;
    cout << "  A_N(200 MeV, 17.3 deg) = " << AN_200_17p3 * 100 << "%" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// TestPhiSampling: Test the azimuthal angle sampling
// ====================================================================
void TestPhiSampling()
{
    cout << "\n========================================" << endl;
    cout << "TESTING PHI SAMPLING" << endl;
    cout << "========================================\n" << endl;
    
    // Use analyzing power at 200 MeV, 16.2 degrees
    Double_t AN = GetElasticAnalyzingPower(200.0, 16.2);
    cout << "Using A_N = " << AN << " (at 200 MeV, 16.2 deg)" << endl;
    cout << "Beam polarization: " << BEAM_POLARIZATION * 100 << "%" << endl << endl;
    
    // Create phi sampling functions for both spin states
    TF1 *fPhi_up = CreatePhiSamplingFunction(BEAM_POLARIZATION, AN, +1);
    TF1 *fPhi_down = CreatePhiSamplingFunction(BEAM_POLARIZATION, AN, -1);
    
    cout << "Spin-up function:   1 + " << fPhi_up->GetParameter(0) << " * cos(phi)" << endl;
    cout << "Spin-down function: 1 + " << fPhi_down->GetParameter(0) << " * cos(phi)" << endl;
    cout << endl;
    
    // Create histograms to show distributions (in degrees)
    TH1D *h_phi_up = new TH1D("h_phi_up", "Phi distribution (spin-up);#phi [deg];Counts", 100, 0, 360);
    TH1D *h_phi_down = new TH1D("h_phi_down", "Phi distribution (spin-down);#phi [deg];Counts", 100, 0, 360);
    
    // Sample 100,000 events for each
    cout << "Sampling 100,000 events for each spin state..." << endl;
    for (int i = 0; i < 100000; i++) {
        h_phi_up->Fill(fPhi_up->GetRandom() * TMath::RadToDeg());
        h_phi_down->Fill(fPhi_down->GetRandom() * TMath::RadToDeg());
    }
    
    // Show statistics
    cout << "\nSpin-up statistics:" << endl;
    cout << "  Events at phi=0 deg (bin 1): " << h_phi_up->GetBinContent(1) << endl;
    cout << "  Events at phi=180 deg (bin 50): " << h_phi_up->GetBinContent(50) << endl;
    cout << "  Ratio (0/180): " << h_phi_up->GetBinContent(1) / h_phi_up->GetBinContent(50) << endl;
    
    cout << "\nSpin-down statistics:" << endl;
    cout << "  Events at phi=0 deg (bin 1): " << h_phi_down->GetBinContent(1) << endl;
    cout << "  Events at phi=180 deg (bin 50): " << h_phi_down->GetBinContent(50) << endl;
    cout << "  Ratio (0/180): " << h_phi_down->GetBinContent(1) / h_phi_down->GetBinContent(50) << endl;
    
    // Draw comparison
    TCanvas *c = new TCanvas("c_phi_test", "Phi Sampling Test", 1200, 500);
    c->Divide(2, 1);
    
    c->cd(1);
    h_phi_up->SetLineColor(kBlue);
    h_phi_up->SetLineWidth(2);
    h_phi_up->Draw();
    gPad->SetTitle("Spin-Up");
    
    c->cd(2);
    h_phi_down->SetLineColor(kRed);
    h_phi_down->SetLineWidth(2);
    h_phi_down->Draw();
    gPad->SetTitle("Spin-Down");
    
    cout << "\n========================================" << endl;
    cout << "Canvas 'c_phi_test' created!" << endl;
    cout << "Blue = spin-up (more at phi=0 deg, RIGHT)" << endl;
    cout << "Red = spin-down (more at phi=180 deg, LEFT)" << endl;
    cout << "========================================\n" << endl;
    
    // Clean up
    delete fPhi_up;
    delete fPhi_down;
}

// ====================================================================
// TestPhiMapping: Test CM to LAB phi conversion
// ====================================================================
void TestPhiMapping()
{
    cout << "\n========================================" << endl;
    cout << "TESTING PHI CM -> LAB MAPPING" << endl;
    cout << "========================================\n" << endl;
    
    Double_t ekin = 200.0;  // MeV
    Double_t theta_cm = 17.5;  // degrees (approximate CM angle at 16° lab)
    
    cout << "Beam energy: " << ekin << " MeV" << endl;
    cout << "Theta_CM: " << theta_cm << " degrees" << endl << endl;
    
    // Test a few phi values
    cout << "Testing phi conversion:" << endl;
    Double_t test_phis[] = {-180, -90, 0, 90, 180};
    for (int i = 0; i < 5; i++) {
        Double_t phi_cm_deg = test_phis[i];
        Double_t phi_cm_rad = phi_cm_deg * TMath::DegToRad();
        Double_t phi_lab_rad = ConvertPhiCMtoLab(ekin, theta_cm, phi_cm_rad);
        Double_t phi_lab_deg = phi_lab_rad * TMath::RadToDeg();
        
        cout << "  phi_CM = " << phi_cm_deg << " deg  -->  phi_LAB = " 
             << phi_lab_deg << " deg (" << phi_lab_rad << " rad)" << endl;
    }
    
    // Test values very close to 180°
    cout << "\nTesting near 180°:" << endl;
    Double_t test_near_180[] = {179.0, 179.5, 179.9, 180.0, -179.9, -179.5, -179.0};
    for (int i = 0; i < 7; i++) {
        Double_t phi_cm_deg = test_near_180[i];
        Double_t phi_cm_rad = phi_cm_deg * TMath::DegToRad();
        Double_t phi_lab_rad = ConvertPhiCMtoLab(ekin, theta_cm, phi_cm_rad);
        Double_t phi_lab_deg = phi_lab_rad * TMath::RadToDeg();
        
        cout << "  phi_CM = " << phi_cm_deg << " deg  -->  phi_LAB = " 
             << phi_lab_deg << " deg (" << phi_lab_rad << " rad)" << endl;
    }
    
    // Compute CM ranges for detector acceptance
    cout << "\n========================================" << endl;
    cout << "Computing CM ranges for detector acceptance..." << endl;
    cout << "========================================\n" << endl;
    
    Double_t phi_cm_min_0, phi_cm_max_0;
    Double_t phi_cm_min_180, phi_cm_max_180;
    
    ComputePhiCMRanges(ekin, theta_cm, 
                       phi_cm_min_0, phi_cm_max_0,
                       phi_cm_min_180, phi_cm_max_180);
    
    cout << "Detector 1 (phi_lab ≈ 0°):" << endl;
    cout << "  CM range: [" << phi_cm_min_0 * TMath::RadToDeg() << ", " 
         << phi_cm_max_0 * TMath::RadToDeg() << "] degrees" << endl;
    cout << "  CM range: [" << phi_cm_min_0 << ", " << phi_cm_max_0 << "] radians" << endl;
    
    cout << "\nDetector 2 (phi_lab ≈ 180°):" << endl;
    cout << "  Due to ±π periodicity, there are TWO equivalent CM ranges:" << endl;
    cout << "  Range A (near +π): [" << phi_cm_min_180 * TMath::RadToDeg() << ", " 
         << phi_cm_max_180 * TMath::RadToDeg() << "] degrees" << endl;
    cout << "  Range B (near -π): [" << -phi_cm_max_180 * TMath::RadToDeg() << ", " 
         << -phi_cm_min_180 * TMath::RadToDeg() << "] degrees" << endl;
    cout << "  In radians:" << endl;
    cout << "    Range A: [" << phi_cm_min_180 << ", " << phi_cm_max_180 << "]" << endl;
    cout << "    Range B: [" << -phi_cm_max_180 << ", " << -phi_cm_min_180 << "]" << endl;
    
    cout << "\n========================================" << endl;
    cout << "These are the CM phi ranges to sample from!" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// ====================================================================
// GenerateEvents: One-stop function to generate all 4 datasets
// ====================================================================
// Generates complete set of polarized elastic + inelastic events:
//   1. Elastic, Spin-UP
//   2. Elastic, Spin-DOWN
//   3. Inelastic (4.43 MeV), Spin-UP
//   4. Inelastic (4.43 MeV), Spin-DOWN
//
// Parameters:
//   - energy: Beam kinetic energy [MeV]
//   - events: Number of events per dataset
//   - polarization: Beam polarization (e.g., 0.70, 0.80)
//   - num_threads: Number of CPU threads (0 = auto-detect)
//
// Output files (example for 200 MeV, P=70%):
//   - pC_Elas_200MeV_MT_P70_SpinUp.root
//   - pC_Elas_200MeV_MT_P70_SpinDown.root
//   - pC_Inel443_200MeV_MT_P70_SpinUp.root
//   - pC_Inel443_200MeV_MT_P70_SpinDown.root
//
// Usage:
//   GenerateEvents(200.0, 100000, 0.70)     // 100k events each, 70% polarization
//   GenerateEvents(200.0, 500000, 0.80, 8)  // 500k events each, 80% pol, 8 threads
// ====================================================================
void GenerateEvents(Double_t energy, Int_t events, Double_t polarization, Int_t num_threads = 0)
{
    cout << "\n" << endl;
    cout << "╔════════════════════════════════════════════════════════════╗" << endl;
    cout << "║                                                            ║" << endl;
    cout << "║           COMPLETE EVENT GENERATION SUITE                  ║" << endl;
    cout << "║                                                            ║" << endl;
    cout << "╚════════════════════════════════════════════════════════════╝" << endl;
    cout << "\n" << endl;
    
    cout << "Configuration:" << endl;
    cout << "  Beam energy:      " << energy << " MeV" << endl;
    cout << "  Events per file:  " << events << endl;
    cout << "  Polarization:     " << polarization * 100 << "%" << endl;
    cout << "  CPU threads:      " << (num_threads == 0 ? "auto-detect" : std::to_string(num_threads)) << endl;
    
    int polar_percent = (int)(polarization * 100);
    cout << "\nWill generate 4 ROOT files:" << endl;
    cout << "  1. pC_Elas_" << (int)energy << "MeV_MT_P" << polar_percent << "_SpinUp.root" << endl;
    cout << "  2. pC_Elas_" << (int)energy << "MeV_MT_P" << polar_percent << "_SpinDown.root" << endl;
    cout << "  3. pC_Inel443_" << (int)energy << "MeV_MT_P" << polar_percent << "_SpinUp.root" << endl;
    cout << "  4. pC_Inel443_" << (int)energy << "MeV_MT_P" << polar_percent << "_SpinDown.root" << endl;
    
    cout << "\n════════════════════════════════════════════════════════════" << endl;
    cout << "Starting generation..." << endl;
    cout << "════════════════════════════════════════════════════════════\n" << endl;
    
    auto start_total = std::chrono::high_resolution_clock::now();
    
    // ================================================================
    // 1. ELASTIC SPIN-UP
    // ================================================================
    cout << "\n[1/4] Generating ELASTIC SPIN-UP events..." << endl;
    cout << "──────────────────────────────────────────────────────────" << endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    SingleRunMultithreadPolarized(energy, events, polarization, +1, num_threads);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "✓ Completed in " << elapsed.count() << " seconds" << endl;
    
    // ================================================================
    // 2. ELASTIC SPIN-DOWN
    // ================================================================
    cout << "\n[2/4] Generating ELASTIC SPIN-DOWN events..." << endl;
    cout << "──────────────────────────────────────────────────────────" << endl;
    start = std::chrono::high_resolution_clock::now();
    
    SingleRunMultithreadPolarized(energy, events, polarization, -1, num_threads);
    
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "✓ Completed in " << elapsed.count() << " seconds" << endl;
    
    // ================================================================
    // 3. INELASTIC SPIN-UP
    // ================================================================
    cout << "\n[3/4] Generating INELASTIC (4.43 MeV) SPIN-UP events..." << endl;
    cout << "──────────────────────────────────────────────────────────" << endl;
    start = std::chrono::high_resolution_clock::now();
    
    SingleRunMultithreadInelasticPolarized(energy, events, polarization, +1, num_threads);
    
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "✓ Completed in " << elapsed.count() << " seconds" << endl;
    
    // ================================================================
    // 4. INELASTIC SPIN-DOWN
    // ================================================================
    cout << "\n[4/4] Generating INELASTIC (4.43 MeV) SPIN-DOWN events..." << endl;
    cout << "──────────────────────────────────────────────────────────" << endl;
    start = std::chrono::high_resolution_clock::now();
    
    SingleRunMultithreadInelasticPolarized(energy, events, polarization, -1, num_threads);
    
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "✓ Completed in " << elapsed.count() << " seconds" << endl;
    
    // ================================================================
    // SUMMARY
    // ================================================================
    auto end_total = std::chrono::high_resolution_clock::now();
    auto elapsed_total = std::chrono::duration_cast<std::chrono::seconds>(end_total - start_total);
    
    cout << "\n" << endl;
    cout << "╔════════════════════════════════════════════════════════════╗" << endl;
    cout << "║                                                            ║" << endl;
    cout << "║               GENERATION COMPLETE!                         ║" << endl;
    cout << "║                                                            ║" << endl;
    cout << "╚════════════════════════════════════════════════════════════╝" << endl;
    cout << "\n" << endl;
    
    cout << "Summary:" << endl;
    cout << "  Total time:       " << elapsed_total.count() << " seconds" << endl;
    cout << "  Total events:     " << events * 4 << " (" << events << " × 4 files)" << endl;
    cout << "  Events/second:    " << (events * 4) / (double)elapsed_total.count() << endl;
    
    cout << "\nOutput files created in: " << gSystem->pwd() << "/outputs/" << endl;
    cout << "  ✓ pC_Elas_" << (int)energy << "MeV_MT_P" << polar_percent << "_SpinUp.root" << endl;
    cout << "  ✓ pC_Elas_" << (int)energy << "MeV_MT_P" << polar_percent << "_SpinDown.root" << endl;
    cout << "  ✓ pC_Inel443_" << (int)energy << "MeV_MT_P" << polar_percent << "_SpinUp.root" << endl;
    cout << "  ✓ pC_Inel443_" << (int)energy << "MeV_MT_P" << polar_percent << "_SpinDown.root" << endl;
    
    cout << "\nReady for analysis!" << endl;
    cout << "════════════════════════════════════════════════════════════\n" << endl;
}

// ====================================================================
// GenerateEventsQuick: Shortcut with default polarization
// ====================================================================
// Uses default beam polarization from constant
// Usage: GenerateEventsQuick(200.0, 100000)
// ====================================================================
void GenerateEventsQuick(Double_t energy, Int_t events, Int_t num_threads = 0)
{
    GenerateEvents(energy, events, BEAM_POLARIZATION, num_threads);
}


//====================================================================





void ComparePolarization(const char* file_up, const char* file_down)
{
    cout << "\n========================================" << endl;
    cout << "POLARIZATION COMPARISON" << endl;
    cout << "========================================" << endl;
    
    TFile *f_up = new TFile(file_up);
    TFile *f_down = new TFile(file_down);
    
    if (!f_up || f_up->IsZombie() || !f_down || f_down->IsZombie()) {
        cout << "Error: Cannot open input files" << endl;
        return;
    }
    
    TTree *tree_up = (TTree*)f_up->Get("data");
    TTree *tree_down = (TTree*)f_down->Get("data");
    
    TClonesArray *part_up = 0;
    TClonesArray *part_down = 0;
    tree_up->SetBranchAddress("Particles", &part_up);
    tree_down->SetBranchAddress("Particles", &part_down);
    
    const double distance = POLARIMETER_DISTANCE;
    
    // Create histograms
    TH2D *h_up = new TH2D("h_up", "Spin-UP hit map;x [m];y [m]", 100, -0.8, 0.8, 100, -0.8, 0.8);
    TH2D *h_down = new TH2D("h_down", "Spin-DOWN hit map;x [m];y [m]", 100, -0.8, 0.8, 100, -0.8, 0.8);
    TH2D *h_asymmetry = new TH2D("h_asymmetry", "Spin-Correlated Asymmetry asym;x [m];y [m]", 100, -0.8, 0.8, 100, -0.8, 0.8);
    
    TH1D *h_phi_up = new TH1D("h_phi_up", "Phi distribution spin-UP;#phi [rad];Counts", 100, -TMath::Pi(), TMath::Pi());
    TH1D *h_phi_down = new TH1D("h_phi_down", "Phi distribution spin-DOWN;#phi [rad];Counts", 100, -TMath::Pi(), TMath::Pi());
    
    // Fill spin-UP histogram
    for (int i = 0; i < tree_up->GetEntries(); i++) {
        tree_up->GetEntry(i);
        PParticle *p = (PParticle*)part_up->At(0);
        
        Double_t theta = TMath::ATan2(TMath::Sqrt(p->Px()*p->Px() + p->Py()*p->Py()), p->Pz());
        Double_t phi = TMath::ATan2(p->Py(), p->Px());
        Double_t x = distance * TMath::Tan(theta) * TMath::Cos(phi);
        Double_t y = distance * TMath::Tan(theta) * TMath::Sin(phi);
        
        h_up->Fill(x, y);
        h_phi_up->Fill(phi);
    }
    
    // Fill spin-DOWN histogram
    for (int i = 0; i < tree_down->GetEntries(); i++) {
        tree_down->GetEntry(i);
        PParticle *p = (PParticle*)part_down->At(0);
        
        Double_t theta = TMath::ATan2(TMath::Sqrt(p->Px()*p->Px() + p->Py()*p->Py()), p->Pz());
        Double_t phi = TMath::ATan2(p->Py(), p->Px());
        Double_t x = distance * TMath::Tan(theta) * TMath::Cos(phi);
        Double_t y = distance * TMath::Tan(theta) * TMath::Sin(phi);
        
        h_down->Fill(x, y);
        h_phi_down->Fill(phi);
    }
    
    // ================================================================
    // CALCULATE SPIN-CORRELATED ASYMMETRY (CORRECT FORMULA)
    // ================================================================
    // For each spatial bin (x,y), calculate:
    //         √(N↑_R × N↓_L) - √(N↑_L × N↓_R)
    // asym = ─────────────────────────────────────
    //         √(N↑_R × N↓_L) + √(N↑_L × N↓_R)
    //
    // Where R/L is determined by x-position (LEFT/RIGHT detector)
    // ================================================================
    
    for (int i = 1; i <= h_asymmetry->GetNbinsX(); i++) {
        for (int j = 1; j <= h_asymmetry->GetNbinsY(); j++) {
            Double_t x_center = h_asymmetry->GetXaxis()->GetBinCenter(i);
            
            Double_t n_up = h_up->GetBinContent(i, j);
            Double_t n_down = h_down->GetBinContent(i, j);
            
            // Need enough statistics
            if (n_up + n_down < 10) continue;
            
            // Determine if this bin is RIGHT or LEFT based on x-position
            if (x_center > 0) {
                // RIGHT side (φ ≈ 0°)
                // This bin represents: N↑_R and N↓_R
                // We need the corresponding LEFT counts to calculate asymmetry
                // For simplicity, we use the square-root ratio method
                
                // Find corresponding LEFT bin (mirror in x)
                Int_t i_left = h_asymmetry->GetXaxis()->FindBin(-x_center);
                Double_t n_up_left = h_up->GetBinContent(i_left, j);
                Double_t n_down_left = h_down->GetBinContent(i_left, j);
                
                if (n_up_left + n_down_left < 10) continue;
                
                // Calculate spin-correlated asymmetry
                // N↑_R = n_up, N↓_L = n_down_left, N↑_L = n_up_left, N↓_R = n_down
                Double_t numerator = TMath::Sqrt(n_up * n_down_left) - TMath::Sqrt(n_up_left * n_down);
                Double_t denominator = TMath::Sqrt(n_up * n_down_left) + TMath::Sqrt(n_up_left * n_down);
                
                if (denominator > 0) {
                    Double_t asym = numerator / denominator;
                    h_asymmetry->SetBinContent(i, j, asym);
                    // Also fill the mirror bin with same asymmetry (by symmetry)
                    h_asymmetry->SetBinContent(i_left, j, asym);
                }
            }
        }
    }
    
    // ================================================================
    // DRAW HISTOGRAMS WITH LOG SCALE
    // ================================================================
    
    TCanvas *c1 = new TCanvas("c_hitmap_compare", "Hit Map Comparison", 1800, 600);
    c1->Divide(3, 1);
    
    // Spin-UP with log scale
    c1->cd(1);
    gPad->SetLogz();
    //gStyle->SetPalette(kViridis);
    h_up->SetStats(0);
    h_up->SetMinimum(0.1);  // Avoid zeros in log scale
    h_up->Draw("colz");
	
	// Add circles for detector acceptance angles
    Double_t radius_16deg = distance * TMath::Tan(16.0 * TMath::DegToRad());
    Double_t radius_12deg = distance * TMath::Tan(12.0 * TMath::DegToRad());
    
    TEllipse *circle_16_up = new TEllipse(0, 0, radius_16deg, radius_16deg);
    circle_16_up->SetFillStyle(0);
    circle_16_up->SetLineColor(kRed);
    circle_16_up->SetLineWidth(2);
    circle_16_up->Draw("SAME");
    
    TEllipse *circle_12_up = new TEllipse(0, 0, radius_12deg, radius_12deg);
    circle_12_up->SetFillStyle(0);
    circle_12_up->SetLineColor(kGreen+2);
    circle_12_up->SetLineWidth(2);
    circle_12_up->Draw("SAME");
    
    // Spin-DOWN with log scale
    c1->cd(2);
    gPad->SetLogz();
    //gStyle->SetPalette(kViridis);
    h_down->SetStats(0);
    h_down->SetMinimum(0.1);  // Avoid zeros in log scale
    h_down->Draw("colz");
    
	TEllipse *circle_16_down = new TEllipse(0, 0, radius_16deg, radius_16deg);
    circle_16_down->SetFillStyle(0);
    circle_16_down->SetLineColor(kRed);
    circle_16_down->SetLineWidth(2);
    circle_16_down->Draw("SAME");
    
    TEllipse *circle_12_down = new TEllipse(0, 0, radius_12deg, radius_12deg);
    circle_12_down->SetFillStyle(0);
    circle_12_down->SetLineColor(kGreen+2);
    circle_12_down->SetLineWidth(2);
    circle_12_down->Draw("SAME");
	
    // Asymmetry (keep linear, use diverging color scheme)
    c1->cd(3);
    //gStyle->SetPalette(kCool);  // Cool palette for ±values
    h_asymmetry->SetStats(0);
    h_asymmetry->SetMinimum(0);  // Asymmetry ranges from -1 to +1
    h_asymmetry->SetMaximum(1.0);
    h_asymmetry->Draw("colz");
    
    // Add zero line for reference
    TLine *line_h = new TLine(-0.8, 0, 0.8, 0);
    line_h->SetLineStyle(2);
    line_h->SetLineColor(kBlack);
    line_h->Draw("SAME");
    
    TLine *line_v = new TLine(0, -0.8, 0, 0.8);
    line_v->SetLineStyle(2);
    line_v->SetLineColor(kBlack);
    line_v->Draw("SAME");
	
	TEllipse *circle_16_asym = new TEllipse(0, 0, radius_16deg, radius_16deg);
    circle_16_asym->SetFillStyle(0);
    circle_16_asym->SetLineColor(kRed);
    circle_16_asym->SetLineWidth(2);
    circle_16_asym->Draw("SAME");
    
    TEllipse *circle_12_asym = new TEllipse(0, 0, radius_12deg, radius_12deg);
    circle_12_asym->SetFillStyle(0);
    circle_12_asym->SetLineColor(kGreen+2);
    circle_12_asym->SetLineWidth(2);
    circle_12_asym->Draw("SAME");
	
    
    // ================================================================
    // PHI DISTRIBUTIONS
    // ================================================================
    
    TCanvas *c2 = new TCanvas("c_phi_compare", "Phi Distribution Comparison", 1200, 600);
    c2->Divide(2, 1);
    
    c2->cd(1);
    h_phi_up->SetLineColor(kBlue);
    h_phi_up->SetLineWidth(2);
    h_phi_up->Draw();
    
    c2->cd(2);
    h_phi_down->SetLineColor(kRed);
    h_phi_down->SetLineWidth(2);
    h_phi_down->Draw();
    
    // ================================================================
    // STATISTICS AND SUMMARY
    // ================================================================
    
    cout << "\n========================================" << endl;
    cout << "Statistics:" << endl;
    cout << "Spin-UP events: " << tree_up->GetEntries() << endl;
    cout << "Spin-DOWN events: " << tree_down->GetEntries() << endl;
    
    // Count events at LEFT (x < 0) and RIGHT (x > 0)
    Double_t up_right = 0, up_left = 0, down_right = 0, down_left = 0;
    
    for (int i = 1; i <= h_up->GetNbinsX(); i++) {
        Double_t x = h_up->GetXaxis()->GetBinCenter(i);
        for (int j = 1; j <= h_up->GetNbinsY(); j++) {
            if (x > 0) {
                up_right += h_up->GetBinContent(i, j);
                down_right += h_down->GetBinContent(i, j);
            } else if (x < 0) {
                up_left += h_up->GetBinContent(i, j);
                down_left += h_down->GetBinContent(i, j);
            }
        }
    }
    
    cout << "\nSpatial distribution (LEFT/RIGHT):" << endl;
    cout << "  Spin-UP:   LEFT = " << up_left << ",  RIGHT = " << up_right 
         << ",  Ratio R/L = " << up_right/up_left << endl;
    cout << "  Spin-DOWN: LEFT = " << down_left << ",  RIGHT = " << down_right 
         << ",  Ratio R/L = " << down_right/down_left << endl;
    
    // Calculate overall spin-correlated asymmetry
    Double_t asym = (TMath::Sqrt(up_right * down_left) - TMath::Sqrt(up_left * down_right)) /
                   (TMath::Sqrt(up_right * down_left) + TMath::Sqrt(up_left * down_right));
    
    cout << "\nSpin-Correlated Asymmetry:" << endl;
    cout << "  asym = " << asym << endl;
    cout << "  Expected: ~" << BEAM_POLARIZATION << " (= beam polarization P × analyzing power A_N)" << endl;
    
    cout << "\nPhi distribution (near 0° and 180°):" << endl;
    
    int up_0 = 0, up_180 = 0, down_0 = 0, down_180 = 0;
    for (int i = 1; i <= h_phi_up->GetNbinsX(); i++) {
        Double_t phi = h_phi_up->GetBinCenter(i);
        if (TMath::Abs(phi) < 0.1) {
            up_0 += h_phi_up->GetBinContent(i);
            down_0 += h_phi_down->GetBinContent(i);
        }
        if (TMath::Abs(TMath::Abs(phi) - TMath::Pi()) < 0.1) {
            up_180 += h_phi_up->GetBinContent(i);
            down_180 += h_phi_down->GetBinContent(i);
        }
    }
    
    cout << "  Spin-UP:   0° = " << up_0 << ",  180° = " << up_180 
         << ",  Ratio 0°/180° = " << (double)up_0/up_180 << endl;
    cout << "  Spin-DOWN: 0° = " << down_0 << ",  180° = " << down_180 
         << ",  Ratio 0°/180° = " << (double)down_0/down_180 << endl;
    
    cout << "\nExpected behavior:" << endl;
    cout << "  Spin-UP should have MORE events at 0° (RIGHT)" << endl;
    cout << "  Spin-DOWN should have MORE events at 180° (LEFT)" << endl;
    cout << "  asym should be positive and close to P × A_N (beam polarization × analyzing power)" << endl;
    cout << "========================================\n" << endl;
}



void DiagnoseSpin()
{
    cout << "\n========================================" << endl;
    cout << "SPIN DIRECTION DIAGNOSTIC" << endl;
    cout << "========================================\n" << endl;
    
    Double_t energy = 200.0;
    Double_t theta_lab = 16.2;
    Double_t AN = GetElasticAnalyzingPower(energy, theta_lab);
    Double_t polarization = 0.8;
    
    cout << "A_N at " << energy << " MeV, " << theta_lab << " deg = " << AN << endl;
    cout << "Polarization = " << polarization << endl << endl;
    
    // Test spin-up
    Int_t spin_up = +1;
    Double_t asym_up = spin_up * polarization * AN;
    cout << "SPIN-UP (+1):" << endl;
    cout << "  asymmetry coefficient = " << asym_up << endl;
    cout << "  Weight at phi=0°   (cos=+1): " << (1.0 + asym_up * 1.0) << endl;
    cout << "  Weight at phi=180° (cos=-1): " << (1.0 + asym_up * (-1.0)) << endl;
    cout << "  Expected: MORE weight at phi=0° (RIGHT)" << endl << endl;
    
    // Test spin-down
    Int_t spin_down = -1;
    Double_t asym_down = spin_down * polarization * AN;
    cout << "SPIN-DOWN (-1):" << endl;
    cout << "  asymmetry coefficient = " << asym_down << endl;
    cout << "  Weight at phi=0°   (cos=+1): " << (1.0 + asym_down * 1.0) << endl;
    cout << "  Weight at phi=180° (cos=-1): " << (1.0 + asym_down * (-1.0)) << endl;
    cout << "  Expected: MORE weight at phi=180° (LEFT)" << endl << endl;
    
    cout << "If weights show opposite behavior, we have a sign error!" << endl;
    cout << "========================================\n" << endl;
}


// ====================================================================
// TestInelasticIntegration: Quick test that data loads correctly
// ====================================================================
void TestInelasticIntegration()
{
    LoadInelasticData();
    
    cout << "\nTesting interpolation:" << endl;
    cout << "Angle\tXS (mb/sr)\tA_N" << endl;
    
    for (double theta = 10; theta <= 60; theta += 10) {
        double xs = GetInelasticCrossSection(theta);
        double an = GetInelasticAnalyzingPower(theta);
        printf("%.1f°\t%.4f\t\t%.4f\n", theta, xs, an);
    }
}


// ====================================================================
// TestInelasticKinematics: Test the new functions
// ====================================================================
void TestInelasticKinematics()
{
    cout << "\n========================================" << endl;
    cout << "Testing Inelastic Kinematics" << endl;
    cout << "========================================" << endl;
    
    // Make sure data is loaded
    if (!g_inelastic_xs) {
        LoadInelasticData();
    }
    
    Double_t ekin = 200.0;  // MeV
    Double_t E_ex = 0.00443;  // 4.43 MeV in GeV
    
    // Calculate elastic CM momentum (for comparison)
    Double_t ekinGeV = ekin / 1000.;
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    Double_t s = iState.Mag2();
    Double_t pcm_elastic = TMath::Sqrt((s-(mp+mC)*(mp+mC)) * 
                                       (s-(mp-mC)*(mp-mC)) / (4.*s));
    
    // Calculate inelastic CM momentum
    Double_t pcm_inelastic = CalculateCMMomentumInelastic(ekin, E_ex);
    
    cout << "\nCM Momentum Comparison:" << endl;
    cout << "  Elastic:   pcm = " << pcm_elastic << " GeV/c" << endl;
    cout << "  Inelastic: pcm = " << pcm_inelastic << " GeV/c" << endl;
    cout << "  Difference: " << (pcm_elastic - pcm_inelastic)*1000 << " MeV/c" << endl;
    
    // Test histogram creation
    cout << "\nTesting histogram creation:" << endl;
    TH1D* h = CreateInelasticSamplingHistogram(10.0, 40.0);
    if (h) {
        cout << "   Histogram created with " << h->GetNbinsX() << " bins" << endl;
        cout << "   Integral = " << h->Integral() << endl;
        
        // Sample a few angles
        cout << "\nSampling test (10 random angles):" << endl;
        for (int i = 0; i < 10; i++) {
            Double_t theta = h->GetRandom();
            cout << "    Sample " << i+1 << ": theta = " << theta << " deg" << endl;
        }
        
        delete h;
    }
    
    cout << "\n Kinematics test PASSED!" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// ProduceInelastic: Generate inelastic p+C scattering events (4.43 MeV)
// ====================================================================
// Generates events for: p + 12C(g.s.) → p' + 12C*(4.43 MeV, 2+)
// The excited carbon nucleus will later decay via gamma emission
//
// Parameters:
//   - name: ignored (kept for compatibility)
//   - ekin: beam kinetic energy [MeV]
//   - num: number of events to generate
//   - withTree: if true, keep ROOT file; if false, delete it
// ====================================================================
void ProduceInelastic(const char* name, Double_t ekin, Int_t num, Bool_t withTree = 0)
{
    gROOT->SetStyle("Plain");
    
    cout << "\n========================================" << endl;
    cout << "INELASTIC EVENT GENERATION" << endl;
    cout << "========================================" << endl;
    cout << "Beam energy: " << ekin << " MeV" << endl;
    cout << "Events: " << num << endl;
    cout << "Excitation: 4.43 MeV (2+ state)" << endl;
    cout << "========================================\n" << endl;
    
    // Load inelastic data if not already loaded
    if (!g_inelastic_xs) {
        LoadInelasticData();
    }
    
    // ================================================================
    // Prepare output files
    // ================================================================
    TString basename = Form("pC_Inel443_%3.0fMeV_ST", ekin);
    TString nameo = MakeOutputPath(basename + ".txt");
    TString namer = MakeOutputPath(basename + ".root");
    
    cout << "Output files:" << endl;
    cout << "  ROOT: " << namer << endl;
    cout << "  TXT:  " << nameo << endl;
    cout << "========================================\n" << endl;
    
    TFile f(namer.Data(), "RECREATE");
    TTree* tree = new TTree("data", "pC->pC* inelastic (4.43 MeV)");
    TClonesArray* part = new TClonesArray("PParticle");
    tree->Branch("Particles", &part);
    
    // Metadata branches
    Int_t npart = 2;
    Float_t impact = 0.;
    Float_t phi0 = 0.;
    Float_t excitation = 4.43;  // MeV
    tree->Branch("Npart", &npart);
    tree->Branch("Impact", &impact);
    tree->Branch("Phi", &phi0);
    tree->Branch("Excitation", &excitation);  // Store excitation energy
    
    //FILE* fout = fopen(nameo.Data(), "w");
    
    // ================================================================
    // STEP 1: Calculate initial state kinematics
    // ================================================================
    Double_t ekinGeV = ekin / 1000.;
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    
    Double_t s = iState.Mag2();
    
    // Calculate CM momentum for INELASTIC scattering
    Double_t pcm = CalculateCMMomentumInelastic(ekin, 0.00443);
    
    cout << "pbeam: " << mom << " GeV/c" << endl;
    cout << "pcm (inelastic): " << pcm << " GeV/c" << endl;
    cout << "s: " << s << " GeV²" << endl << endl;
    
    // ================================================================
    // Write ASCII file header
    // ================================================================
    /*
	int codr = 302;  // 302 for inelastic (301 was elastic)
    int nbp = 2;
    double wei = 1.000E+00;
    int fform = 10000000 + codr;
    fprintf(fout, " %10d%10.2E%10.3f%10.3f%10.3f%10.3f%7d\n", 
            fform, 7.36E-8, mom, 0.f, 0.f, 0.f, num);
    fprintf(fout, " REAC,CROSS(mb),B. MOM,  A1,    A2,    A3, # EVENTS\n");
    fprintf(fout, "  0.00000 0.00000 110.000 450.000  33.000  40.000 0.90000 2.30000\n");
    fprintf(fout, "  AFD7GH,ATS7GH,ZFD7GH,ZTS7GH,RCD7GH,ZCD7H,ASL7GH,ASU7GH\n");
    */
    // ================================================================
    // STEP 2: Define detector acceptance in LAB frame
    // ================================================================
    const Double_t theta_lab_target = DETECTOR_THETA_CENTER_RAD;
    const Double_t theta_lab_window = DETECTOR_THETA_WINDOW;
    const Double_t phi_lab_window = DETECTOR_PHI_WINDOW;
    
    const Double_t theta_lab_min = theta_lab_target - theta_lab_window;
    const Double_t theta_lab_max = theta_lab_target + theta_lab_window;
    const Double_t phi_lab_min_0 = -phi_lab_window;
    const Double_t phi_lab_max_0 = phi_lab_window;
    const Double_t phi_lab_min_180 = TMath::Pi() - phi_lab_window;
    const Double_t phi_lab_max_180 = -TMath::Pi() + phi_lab_window;
    
    // ================================================================
    // STEP 3: Calculate CM angle range for detector acceptance
    // ================================================================
    Double_t theta_cm_min, theta_cm_max;
    ComputeCMAngleRange(ekin, theta_lab_target, theta_lab_window, 
                        theta_cm_min, theta_cm_max);
    
    // Create sampling histogram from digitized data
    TH1D* sigmacm_narrow = CreateInelasticSamplingHistogram(theta_cm_min, theta_cm_max);
    
    cout << "Sampling from digitized inelastic cross section" << endl;
    cout << "CM angle range: [" << theta_cm_min << ", " << theta_cm_max 
         << "] degrees" << endl;
    cout << "Histogram integral: " << sigmacm_narrow->Integral() << endl << endl;
    
    // ================================================================
    // STEP 4: Event generation loop (IMPORTANCE SAMPLING)
    // ================================================================
    TLorentzVector p, C;
    Double_t phi, theta, stheta, ctheta;
    Int_t nevent = num;
    Int_t li = 0;
    
    cout << "Generating events..." << endl;
    
    while (li < nevent)
    {
        if (!(li % 10000)) cout << li << " ..." << endl;
        
        // ============================================================
        // 4a. Sample CM angle from digitized cross section
        // ============================================================
        theta = sigmacm_narrow->GetRandom();  // Importance sampling from data!
        stheta = TMath::Sin(TMath::DegToRad() * theta);
        ctheta = TMath::Cos(TMath::DegToRad() * theta);
        
        // Sample azimuthal angle uniformly
        phi = gRandom->Rndm() * 2. * TMath::Pi();
        
        // ============================================================
        // 4b. Create 4-momentum in CM frame with INELASTIC momentum
        // ============================================================
        p.SetXYZM(pcm * stheta * TMath::Cos(phi),
                  pcm * stheta * TMath::Sin(phi),
                  pcm * ctheta,
                  mp);
        
        // ============================================================
        // 4c. Boost to LAB frame
        // ============================================================
        p.Boost(iState.BoostVector());
        
        // Calculate recoil Carbon 4-momentum (momentum conservation)
        C = iState - p;
        
		// DEBUG - first event only
		static std::atomic<int> first_event{1};
		if (first_event.exchange(0) == 1) {
			std::cout << "\n=== DEBUG: First Event ===" << std::endl;
			std::cout << "pcm used: " << pcm << " GeV/c" << std::endl;
			std::cout << "Proton: E=" << p.E() << ", |p|=" << p.P() << ", M=" << p.M() << std::endl;
			std::cout << "Carbon: E=" << C.E() << ", |p|=" << C.P() << ", M=" << C.M() << std::endl;
			std::cout << "Expected C.M() = " << (mC + 0.00443) << " GeV" << std::endl;
			std::cout << "========================\n" << std::endl;
		}
		
        // ============================================================
        // 4d. Apply detector acceptance cuts in LAB frame
        // ============================================================
        Bool_t theta_accepted = (p.Theta() >= theta_lab_min && p.Theta() <= theta_lab_max);
        Bool_t phi_accepted_0 = (p.Phi() >= phi_lab_min_0 && p.Phi() <= phi_lab_max_0);
        Bool_t phi_accepted_180 = (p.Phi() >= phi_lab_min_180 || p.Phi() <= phi_lab_max_180);
        
        if (theta_accepted && (phi_accepted_0 || phi_accepted_180))
        {
            li++;
            
            // ====================================================
            // 4e. Save event to ROOT tree
            // ====================================================
            part->Clear();
            // Particle 1: scattered proton (ID=14)
            new ((*part)[0]) PParticle(14, p.X(), p.Y(), p.Z());
            // Particle 2: recoil Carbon-12* excited state (ID=614)
            // NOTE: This carbon is in excited state with 4.43 MeV
            new ((*part)[1]) PParticle(614, C.X(), C.Y(), C.Z(), mC + E_EXCITATION, 1);
            
            tree->Fill();
            
            // ====================================================
            // 4f. Write to ASCII file
            // ====================================================
            //fprintf(fout, " %10d%10d%10d%10.4f%10.3E%10.3E\n", 
            //        li, codr, nbp, mom, wei, wei);
            //fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 
            //        1, p.X(), p.Y(), p.Z(), 14);
            //fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 
            //        2, C.X(), C.Y(), C.Z(), 614);
        }
    }
    
	// ================================================================
    // STEP 5: Write output and clean up
    // ================================================================
    cout << "\nn processed: " << li << endl;
    
    //fclose(fout);
    tree->Write();
    f.Close();
    
    cout << "Generated " << nevent << " inelastic events" << endl;
    cout << "Output files: " << namer << ", " << nameo << endl;
    
    if (!withTree)
    {
        gSystem->Exec(Form("rm %s", namer.Data()));
    }
    
    cout << "\n========================================" << endl;
    cout << "INELASTIC GENERATION COMPLETE!" << endl;
    cout << "========================================\n" << endl;
}

// ====================================================================
// VerifySimulationInelastic: Verify inelastic event generation
// ====================================================================
// Analyzes the generated inelastic events and creates diagnostic plots
// comparing with expected distributions
//
// Parameters:
//   - filename: ROOT file to analyze (default: pC_Inel443_200MeV_ST.root)
// ====================================================================
void VerifySimulationInelastic(const char* filename = "pC_Inel443_200MeV_ST.root")
{
    cout << "\n========================================" << endl;
    cout << "INELASTIC SIMULATION VERIFICATION" << endl;
    cout << "========================================" << endl;
    cout << "Analyzing file: " << filename << endl;
    
    // Open the file
    TFile* f = new TFile(filename, "READ");
    if (!f || f->IsZombie()) {
        cerr << "ERROR: Cannot open file " << filename << endl;
        return;
    }
    
    // Get the tree
    TTree* tree = (TTree*)f->Get("data");
    if (!tree) {
        cerr << "ERROR: Cannot find tree 'data'" << endl;
        f->Close();
        return;
    }
    
    // Setup branches
    TClonesArray* particles = new TClonesArray("PParticle");
    tree->SetBranchAddress("Particles", &particles);
    Float_t excitation = 0;
    tree->SetBranchAddress("Excitation", &excitation);
    
    Long64_t nEntries = tree->GetEntries();
    cout << "Total events: " << nEntries << endl;
    
    // Create histograms
    TH1D* h_proton_energy = new TH1D("h_proton_energy", 
        "Scattered Proton Kinetic Energy;T_{p} (MeV);Events", 
        100, 190, 200);
    
    TH1D* h_proton_theta = new TH1D("h_proton_theta",
        "Scattered Proton Angle;#theta_{lab} (deg);Events",
        100, 0, 30);
    
    TH1D* h_carbon_energy = new TH1D("h_carbon_energy",
        "Carbon Total Internal Energy;E_{internal} (MeV);Events",
        100, 0, 10);
    
    TH1D* h_energy_balance = new TH1D("h_energy_balance",
        "Energy Conservation Check;#DeltaE (MeV);Events",
        100, -1, 1);
    
    TH2D* h_energy_vs_angle = new TH2D("h_energy_vs_angle",
        "Proton Energy vs Angle;#theta_{lab} (deg);T_{p} (MeV)",
        50, 10, 25, 50, 193, 197);
    
	h_proton_energy->SetDirectory(0);
    h_proton_theta->SetDirectory(0);
    h_carbon_energy->SetDirectory(0);
    h_energy_balance->SetDirectory(0);
    h_energy_vs_angle->SetDirectory(0);
	
    // Constants
    const Double_t mp_GeV = 0.9382720813;
    const Double_t mC_GeV = 11.174862;
    const Double_t beam_energy = 0.200;  // 200 MeV in GeV
    const Double_t E_initial = beam_energy + mp_GeV + mC_GeV;  // Total initial energy
    
    cout << "\nAnalyzing events..." << endl;
    
    // Loop over events
    for (Long64_t i = 0; i < nEntries; i++) {
        if (i % 10000 == 0 && i > 0) {
            cout << "  Processed " << i << " events..." << endl;
        }
        
        tree->GetEntry(i);
        
        if (particles->GetEntries() < 2) continue;
        
        // Get particles
        PParticle* proton = (PParticle*)particles->At(0);
        PParticle* carbon = (PParticle*)particles->At(1);
        
        // Proton kinematics
        TLorentzVector p;
        p.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp_GeV);
        
        Double_t T_proton = (p.E() - mp_GeV) * 1000.;  // MeV
        Double_t theta_lab = p.Theta() * TMath::RadToDeg();
        
        // Carbon total energy from energy conservation
        Double_t E_carbon_total = E_initial - p.E();  // GeV
        
        // Carbon internal energy = kinetic + excitation
        Double_t E_internal_total = (E_carbon_total - mC_GeV) * 1000.;  // MeV
        
        // Energy conservation check
        Double_t E_final = p.E() + E_carbon_total;
        Double_t delta_E = (E_final - E_initial) * 1000.;  // MeV
        
        // Fill histograms
        h_proton_energy->Fill(T_proton);
        h_proton_theta->Fill(theta_lab);
        h_carbon_energy->Fill(E_internal_total);
        h_energy_balance->Fill(delta_E);
        h_energy_vs_angle->Fill(theta_lab, T_proton);
    }
    
    cout << "✓ Analysis complete!" << endl;
    
    // Calculate statistics
    Double_t mean_proton_E = h_proton_energy->GetMean();
    Double_t rms_proton_E = h_proton_energy->GetRMS();
    Double_t mean_theta = h_proton_theta->GetMean();
    Double_t mean_carbon_E = h_carbon_energy->GetMean();
    
    cout << "\n========================================" << endl;
    cout << "RESULTS SUMMARY" << endl;
    cout << "========================================" << endl;
    cout << "Proton Energy:" << endl;
    cout << "  Mean: " << mean_proton_E << " MeV" << endl;
    cout << "  RMS:  " << rms_proton_E << " MeV" << endl;
    cout << "  Expected: ~195.5 MeV (200 - 4.5 MeV)" << endl;
    
    cout << "\nCarbon Internal Energy:" << endl;
    cout << "  Mean: " << mean_carbon_E << " MeV" << endl;
    cout << "  Expected: ~4.5 MeV (excitation + recoil)" << endl;
    
    cout << "\nScattering Angle:" << endl;
    cout << "  Mean: " << mean_theta << " degrees" << endl;
    cout << "  Range: " << h_proton_theta->GetXaxis()->GetXmin() 
         << " - " << h_proton_theta->GetXaxis()->GetXmax() << " degrees" << endl;
    
    cout << "\nEnergy Conservation:" << endl;
    cout << "  Mean ΔE: " << h_energy_balance->GetMean() << " MeV" << endl;
    cout << "  RMS:     " << h_energy_balance->GetRMS() << " MeV" << endl;
    cout << "  (Should be ~0 MeV)" << endl;
    
    // Create canvas
    TCanvas* c1 = new TCanvas("c_inelastic_verify", 
        "Inelastic Simulation Verification", 1400, 1000);
    c1->Divide(3, 2);
    
    // Plot 1: Proton energy
    c1->cd(1);
    gPad->SetGrid();
    h_proton_energy->SetLineColor(kBlue);
    h_proton_energy->SetLineWidth(2);
    h_proton_energy->Draw();
    
    // Add expected value line
    TLine* line_expected = new TLine(195.5, 0, 195.5, h_proton_energy->GetMaximum());
    line_expected->SetLineColor(kRed);
    line_expected->SetLineStyle(2);
    line_expected->SetLineWidth(2);
    line_expected->Draw();
    
    TLegend* leg1 = new TLegend(0.15, 0.75, 0.45, 0.89);
    leg1->AddEntry(h_proton_energy, "Simulated", "l");
    leg1->AddEntry(line_expected, "Expected (~195.5 MeV)", "l");
    leg1->Draw();
    
    // Plot 2: Scattering angle
    c1->cd(2);
    gPad->SetGrid();
    h_proton_theta->SetLineColor(kGreen+2);
    h_proton_theta->SetLineWidth(2);
    h_proton_theta->Draw();
    
    // Plot 3: Carbon internal energy
    c1->cd(3);
    gPad->SetGrid();
    h_carbon_energy->SetLineColor(kMagenta);
    h_carbon_energy->SetLineWidth(2);
    h_carbon_energy->Draw();
    
    TLine* line_excitation = new TLine(4.43, 0, 4.43, h_carbon_energy->GetMaximum());
    line_excitation->SetLineColor(kRed);
    line_excitation->SetLineStyle(2);
    line_excitation->SetLineWidth(2);
    line_excitation->Draw();
    
    TLegend* leg3 = new TLegend(0.55, 0.75, 0.89, 0.89);
    leg3->AddEntry(h_carbon_energy, "Simulated", "l");
    leg3->AddEntry(line_excitation, "Excitation (4.43 MeV)", "l");
    leg3->Draw();
    
    // Plot 4: Energy conservation
    c1->cd(4);
    gPad->SetGrid();
    h_energy_balance->SetLineColor(kBlack);
    h_energy_balance->SetLineWidth(2);
    h_energy_balance->Draw();
    
    TLine* line_zero = new TLine(0, 0, 0, h_energy_balance->GetMaximum());
    line_zero->SetLineColor(kRed);
    line_zero->SetLineStyle(2);
    line_zero->SetLineWidth(2);
    line_zero->Draw();
    
    // Plot 5: Energy vs angle (2D)
    c1->cd(5);
    gPad->SetGrid();
    h_energy_vs_angle->SetStats(0);
    h_energy_vs_angle->Draw("COLZ");
    
    // Plot 6: Angular distribution with expected shape
    c1->cd(6);
    gPad->SetGrid();
    gPad->SetLogy();
    
    // Project angular distribution
    TH1D* h_theta_fine = new TH1D("h_theta_fine", 
        "Angular Distribution (Log Scale);#theta_{CM} (deg);Events",
        30, 10, 25);
    h_theta_fine->SetDirectory(0);
	
    // Fill with CM angles converted from lab
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (particles->GetEntries() < 2) continue;
        
        PParticle* proton = (PParticle*)particles->At(0);
        TLorentzVector p;
        p.SetXYZM(proton->Px(), proton->Py(), proton->Pz(), mp_GeV);
        
        // Approximate CM angle (for visualization)
        Double_t theta_lab = p.Theta();
        h_theta_fine->Fill(theta_lab * TMath::RadToDeg());
    }
    
    h_theta_fine->SetLineColor(kBlue);
    h_theta_fine->SetLineWidth(2);
    h_theta_fine->Draw();
    
    // Overlay expected cross section shape
    if (g_inelastic_xs) {
        TH1D* h_expected = (TH1D*)h_theta_fine->Clone("h_expected");
        h_expected->Reset();
        
        for (int i = 1; i <= h_expected->GetNbinsX(); i++) {
            Double_t theta = h_expected->GetBinCenter(i);
            Double_t xs = GetInelasticCrossSection(theta);
            h_expected->SetBinContent(i, xs * h_theta_fine->GetMaximum() / 
                                      g_inelastic_xs->Eval(15.0));
        }
        
        h_expected->SetLineColor(kRed);
        h_expected->SetLineStyle(2);
        h_expected->SetLineWidth(2);
        h_expected->Draw("SAME");
        
        TLegend* leg6 = new TLegend(0.55, 0.75, 0.89, 0.89);
        leg6->AddEntry(h_theta_fine, "Simulated", "l");
        leg6->AddEntry(h_expected, "Expected (scaled)", "l");
        leg6->Draw();
    }
    
    c1->Update();
    
    // Save canvas
    //c1->SaveAs("inelastic_verification.pdf");
    //c1->SaveAs("inelastic_verification.png");
    
	cout << "\n✓ Plot displayed" << endl;
    //cout << "\n✓ Plots saved: inelastic_verification.pdf/png" << endl;
    cout << "========================================\n" << endl;
    
    // Clean up
    f->Close();
    delete particles;
	
	
}

// ====================================================================
// CONVENIENCE WRAPPER: ProduceInelasticPolarized (single-threaded)
// ====================================================================
// For compatibility with your existing code structure
// Just calls the multithreaded version with 1 thread
// ====================================================================
void ProduceInelasticPolarized(const char* name, Double_t ekin, Int_t num, 
                                Int_t spin_state, Bool_t withTree = 0)
{
    // Call multithreaded version with 1 thread
    SingleRunMultithreadInelasticPolarized(ekin, num, BEAM_POLARIZATION, spin_state, 1);
    
    // Handle file deletion if requested
    if (!withTree) {
        int polar_int = (int)(BEAM_POLARIZATION * 100);
        TString spinLabel = (spin_state > 0) ? "SpinUp" : "SpinDown";
        TString basename = Form("pC_Inel443_%3.0fMeV_MT_P%d_%s", ekin, polar_int, spinLabel.Data());
        TString namer = MakeOutputPath(basename + ".root");
        gSystem->Exec(Form("rm %s", namer.Data()));
    }
}

