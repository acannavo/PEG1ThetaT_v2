// EventGenerator.C
// Implementation of multithreaded event generation

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
	
	TF1 *ffps_narrow = new TF1(Form("ffps_narrow_t%d", thread_id), xsFormula, 
	                            theta_cm_min*TMath::DegToRad(), 
	                            theta_cm_max*TMath::DegToRad());
	
	TH1D* sigmacm_narrow = new TH1D(Form("sigmacm_narrow_t%d", thread_id), 
	                                 Form("narrow range t%d", thread_id), 
	                                 10000, theta_cm_min, theta_cm_max);
	for (int ibin = 1; ibin <= 10000; ibin++) {
		Double_t theta_deg = theta_cm_min + (theta_cm_max - theta_cm_min) * (ibin-0.5) / 10000.;
		sigmacm_narrow->SetBinContent(ibin, ffps_narrow->Eval(theta_deg * TMath::DegToRad()));
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
// SingleRunMultithreadPolarized - Polarized parallel event generation
// ====================================================================
void SingleRunMultithreadPolarized(Double_t energy, Int_t number_total, 
                                    Double_t polarization, Int_t spin_state,
                                    Int_t num_threads)
{
	std::cout << "\n=== Multithreaded POLARIZED Event Generation ===" << std::endl;
	std::cout << "Energy: " << energy << " MeV" << std::endl;
	std::cout << "Total events: " << number_total << std::endl;
	std::cout << "Polarization: " << polarization * 100 << "%" << std::endl;
	std::cout << "Spin state: " << (spin_state > 0 ? "UP" : "DOWN") << std::endl;
	
	auto start = std::chrono::high_resolution_clock::now();
	
	if (num_threads == 0) {
		num_threads = std::thread::hardware_concurrency();
	}
	std::cout << "Number of threads: " << num_threads << std::endl;
	
	Int_t events_per_thread = number_total / num_threads;
	
	const Double_t theta_lab_target = DETECTOR_THETA_CENTER_RAD;
	const Double_t theta_lab_window = DETECTOR_THETA_WINDOW;
	const Double_t phi_lab_window   = DETECTOR_PHI_WINDOW;
	
	const char* xsFormula = GetXSFormula(energy);
	
	ROOT::TThreadExecutor executor(num_threads);
	
	auto thread_work = [energy, events_per_thread, theta_lab_target, 
	                     theta_lab_window, phi_lab_window, xsFormula,
	                     polarization, spin_state]
	                    (int tid) {
		return GenerateThreadEventsPolarized(tid, energy, events_per_thread, 
		                                      theta_lab_target, theta_lab_window, 
		                                      phi_lab_window, xsFormula,
		                                      polarization, spin_state);
	};
	
	std::cout << "Starting multithreaded generation..." << std::endl;
	auto results = executor.Map(thread_work, ROOT::TSeq<int>(0, num_threads));
	
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	
	std::cout << "\n=== Writing to ROOT file ===" << std::endl;
	
	
	// Prepare output files with standardized names
	int polar_int = (int)(polarization * 100);
	TString spinLabel = (spin_state > 0) ? "SpinUp" : "SpinDown";
	TString basename = Form("pC_Elas_%3.0fMeV_MT_P%d_%s_%dp%d", energy, polar_int, spinLabel.Data(), 
                        (int)DETECTOR_THETA_CENTER, (int)(DETECTOR_THETA_CENTER*10)%10);
	TString nameo = MakeOutputPath(basename + ".txt");
	TString namer = MakeOutputPath(basename + ".root");

	std::cout << "Output files:" << std::endl;
	std::cout << "  ROOT: " << namer << std::endl;
	std::cout << "  TXT:  " << nameo << std::endl;
	std::cout << "========================================\n" << std::endl;
		
	TFile f(namer.Data(), "RECREATE");
	TTree* tree = new TTree("data", "pC->pC cross section");
	TClonesArray* part = new TClonesArray("PParticle");
	tree->Branch("Particles", &part, 32000, 0);
	
	Int_t npart = 2;
	Float_t impact = 0.;
	Float_t phi0 = 0.;
	tree->Branch("Npart", &npart);
	tree->Branch("Impact", &impact);
	tree->Branch("Phi", &phi0);
	
	FILE* fout = fopen(nameo.Data(), "w");
	
	Double_t ekinGeV = energy/1000.;
	Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
	int codr = 301;
	int nbp = 2;
	double wei = 1.000E+00;
	int fform = 10000000 + codr;
	fprintf(fout, " %10d%10.2E%10.3f%10.3f%10.3f%10.3f%7d\n", fform, 7.36E-8, mom, 0.f, 0.f, 0.f, number_total);
	fprintf(fout, " REAC,CROSS(mb),B. MOM,  A1,    A2,    A3, # EVENTS\n");
	fprintf(fout, "  0.00000 0.00000 110.000 450.000  33.000  40.000 0.90000 2.30000\n");
	fprintf(fout, "  AFD7GH,ATS7GH,ZFD7GH,ZTS7GH,RCD7GH,ZCD7H,ASL7GH,ASU7GH\n");
	
	int global_event_id = 0;
	for (auto& result : results) {
		for (int i = 0; i < result.count; i++) {
			global_event_id++;
			
			part->Clear();
			new ((*part)[0]) PParticle(14, result.px[i], result.py[i], result.pz[i]);
			new ((*part)[1]) PParticle(614, result.cx[i], result.cy[i], result.cz[i], mC, 1);
			tree->Fill();
			
			fprintf(fout, " %10d%10d%10d%10.4f%10.3E%10.3E\n", global_event_id, codr, nbp, mom, wei, wei);
			fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 1, result.px[i], result.py[i], result.pz[i], 14);
			fprintf(fout, "%4d%10.4f%10.4f%10.4f%3d\n", 2, result.cx[i], result.cy[i], result.cz[i], 614);
		}
	}
	
	fclose(fout);
	tree->Write();
	f.Close();
	
	std::cout << "Total events written: " << global_event_id << std::endl;
	std::cout << "Output file: " << namer << std::endl;
	std::cout << "Time elapsed: " << elapsed.count() << " seconds" << std::endl;
	if (elapsed.count() > 0) {
		std::cout << "Events/second: " << global_event_id / (double)elapsed.count() << std::endl;
	}
	std::cout << "===================================\n" << endl;
}


    
    // cout << "  Spin-UP:   0° = " << up_0 << ",  180° = " << up_180 << ",  Ratio = " << (double)up_0/up_180 << endl;
    // cout << "  Spin-DOWN: 0° = " << down_0 << ",  180° = " << down_180 << ",  Ratio = " << (double)down_0/down_180 << endl;
    // cout << "\nExpected behavior:" << endl;
    // cout << "  Spin-UP should have MORE events at 0° (right)" << endl;
    // cout << "  Spin-DOWN should have MORE events at 180° (left)" << endl;
    // cout << "========================================\n" << endl;
// }





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
                                             Int_t num_threads)
{
    std::cout << "\n=== Multithreaded POLARIZED INELASTIC Event Generation ===" << std::endl;
    std::cout << "Energy: " << energy << " MeV" << std::endl;
    std::cout << "Total events: " << number_total << std::endl;
    std::cout << "Excitation: 4.43 MeV (2+ state)" << std::endl;
    std::cout << "Polarization: " << polarization * 100 << "%" << std::endl;
    std::cout << "Spin state: " << (spin_state > 0 ? "UP" : "DOWN") << std::endl;
    
    // Load inelastic data if not already loaded
    if (!g_inelastic_xs || !g_inelastic_AN) {
        LoadInelasticData();
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
    auto results = executor.Map(thread_work, ROOT::TSeq<int>(0, num_threads));
    
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    
    std::cout << "\n=== Writing to ROOT file ===" << std::endl;
    
    // Prepare output files with standardized names
    int polar_int = (int)(polarization * 100);
    TString spinLabel = (spin_state > 0) ? "SpinUp" : "SpinDown";
    TString basename = Form("pC_Inel443_%3.0fMeV_MT_P%d_%s_%dp%d", energy, polar_int, spinLabel.Data(), 
                        (int)DETECTOR_THETA_CENTER, (int)(DETECTOR_THETA_CENTER*10)%10);
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
