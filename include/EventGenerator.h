// EventGenerator.h
// Multithreaded event generation for elastic and inelastic scattering
//
// ====================================================================
// SAMPLING MODE DEFAULTS (see EventGenerator.C for full guide)
// ====================================================================
// DEFAULT for both elastic and inelastic: T_SAMPLING_MEYER
//
//   T_SAMPLING_MEYER  — valid 160–200 MeV, energy-interpolated XS + A_N
//   THETA_CM_SAMPLING — legacy; inelastic A_N ONLY valid at 200 MeV
//   T_SAMPLING        — elastic only, energy-independent reference
// ====================================================================

#ifndef EVENT_GENERATOR_H
#define EVENT_GENERATOR_H

#include "ThreadUtils.h"
#include "TString.h"
#include "TH1.h"       // TH1D — required for Meyer envelope pointer in declarations

// ====================================================================
// Sampling mode selector
// ====================================================================
enum SamplingMode {
    THETA_CM_SAMPLING = 0,   // Legacy: sample θ_CM from energy-specific XS
    T_SAMPLING = 1,          // Current: sample t (energy-independent splines)
    T_SAMPLING_MEYER = 2     // NEW: sample t (Meyer energy-dependent splines)
};

// ====================================================================
// GenerateThreadEventsElasticPolarized
// ====================================================================
// Worker function for elastic polarized event generation
// Each thread independently generates events using elastic cross sections
//
// Inputs:
//   - thread_id: unique thread identifier (for RNG seeding)
//   - ekin: beam kinetic energy [MeV]
//   - events_per_thread: number of events this thread should generate
//   - theta_lab_target: central detector angle [radians]
//   - theta_lab_window: detector angular acceptance [radians]
//   - phi_lab_window: detector azimuthal acceptance [radians]
//   - xsFormula: TF1 formula string for elastic cross section
//   - polarization: beam polarization (0 to 1)
//   - spin_state: +1 for spin-up, -1 for spin-down
//
// Returns:
//   - ThreadData structure with generated events
// ====================================================================
ThreadData GenerateThreadEventsElasticPolarized(
    int thread_id,
    Double_t ekin,
    Int_t events_per_thread,
    Double_t theta_lab_target,
    Double_t theta_lab_window,
    Double_t phi_lab_window,
    const char* xsFormula,
    Double_t polarization,
    Int_t spin_state
);

// ====================================================================
// GenerateThreadEventsInelasticPolarized
// ====================================================================
// Worker function for inelastic polarized event generation
// Each thread independently generates events using inelastic cross sections
//
// Inputs:
//   - thread_id: unique thread identifier (for RNG seeding)
//   - ekin: beam kinetic energy [MeV]
//   - events_per_thread: number of events this thread should generate
//   - theta_lab_target: central detector angle [radians]
//   - theta_lab_window: detector angular acceptance [radians]
//   - phi_lab_window: detector azimuthal acceptance [radians]
//   - polarization: beam polarization (0 to 1)
//   - spin_state: +1 for spin-up, -1 for spin-down
//
// Returns:
//   - ThreadData structure with generated events
//
// Note: Requires LoadInelasticData() to be called first!
// ====================================================================
ThreadData GenerateThreadEventsInelasticPolarized(
    int thread_id,
    Double_t ekin,
    Int_t events_per_thread,
    Double_t theta_lab_target,
    Double_t theta_lab_window,
    Double_t phi_lab_window,
    Double_t polarization,
    Int_t spin_state
);

// ====================================================================
// GenerateThreadEventsPolarized_TSampling
// ====================================================================
// Worker function for elastic polarized event generation using t-sampling
// 
// KEY DIFFERENCES from theta-based version:
//   - Samples from Mandelstam t instead of θ_CM
//   - Uses XStSpline(t) for cross section (energy-independent!)
//   - Uses APSpline(t) for analyzing power (energy-independent!)
//   - Converts t → θ_CM using actual beam energy
//
// Inputs:
//   - thread_id: unique thread identifier (for RNG seeding)
//   - ekin: beam kinetic energy [MeV]
//   - events_per_thread: number of events this thread should generate
//   - t_min, t_max: momentum transfer range [GeV²]
//   - theta_lab_target: central detector angle [radians]
//   - theta_lab_window: detector angular acceptance [radians]
//   - phi_lab_window: detector azimuthal acceptance [radians]
//   - polarization: beam polarization (0 to 1)
//   - spin_state: +1 for spin-up, -1 for spin-down
//
// Returns:
//   - ThreadData structure with generated events
// ====================================================================
ThreadData GenerateThreadEventsPolarized_TSampling(
    int thread_id,
    Double_t ekin,
    Int_t events_per_thread,
    Double_t t_min,
    Double_t t_max,
    Double_t theta_lab_target,
    Double_t theta_lab_window,
    Double_t phi_lab_window,
    Double_t polarization,
    Int_t spin_state
);

// ====================================================================
// GenerateThreadEventsPolarized_MeyerSampling
// ====================================================================
// Worker function for ELASTIC polarized event generation using the
// Meyer energy-dependent envelope sampling method.
//
// KEY DIFFERENCES from T_SAMPLING version:
//   - Receives a pre-built envelope TH1D* (MAX of XS_160/XS_200 for
//     all t) instead of building its own sampling histogram.
//   - Two-level accept/reject:
//       (1) Sample t from envelope using GetRandom()
//       (2) Accept with P = MeyerXS_Elastic(t, ekin) / envelope_bin
//   - Uses MeyerAP_Elastic(t, ekin) for energy-interpolated A_N.
//   - Fully thread-safe: uses only a thread-local TRandom3; the
//     envelope is treated as read-only shared data.
//
// Inputs:
//   - thread_id:          unique thread identifier (for RNG seeding)
//   - ekin:               beam kinetic energy [MeV]
//   - events_per_thread:  number of accepted events to produce
//   - envelope:           pre-built envelope histogram (READ-ONLY)
//                         caller owns it; do NOT delete inside worker
//   - theta_lab_target:   central detector polar angle [radians]
//   - theta_lab_window:   detector polar acceptance half-width [radians]
//   - phi_lab_window:     detector azimuthal acceptance half-width [radians]
//   - polarization:       beam polarization fraction (0 to 1)
//   - spin_state:         +1 for spin-up, -1 for spin-down
//
// Returns:
//   - ThreadData with accepted events
//
// Usage:
//   Build envelope once with CreateMeyerEnvelope_Elastic() before
//   launching threads; pass the same pointer to every thread;
//   delete envelope only after executor.Map() returns.
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
    Int_t spin_state
);

// ====================================================================
// GenerateThreadEventsInelastic_MeyerSampling
// ====================================================================
// Worker function for INELASTIC polarized event generation using the
// Meyer energy-dependent envelope sampling method.
//
// Mirrors GenerateThreadEventsPolarized_MeyerSampling but uses:
//   - MeyerXS_Inelastic(t, ekin)  for cross section sampling
//   - MeyerAP_Inelastic(t, ekin)  for analyzing power
//   - CalculateCMMomentumInelastic(ekin, 0.00443) for reduced pcm
//   - Carbon recoil stored with (mC + E_EXCITATION) mass, status=1
//
// Inputs:
//   - thread_id:          unique thread identifier (for RNG seeding)
//   - ekin:               beam kinetic energy [MeV]
//   - events_per_thread:  number of accepted events to produce
//   - envelope:           pre-built inelastic envelope (READ-ONLY)
//                         build with CreateMeyerEnvelope_Inelastic()
//   - theta_lab_target:   central detector polar angle [radians]
//   - theta_lab_window:   detector polar acceptance half-width [radians]
//   - phi_lab_window:     detector azimuthal acceptance half-width [radians]
//   - polarization:       beam polarization fraction (0 to 1)
//   - spin_state:         +1 for spin-up, -1 for spin-down
//
// Returns:
//   - ThreadData with accepted events
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
    Int_t spin_state
);

// ====================================================================
// SingleRunMultithreadPolarized
// ====================================================================
// High-level orchestrator for elastic polarized event generation.
// Dispatches to the correct worker based on SamplingMode:
//   THETA_CM_SAMPLING  → GenerateThreadEventsPolarized (legacy)
//   T_SAMPLING         → GenerateThreadEventsPolarized_TSampling
//   T_SAMPLING_MEYER   → GenerateThreadEventsPolarized_MeyerSampling
//
// For T_SAMPLING_MEYER, the envelope is built once here and shared
// read-only across all worker threads; it is deleted after Map().
//
// Inputs:
//   - energy:        beam kinetic energy [MeV]
//   - number_total:  total number of events to generate
//   - polarization:  beam polarization (0 to 1)
//   - spin_state:    +1 for spin-up, -1 for spin-down
//   - num_threads:   number of threads (0 = auto-detect)
//   - mode:          sampling mode selector (default T_SAMPLING)
//
// Output files:
//   - pC_Elas_<energy>MeV_MT_P<pol>_<Spin>.root
//   - pC_Elas_<energy>MeV_MT_P<pol>_<Spin>.txt
// ====================================================================
void SingleRunMultithreadPolarized(
    Double_t energy,
    Int_t number_total,
    Double_t polarization,
    Int_t spin_state,
    Int_t num_threads = 0,
    SamplingMode mode = T_SAMPLING_MEYER  // Default: Meyer energy-interpolated
);

// ====================================================================
// SingleRunMultithreadInelasticPolarized
// ====================================================================
// High-level orchestrator for inelastic polarized event generation.
// Accepts an optional SamplingMode parameter:
//   THETA_CM_SAMPLING  → GenerateThreadEventsInelasticPolarized (legacy)
//   T_SAMPLING_MEYER   → GenerateThreadEventsInelastic_MeyerSampling
//
// For T_SAMPLING_MEYER, the inelastic envelope is built once here
// and shared read-only across all worker threads.
//
// Inputs:
//   - energy:        beam kinetic energy [MeV]
//   - number_total:  total number of events to generate
//   - polarization:  beam polarization (0 to 1)
//   - spin_state:    +1 for spin-up, -1 for spin-down
//   - num_threads:   number of threads (0 = auto-detect)
//   - mode:          sampling mode selector (default T_SAMPLING_MEYER
//                    — the only correct choice away from 200 MeV)
//
// Output files:
//   - pC_Inel443_<energy>MeV_MT_P<pol>_<Spin>.root
//   - pC_Inel443_<energy>MeV_MT_P<pol>_<Spin>.txt
//
// Note: Legacy mode automatically calls LoadInelasticData() if needed.
// ====================================================================
void SingleRunMultithreadInelasticPolarized(
    Double_t energy,
    Int_t number_total,
    Double_t polarization,
    Int_t spin_state,
    Int_t num_threads = 0,
    SamplingMode mode = T_SAMPLING_MEYER   // Default: Meyer (THETA_CM only valid at 200 MeV)
);

#endif // EVENT_GENERATOR_H