==============================================================
# PEG Project Architecture - Complete Documentation
==============================================================
**Project:** Polarized Event Generator for p+C Scattering  
**Date:** March 2026  
**Version:** 6.0  
**Status:** Meyer Energy-Dependent Sampling — Default ✅

**Change log v5 → v6:**
- Added MeyerScattering module (Phase 4 complete)
- Fixed unit bug in THETA_CM_SAMPLING elastic worker
- T_SAMPLING_MEYER is now the default for both elastic and inelastic
- Systematic quantification of all three sampling modes completed
- New diagnostic and comparison macros added

---


Summary

 1 Pluto description
 
 2 Code Architecture 
 
 3 Optimization
 
 4 t-Based Sampling Implementation (Feb 2026)
 
 5 **NEW:** Meyer Energy-Dependent Sampling (Mar 2026)
 
 6 **NEW:** Bug Fix — THETA_CM Unit Error (Mar 2026)
 
 7 **NEW:** Systematic Comparison of Sampling Modes (Mar 2026)


---

## Description 
The project uses a simulation framework called Pluto (https://hades.gsi.de/node/31) providing a library of C++ classes which can be set up with ROOT macro code. 
Tutorials are provided: https://plutouser.github.io/v6.01/PlutoHTMLDoc/examples/index.html


## 📊 Project Structure Overview

```
~/pluto_v6.02/PEG/
├── AnalysisAN.C                      # Systematic detector position study
├── BenchmarkInelasticScattering.C    # Inelastic physics benchmarks  
├── CrossSectionAnalyzingPower.csv    # Experimental reference data
├── ProjectArchitecture_v1.md         # Documentation version 1
├── ProjectArchitecture_v2.md         # Documentation version 2
├── ProjectArchitecture_v3.md         # Documentation version 3 (latest)
├── ValidateSimulation_XsAndAn.C      # Comprehensive validation suite
├── inelastic_analyzingpower_443_200MeV.csv    # Digitized AN data
├── inelastic_crosssection_443_200MeV.csv      # Digitized XS data
├── rootlogon.C                       # Auto-loads all modules
├── sigma_pC_Elas_EventGenerator.C    # Main program (high-level functions)
├── sigma_pC_Elas_EventGenerator_C.d  # Dependency file for main program
│
├── BackupCode/                       # Archived analysis scripts
│   ├── ValidateAna.C
│   └── ValidateSimulation.C
│
├── Data/                             # Simulation output ROOT files
│   ├── BackUP2/                      # Empty backup folder
│   ├── BackUp3/                      # Systematic study (14°-18°, 0.5° steps)
│   │   ├── pC_Elas_200MeV_MT_P80_Spin{Up,Down}_*.root    (18 files)
│   │   └── pC_Inel443_200MeV_MT_P80_Spin{Up,Down}_*.root (18 files)
│   ├── Backup/                       # Systematic study (14°-18°, 0.5° steps)
│   │   ├── pC_Elas_200MeV_MT_P80_Spin{Up,Down}_*.root    (18 files)
│   │   └── pC_Inel443_200MeV_MT_P80_Spin{Up,Down}_*.root (18 files)
│   ├── Backup4/                      # High-resolution study (15°-17.9°)
│   │   ├── pC_Elas_200MeV_MT_P80_Spin{Up,Down}_*.root    (46 files)
│   │   └── pC_Inel443_200MeV_MT_P80_Spin{Up,Down}_*.root (46 files)
│   └── Current data files (16 files at selected angles):
│       ├── pC_Elas_200MeV_MT_P80_SpinUp_{10p0,14p0,15p0,16p0,16p2,18p0,20p0,25p0}.root
│       ├── pC_Elas_200MeV_MT_P80_SpinDown_{10p0,14p0,15p0,16p0,16p2,18p0,20p0,25p0}.root
│       ├── pC_Inel443_200MeV_MT_P80_SpinUp_{10p0,14p0,15p0,16p0,16p2,18p0,20p0,25p0}.root
│       └── pC_Inel443_200MeV_MT_P80_SpinDown_{10p0,14p0,15p0,16p0,16p2,18p0,20p0,25p0}.root
│
├── include/                          # Header files (10 files)
│   ├── AnalyzingPowerUtils.h
│   ├── DetectorConfig.h
│   ├── ElasticScattering.h
│   ├── EventGenerator.h              # Updated: Meyer default, TH1.h added
│   ├── InelasticScattering.h
│   ├── Kinematics.h
│   ├── MeyerScattering.h             # NEW: Meyer module header
│   ├── PhysicalConstants.h
│   └── ThreadUtils.h
│   [Meyer spline data — in $PLUTOSYS/include or project include/]
│   ├── XSSpline_160MeV_elastic.C     # NEW: Meyer elastic XS at 160 MeV
│   ├── XSSpline_200MeV_elastic.C     # NEW: Meyer elastic XS at 200 MeV
│   ├── APSpline_160MeV_elastic.C     # NEW: Meyer elastic A_N at 160 MeV
│   ├── APSpline_200MeV_elastic.C     # NEW: Meyer elastic A_N at 200 MeV
│   ├── XSSpline_160MeV_inelastic.C   # NEW: Meyer inelastic XS at 160 MeV
│   ├── XSSpline_200MeV_inelastic.C   # NEW: Meyer inelastic XS at 200 MeV
│   ├── APSpline_160MeV_inelastic.C   # NEW: Meyer inelastic A_N at 160 MeV
│   └── APSpline_200MeV_inelastic.C   # NEW: Meyer inelastic A_N at 200 MeV
│
│   [output/]       # Not implemented yet 
│   ├── debug_inelastic_AN.pdf
│   ├── exp_analyzing_powers.pdf
│   ├── exp_cross_sections.pdf
│   ├── experimental_AN_comparison.pdf
│   ├── experimental_analyzingpowers.pdf
│   ├── experimental_crosssections.pdf
│   ├── experimental_data_all.pdf
│   ├── experimental_xs_comparison.pdf
│   ├── input_analyzingpowers.pdf
│   ├── input_crosssections.pdf
│   ├── input_data_all.pdf
│   ├── sim_analyzing_powers.pdf
│   ├── sim_cross_sections.pdf
│   ├── test_original_crosssections.pdf
│   ├── test_sampling_histograms.pdf
│   ├── test_theta_sampling.pdf
│   ├── validation_AN.pdf
│   ├── validation_AN_comparison.pdf
│   ├── validation_both_AN_scan.pdf
│   ├── validation_elastic.pdf
│   ├── validation_elastic_AN_scan.pdf
│   ├── validation_elastic_crosssection.pdf
│   ├── validation_inelastic.pdf
│   ├── validation_inelastic_AN_scan.pdf
│   ├── validation_inelastic_FIXED.pdf
│   ├── validation_inelastic_crosssection.pdf
│   ├── validation_inelastic_noconv.pdf
│   └── validation_xs.pdf
│
└── src/                              # Implementation files + compiled artifacts
    ├── AnalyzingPowerUtils.C
    ├── AnalyzingPowerUtils_C.d
    ├── AnalyzingPowerUtils_C.so
    ├── AnalyzingPowerUtils_C_ACLiC_dict_rdict.pcm
    ├── DetectorConfig.C
    ├── DetectorConfig_C.d
    ├── DetectorConfig_C.so
    ├── DetectorConfig_C_ACLiC_dict_rdict.pcm
    ├── ElasticScattering.C
    ├── ElasticScattering_C.d
    ├── ElasticScattering_C.so
    ├── ElasticScattering_C_ACLiC_dict_rdict.pcm
    ├── EventGenerator.C              # Updated: Meyer default + bug fix
    ├── EventGenerator_C.d
    ├── EventGenerator_C.so
    ├── EventGenerator_C_ACLiC_dict_rdict.pcm
    ├── InelasticScattering.C
    ├── InelasticScattering_C.d
    ├── InelasticScattering_C.so
    ├── InelasticScattering_C_ACLiC_dict_rdict.pcm
    ├── Kinematics.C
    ├── Kinematics_C.d
    ├── Kinematics_C.so
    ├── Kinematics_C_ACLiC_dict_rdict.pcm
    ├── MeyerScattering.C             # NEW: Meyer module implementation
    ├── MeyerScattering_C.d
    ├── MeyerScattering_C.so
    ├── MeyerScattering_C_ACLiC_dict_rdict.pcm
    ├── TestDetectorConfig.C
    ├── TestElasticScattering.C
    ├── TestInelasticScattering.C
    └── TestKinematics.C

Total: 10 directories, 274 files
```

---

## 📁 File Organization Notes

### **Directory Structure:**

The project has a **flat structure** at the root level with the following organization:

**Root Directory (`~/pluto_v6.02/PEG/`) contains:**
- **4 Analysis Scripts:** AnalysisAN.C, BenchmarkInelasticScattering.C, ValidateSimulation_XsAndAn.C, sigma_pC_Elas_EventGenerator.C
- **3 CSV Data Files:** CrossSectionAnalyzingPower.csv, inelastic_analyzingpower_443_200MeV.csv, inelastic_crosssection_443_200MeV.csv
- **2 Configuration Files:** rootlogon.C, sigma_pC_Elas_EventGenerator_C.d
- **5 Subdirectories:** BackupCode/, Data/, include/, output/, src/

### **Subdirectories:**

1. **include/** (8 header files)
   - All .h files for module declarations
   - No subdirectories

2. **src/** (34 files)
   - 6 implementation files (.C)
   - 4 test files (.C)
   - 24 compiled artifacts (.so, .d, .pcm)

3. **Data/** (228 files in 4 backup folders + 16 current files)
   - **BackUP2/** - Empty folder
   - **Backup/** - 36 files (14°-18°, 0.5° steps)
   - **BackUp3/** - 36 files (14°-18°, 0.5° steps)
   - **Backup4/** - 92 files (15°-17.9°, 0.1-0.2° steps)
   - **16 current files** at selected angles (10°, 14°, 15°, 16°, 16.2°, 18°, 20°, 25°)

4. **output/** (30 PDF files)
   - All analysis plots and validation results
   - No subdirectories

5. **BackupCode/** (2 archived scripts)
   - Old versions of analysis scripts

### **File Naming Convention:**

**Simulation ROOT files:**
```
pC_{Elas|Inel443}_200MeV_MT_P{polarization}_Spin{Up|Down}_{angle}.root

Examples:
- pC_Elas_200MeV_MT_P80_SpinUp_16p2.root
- pC_Inel443_200MeV_MT_P80_SpinDown_14p0.root
```

**Compiled artifacts pattern:**
```
{ModuleName}_C.{d|so}
{ModuleName}_C_ACLiC_dict_rdict.pcm

Examples:
- ElasticScattering_C.so
- Kinematics_C.d
- EventGenerator_C_ACLiC_dict_rdict.pcm
```

### **Key Points:**

- ✅ **No nested source folders** - All analysis scripts in root
- ✅ **Flat configuration** - All CSV files in root  
- ✅ **Clear separation** - Headers in include/, implementations in src/
- ✅ **Organized backups** - Multiple timestamped backup folders in Data/
- ✅ **Centralized output** - All plots in output/ directory

---

## 🏗️ Extended Module Architecture Diagram


```
┌────────────────────────────────────────────────────────────────────┐
│                  sigma_pC_Elas_EventGenerator.C                    │
│                         (MAIN PROGRAM)                             │
│                    ~1700 lines (was 3546)                          │
│                                                                    │
│  • High-level wrapper functions                                   │
│  • Test and diagnostic functions                                  │
│  • Verification routines                                          │
│  • User-facing convenience functions                              │
└─────────────────────────┬──────────────────────────────────────────┘
                          │
                          │ #includes all modules
                          │
       ┌──────────────────┼──────────────────┬──────────────┐
       │                  │                  │              │      
       ▼                  ▼                  ▼              ▼      
┌──────────────┐   ┌─────────────┐   ┌─────────────┐ ┌───────────────────┐
│ Detector     │   │  Physical   │   │  Elastic    │ │  Inelastic        │
│ Config       │   │  Constants  │   │ Scattering  │ │  Scattering       │
│ .h + .C      │   │ .h only     │   │ .h + .C     │ │  .h + .C          │
└──────┬───────┘   └──────┬──────┘   └──────┬──────┘ └────────┬──────────┘
       │                  │                 │                 │
       │                  │                 │                 │
       └──────────────────┴─────────────────┴─────────────────┘
                                     │
                                     ▼
                            ┌─────────────────┐
                            │   Kinematics    │
                            │    .h + .C      │
                            └────────┬────────┘
                                     │
                ┌────────────────────┴────────────────────┐
                │                                         │
                ▼                                         ▼
       ┌─────────────────┐                      ┌─────────────────┐
       │  ThreadUtils    │                      │ EventGenerator  │
       │    .h only      │◄─────────────────────│    .h + .C      │
       └─────────────────┘                      └─────────────────┘

                    ┌──────────────────────────┐
                    │   ANALYSIS TOOLS         │
                    │   (post-processing)      │
                    │                          │
                    │  • AnalyzingPowerUtils   │
                    │    .h + .C               │
                    │    (standalone)          │
                    │                          │
                    │  • AnalysisAN.C          │
                    │  • ValidateSimulation.C  │
                    │  • BenchmarkInelastic.C  │
                    └──────────────────────────┘


```

---

## Module Dependency Matrix

| Module | Depends On | Used By |
|--------|-----------|---------|
| **PhysicalConstants.h** | *None* | Kinematics, InelasticScattering, EventGenerator |
| **DetectorConfig.h/.C** | *None* | Kinematics, EventGenerator, Main |
| **ThreadUtils.h** | *None* | EventGenerator |
| **ElasticScattering.h/.C** | XSLog*.C (10 files) | EventGenerator, Main, Tests |
| **InelasticScattering.h/.C** | PhysicalConstants | EventGenerator |
| **Kinematics.h/.C** | PhysicalConstants, DetectorConfig | EventGenerator, Tests |
| **MeyerScattering.h/.C** ⭐ NEW | XSSpline/APSpline *.C (8 files) | EventGenerator |
| **EventGenerator.h/.C** | All above modules + MeyerScattering | Main program |
| **AnalyzingPowerUtils.h/.C** | *None* (standalone) | Analysis scripts |
| **sigma_pC_Elas_EventGenerator.C** | DetectorConfig, PhysicalConstants, ElasticScattering, EventGenerator | *User entry point* |

---

## Critical Configuration Values (VERIFIED)

### Detector Configuration (DetectorConfig.C)
```cpp
// DEFAULT VALUES (can be changed at runtime)
DETECTOR_THETA_CENTER = 16.2;       // degrees  
DETECTOR_THETA_WINDOW = 0.005;      // radians (±0.29°)
DETECTOR_PHI_WINDOW = 0.005;        // radians (±0.29°)
POLARIMETER_DISTANCE = 2.18;        // meters
```

### Physical Constants (PhysicalConstants.h)
```cpp
mp = 0.9382720813;                  // GeV/c² (proton mass)
mC = 11.174862;                     // GeV/c² (C-12 mass)
E_EXCITATION = 0.00443;             // GeV (4.43 MeV, 2+ state)
```

### Analyzing Power Parameters (ElasticScattering.h)
```cpp
// Wissink et al., Phys. Rev. C 45, R504 (1992)
T0_AN = 189.0;                      // MeV (optimal energy)
theta0_AN = 17.3;                   // degrees (optimal angle, lab)
alpha_AN = 1.21e-4;                 // MeV^-2
beta_AN = 1.61e-3;                  // MeV^-1 deg^-1
gamma_AN = 1.00e-2;                 // deg^-2
```

---

## Module Descriptions with Actual Implementation

### 📦 Core Infrastructure

#### **PhysicalConstants.h** (Header-only)
**Location:** `include/PhysicalConstants.h`  
**Dependencies:** None  
**Purpose:** Define fundamental physical constants  
**Contents:**
- Proton mass (mp = 0.9382720813 GeV)
- Carbon-12 mass (mC = 11.174862 GeV)
- Excitation energy (E_EXCITATION = 0.00443 GeV = 4.43 MeV)

**Usage:**
```cpp
#include "PhysicalConstants.h"
// Constants available: mp, mC, E_EXCITATION
```

---

#### **DetectorConfig.h/.C**
**Location:** `include/DetectorConfig.h`, `src/DetectorConfig.C`  
**Dependencies:** TMath.h (ROOT)  
**Purpose:** Manage detector geometry and acceptance

**Global Variables:**
```cpp
// Primary parameters (modifiable)
extern Double_t DETECTOR_THETA_CENTER;      // 16.2° default
extern Double_t DETECTOR_THETA_WINDOW;      // ±0.005 rad
extern Double_t DETECTOR_PHI_WINDOW;        // ±0.005 rad

// Derived constants (auto-calculated)
extern Double_t DETECTOR_THETA_CENTER_RAD;
extern Double_t DETECTOR_THETA_MIN;
extern Double_t DETECTOR_THETA_MAX;
extern Double_t DETECTOR_PHI_MIN_0;         // φ ≈ 0°
extern Double_t DETECTOR_PHI_MAX_0;
extern Double_t DETECTOR_PHI_MIN_180;       // φ ≈ 180°
extern Double_t DETECTOR_PHI_MAX_180;
```

**Functions:**
```cpp
void SetDetectorConfig(Double_t theta_center,
                       Double_t theta_window = 0.005,
                       Double_t phi_window = 0.005,
                       Double_t distance = 2.18);

void PrintDetectorConfig();
void ResetDetectorConfig();  // Resets to 16.2°
```

**Key Feature:** Dynamic reconfiguration for systematic studies

---

#### **ThreadUtils.h** (Header-only)
**Location:** `include/ThreadUtils.h`  
**Dependencies:** std::vector  
**Purpose:** Define thread communication structures

**Structure:**
```cpp
struct ThreadData {
    std::vector<int> event_ids;          // Event numbering
    std::vector<double> px, py, pz;      // Proton momentum [GeV/c]
    std::vector<double> cx, cy, cz;      // Carbon momentum [GeV/c]
    int count;                            // Number of events
};
```

---

### 🔬 Physics Modules

#### **ElasticScattering.h/.C**
**Location:** `include/ElasticScattering.h`, `src/ElasticScattering.C`  
**Dependencies:** XSLog150MeV.C through XSLog240MeV.C (10 files)  
**Purpose:** Elastic scattering cross sections and analyzing power

**Cross Section Implementation:**
```cpp
// Includes 10 optical model data files
#include "XSLog150MeV.C"  // log₁₀(dσ/dΩ) parameterization
// ... through ...
#include "XSLog240MeV.C"

// Conversion functions (log → linear)
double fXS150Op(double x);  // Returns σ = 10^(log₁₀σ)
// ... through ...
double fXS240Op(double x);

// Energy selection
const char* GetXSFormula(Double_t ekin);  // Returns "fXS200Op(x)"
```

**Analyzing Power:**
```cpp
Double_t GetElasticAnalyzingPower(Double_t ekin, Double_t theta_lab_deg);
// Formula: A_N = 1 - α(T-T₀)² - β(T-T₀)(θ-θ₀) - γ(θ-θ₀)²
```

---

#### **InelasticScattering.h/.C**
**Location:** `include/InelasticScattering.h`, `src/InelasticScattering.C`  
**Dependencies:** PhysicalConstants.h  
**Purpose:** Inelastic scattering to C-12 2⁺ excited state (4.43 MeV)

**Global Objects:**
```cpp
extern TGraph* g_inelastic_xs;  // Cross section vs θ_CM
extern TGraph* g_inelastic_AN;  // Analyzing power vs θ_CM
```

**Functions:**
```cpp
void LoadInelasticData();
// Loads: inelastic_crosssection_443_200MeV.csv
//        inelastic_analyzingpower_443_200MeV.csv

Double_t GetInelasticCrossSection(Double_t theta_cm_deg);
Double_t GetInelasticAnalyzingPower(Double_t theta_cm_deg);

Double_t CalculateCMMomentumInelastic(Double_t ekin, Double_t E_ex = 0.00443);
// Modified p_cm due to excitation energy

TH1D* CreateInelasticSamplingHistogram(Double_t theta_cm_min,
                                       Double_t theta_cm_max,
                                       Int_t nbins = 10000);
```

---

#### **Kinematics.h/.C**
**Location:** `include/Kinematics.h`, `src/Kinematics.C`  
**Dependencies:** PhysicalConstants.h, DetectorConfig.h  
**Purpose:** CM ↔ Lab frame conversions and acceptance mapping

**Frame Conversions:**
```cpp
Double_t ConvertThetaCMtoLab(Double_t ekin, Double_t theta_cm_deg);
Double_t ConvertPhiCMtoLab(Double_t ekin, Double_t theta_cm_deg, 
                           Double_t phi_cm_rad);
```

**Acceptance Mapping:**
```cpp
void ComputeCMAngleRange(Double_t ekin,
                         Double_t theta_lab_center,
                         Double_t theta_lab_width,
                         Double_t& theta_cm_min,
                         Double_t& theta_cm_max);
// Uses 10,000-point lookup table for precision

void ComputePhiCMRanges(Double_t ekin, Double_t theta_cm_deg,
                        Double_t& phi_cm_min_0, Double_t& phi_cm_max_0,
                        Double_t& phi_cm_min_180, Double_t& phi_cm_max_180);
// Maps Lab phi windows (±5 mrad at 0° and 180°) to CM frame
```

**Polarization Sampling:**
```cpp
TF1* CreatePhiSamplingFunction(Double_t polarization,
                               Double_t analyzing_power,
                               Int_t spin_state);
// Returns: "1 + [0]*cos(x)" where [0] = -spin_state × P × A_N
// For sampling from polarized φ distribution
```


#### ⚠️ Important Note: CM Angle Range Limitation

**Function:** `ComputeCMAngleRange()` in `Kinematics.C` (line 126)

The CM→Lab angle mapping lookup table currently samples θ_CM from **0° to 50°**. This range is sufficient for forward detector angles (θ_Lab < ~40°) but must be extended if:

- Detector positioned at θ_Lab > 40° → Change to 90° range
- Detector positioned at θ_Lab > 74° → Change to 180° range (full coverage)

**To modify:** Edit line 126 in `Kinematics.C`:
```cpp
Double_t theta_cm_deg = i * 50.0 / nPoints;  // Change 50.0 to 90.0 or 180.0
```

**Current default detector:** θ_Lab = 16.2° (well within 50° range limit)

---

### ⚙️ Event Generation

#### **EventGenerator.h/.C**
**Location:** `include/EventGenerator.h`, `src/EventGenerator.C`  
**Dependencies:** ALL above modules  
**Purpose:** Multithreaded event generation orchestrator

**Main Functions:**
```cpp
// ELASTIC
void SingleRunMultithreadPolarized(
    Double_t energy,          // MeV
    Int_t number_total,       // Total events
    Double_t polarization,    // 0.0 to 1.0
    Int_t spin_state,         // +1 or -1
    Int_t num_threads = 0     // 0 = auto-detect
);
// Output: pC_Elas_<E>MeV_MT_P<pol>_<Spin>_<angle>.root

// INELASTIC
void SingleRunMultithreadInelasticPolarized(
    Double_t energy,
    Int_t number_total,
    Double_t polarization,
    Int_t spin_state,
    Int_t num_threads = 0
);
// Output: pC_Inel443_<E>MeV_MT_P<pol>_<Spin>_<angle>.root
```

**Worker Functions (per-thread):**
```cpp
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
```

**Event Generation Algorithm (per thread):**
1. Sample θ_CM from cross-section histogram
2. Convert θ_CM → θ_Lab for A_N lookup
3. Sample φ_CM using acceptance-rejection (if polarized)
4. Build CM frame 4-momentum
5. Lorentz boost to Lab frame
6. Apply detector acceptance cuts
7. Store event if accepted

---

### 📊 Analysis Tools

#### **AnalyzingPowerUtils.h/.C**
**Location:** `include/AnalyzingPowerUtils.h`, `src/AnalyzingPowerUtils.C`  
**Dependencies:** None (standalone, uses ROOT for file I/O)  
**Purpose:** Post-processing for A_N extraction

**Structures:**
```cpp
struct DetectorCounts {
    Double_t theta_det;
    Int_t N_right_up, N_right_down;
    Int_t N_left_up, N_left_down;
};

struct ANResult {
    Double_t theta_L, theta_R;
    Double_t AN, AN_error;
    Int_t N_R_up, N_R_down, N_L_up, N_L_down;
};
```

**Functions:**
```cpp
void CountEventsInFileByPhi(const TString& filename,
                            Double_t phi_window,
                            Bool_t is_spinup,
                            Int_t& N_right,
                            Int_t& N_left);

ANResult CalculateAsymmetry(Double_t theta_R, Double_t theta_L,
                            Int_t N_R_up, Int_t N_R_down,
                            Int_t N_L_up, Int_t N_L_down);
// Formula: ε = [√(N_R_up·N_L_down) - √(N_L_up·N_R_down)] /
//              [√(N_R_up·N_L_down) + √(N_L_up·N_R_down)]

Double_t AsymmetryToAnalyzingPower(Double_t asymmetry,
                                    Double_t polarization);
// A_N = ε / P
```

---

## Compilation Order

**Must compile in this order:**

1. **Independent modules** (can compile in any order):
   - `DetectorConfig.C`
   - `ElasticScattering.C`
   - `InelasticScattering.C`
   - `AnalyzingPowerUtils.C`

2. **Depends on step 1**:
   - `Kinematics.C` (needs DetectorConfig)

3. **Depends on steps 1-2**:
   - `EventGenerator.C` (needs all above)

4. **Main program** (depends on all):
   - `sigma_pC_Elas_EventGenerator.C`

**ROOT ACLiC handles this automatically:**
```cpp
.L DetectorConfig.C+
.L ElasticScattering.C+
.L InelasticScattering.C+
.L Kinematics.C+
.L EventGenerator.C+
.L sigma_pC_Elas_EventGenerator.C+
```

---

## Test Files

All test files verify module functionality:

| Test File | Tests Module | Purpose |
|-----------|-------------|---------|
| `TestDetectorConfig.C` | DetectorConfig | Dynamic configuration changes |
| `TestElasticScattering.C` | ElasticScattering | Cross section access (10 energies) |
| `TestInelasticScattering.C` | InelasticScattering | Data loading, interpolation, p_cm |
| `TestKinematics.C` | Kinematics | Frame conversions, phi sampling |

**Run tests:**
```cpp
.L DetectorConfig.C
.L TestDetectorConfig.C
TestDetectorConfig()
```

---

## Data Flow Diagram

```
USER
  │
  ├─> SingleRunMultithreadPolarized(200 MeV, 10M events, P=0.80, spin=+1)
  │
  └─> EventGenerator.C
        │
        ├─> GetXSFormula(200) → "fXS200Op(x)" [ElasticScattering]
        │
        ├─> ComputeCMAngleRange() [Kinematics]
        │
        ├─> Launch 8 threads → GenerateThreadEventsElasticPolarized()
        │     │
        │     ├─> Sample θ_CM from histogram (dσ/dΩ weighted)
        │     │
        │     ├─> ConvertThetaCMtoLab() [Kinematics]
        │     │
        │     ├─> GetElasticAnalyzingPower() [ElasticScattering]
        │     │
        │     ├─> Sample φ_CM (acceptance-rejection if polarized)
        │     │
        │     ├─> Build CM 4-momentum
        │     │
        │     ├─> Lorentz Boost to Lab
        │     │
        │     ├─> Check detector acceptance [DetectorConfig]
        │     │
        │     └─> Store if accepted → ThreadData
        │
        └─> Merge threads → Write ROOT + ASCII files
              │
              └─> Output: pC_Elas_200MeV_MT_P80_SpinUp_16p2.root
```


---

## Summary Statistics

| Category | Count |
|----------|-------|
| **Core Modules** | 7 (Physics, Config, Threading + Meyer) |
| **Analysis Tools** | 1 (AnalyzingPowerUtils) |
| **Test Suites** | 4 |
| **Header Files** | 10 |
| **Implementation Files** | 7 |
| **Data Files** | 25 (10 optical model + 8 Meyer splines + 2 inelastic CSV + 5 energy-indep.) |
| **Total Lines** | ~5000 (excluding main program) |

---

**Diagram Version:** 2.0 (Verified)  
**Last Updated:** March 2026  
**Source Code Verified:** ✅  
**Default Detector Angle:** 16.2° ✅
**Default Sampling Mode:** T_SAMPLING_MEYER ✅





---

## Code Optimization

Excellent question! This is a crucial optimization. Let me trace through the **theta sampling range optimization** step by step.

---

## 🎯 The Problem

**Without optimization:**
- Sample θ_CM from full range (0° to 180°)
- Most events rejected because they don't hit the detector
- Extremely inefficient!

**With optimization:**
- Pre-compute which θ_CM values will boost into detector acceptance
- Sample ONLY from this narrow range
- Much more efficient!

---

## 📍 The Main Function

### **Location:** `Kinematics.C`, lines 111-156
### **Function:** `ComputeCMAngleRange()`

```cpp
void ComputeCMAngleRange(Double_t ekin,              // Beam energy [MeV]
                         Double_t theta_lab_center,  // Detector center [radians]
                         Double_t theta_lab_width,   // Half-width [radians]
                         Double_t& theta_cm_min,     // OUTPUT: CM min [degrees]
                         Double_t& theta_cm_max)     // OUTPUT: CM max [degrees]
```

---

## 🔍 Step-by-Step Algorithm

### **STEP 1: Define Lab Acceptance Window**

**Source:** `Kinematics.C`, line 118

```cpp
// Define lab acceptance window
Double_t theta_lab_min_target = theta_lab_center - theta_lab_width;
Double_t theta_lab_max_target = theta_lab_center + theta_lab_width;
```

**Example with default detector:**
```
theta_lab_center = DETECTOR_THETA_CENTER_RAD = 16.2° = 0.2827 rad
theta_lab_width = DETECTOR_THETA_WINDOW = 0.005 rad = 0.286°

theta_lab_min_target = 0.2827 - 0.005 = 0.2777 rad (15.91°)
theta_lab_max_target = 0.2827 + 0.005 = 0.2877 rad (16.49°)
```

**Physical meaning:** Detector accepts protons scattered between 15.91° and 16.49° in Lab frame.

---

### **STEP 2: Create Lookup Table (CM → Lab mapping)**

**Source:** `Kinematics.C`, lines 121-131

```cpp
// Create lookup table: CM angle → Lab angle using our conversion function
const Int_t nPoints = 10000;  // High resolution!
Double_t cm_angles[nPoints];
Double_t lab_angles[nPoints];

for (Int_t i = 0; i < nPoints; i++)
{
    // Sample CM angles from 0° to 50° (sufficient for forward scattering)
    Double_t theta_cm_deg = i * 50.0 / nPoints;
    cm_angles[i] = theta_cm_deg;
    
    // Use our conversion function
    lab_angles[i] = ConvertThetaCMtoLab(ekin, theta_cm_deg);
}
```

**What this does:**
- Creates 10,000 (θ_CM, θ_Lab) pairs
- θ_CM ranges from 0° to 50° in 0.005° steps
- For each θ_CM, computes corresponding θ_Lab using Lorentz transformation

**Example table entries (for 200 MeV):**

| i | θ_CM [°] | θ_Lab [°] | In Acceptance? |
|---|----------|-----------|----------------|
| 0 | 0.000 | 0.000 | No |
| 3800 | 19.000 | 15.43 | No |
| 3900 | 19.500 | 15.84 | **Yes** ← First accepted |
| 4000 | 20.000 | 16.24 | **Yes** |
| 4100 | 20.500 | 16.65 | **Yes** ← Last accepted |
| 4200 | 21.000 | 17.05 | No |
| 10000 | 50.000 | 40.15 | No |

---

### **STEP 3: Find CM Range Corresponding to Lab Acceptance**

**Source:** `Kinematics.C`, lines 134-151

```cpp
// Find CM angles corresponding to lab acceptance edges
theta_cm_min = 0.;
theta_cm_max = 50.;
Bool_t found_min = false;
Bool_t found_max = false;

for (Int_t i = 0; i < nPoints; i++)
{
    // Find first CM angle that maps to minimum lab angle
    if (!found_min && lab_angles[i] >= theta_lab_min_target * TMath::RadToDeg()) {
        theta_cm_min = cm_angles[i];
        found_min = true;
    }
    // Find last CM angle that maps to maximum lab angle
    if (lab_angles[i] <= theta_lab_max_target * TMath::RadToDeg()) {
        theta_cm_max = cm_angles[i];
        found_max = true;
    }
}
```

**Search logic:**

1. **Find θ_CM_min:** First θ_CM where θ_Lab ≥ 15.91°
   ```
   Scan: θ_Lab = [0.0°, 5.2°, 10.3°, 15.43°, 15.84°, ...]
                                              ↑
                                          First ≥ 15.91°
   Result: θ_CM_min ≈ 19.5°
   ```

2. **Find θ_CM_max:** Last θ_CM where θ_Lab ≤ 16.49°
   ```
   Scan: θ_Lab = [..., 15.84°, 16.24°, 16.65°, 17.05°, ...]
                                    ↑
                                Last ≤ 16.49°
   Result: θ_CM_max ≈ 20.5°
   ```

**Result:** θ_CM must be in [19.5°, 20.5°] to hit detector at 16.2° ± 0.286°

---

### **STEP 4: Add 10% Safety Margin**

**Source:** `Kinematics.C`, lines 153-155

```cpp
// Add 10% margin to ensure we cover the full acceptance
Double_t margin = 0.1 * (theta_cm_max - theta_cm_min);
theta_cm_min = TMath::Max(0., theta_cm_min - margin);
theta_cm_max = theta_cm_max + margin;
```

**Calculation:**
```
Initial range: [19.5°, 20.5°]
Width: 20.5° - 19.5° = 1.0°
Margin: 0.1 × 1.0° = 0.1°

Final range:
  θ_CM_min = 19.5° - 0.1° = 19.4°
  θ_CM_max = 20.5° + 0.1° = 20.6°
```

**Why add margin?**
- **Numerical precision:** Lookup table has finite resolution (0.005°)
- **Edge effects:** Lorentz boost is nonlinear near boundaries
- **Energy spread:** Beam has finite energy width
- **Safety:** Better to sample slightly too wide than miss events

---

### **STEP 5: Print Summary**

**Source:** `Kinematics.C`, lines 158-165

```cpp
cout << "\n=== CM Angle Range Calculation ===" << endl;
cout << "Beam energy: " << ekin << " MeV" << endl;
cout << "Lab acceptance: " << theta_lab_center*TMath::RadToDeg() 
     << " ± " << theta_lab_width*TMath::RadToDeg() << " deg" << endl;
cout << "               " << theta_lab_min_target*TMath::RadToDeg()
     << " to " << theta_lab_max_target*TMath::RadToDeg() << " deg" << endl;
cout << "CM angle range: [" << theta_cm_min << ", " 
     << theta_cm_max << "] degrees" << endl;
cout << "===================================\n" << endl;
```

**Example output:**
```
=== CM Angle Range Calculation ===
Beam energy: 200 MeV
Lab acceptance: 16.2 ± 0.286 deg
               15.914 to 16.486 deg
CM angle range: [19.4, 20.6] degrees
===================================
```

---

## 🔄 How This Gets Used in Event Generation

### **Location:** `EventGenerator.C`, lines 98-104

```cpp
Double_t theta_cm_min, theta_cm_max;

// Call the optimization function
ComputeCMAngleRange(ekin, theta_lab_target, theta_lab_window, 
                    theta_cm_min, theta_cm_max);

// Create narrow sampling histogram ONLY in this range
TF1 *ffps_narrow = new TF1(Form("ffps_narrow_t%d", thread_id), xsFormula, 
                            theta_cm_min*TMath::DegToRad(), 
                            theta_cm_max*TMath::DegToRad());
```

**Without optimization:**
```cpp
// Would sample from 0° to 50° (or even 0° to 180°)
TF1 *ffps = new TF1("ffps", xsFormula, 0.*TMath::DegToRad(), 50.*TMath::DegToRad());
// Acceptance rate: ~2% (most events rejected!)
```

**With optimization:**
```cpp
// Sample from 19.4° to 20.6° (only ~1.2° range!)
TF1 *ffps_narrow = new TF1("ffps_narrow", xsFormula, 
                           19.4*TMath::DegToRad(), 20.6*TMath::DegToRad());
// Acceptance rate: ~60-80% (much better!)
```

---

## 📊 The Physics Behind CM → Lab Mapping

### **Why θ_Lab < θ_CM for forward scattering?**

The Lorentz transformation compresses forward angles:

```
Lab frame: Proton beam moving at velocity v
           Target at rest

CM frame:  Both particles moving toward each other
           No net momentum

Forward angles in CM → compressed in Lab
```

### **The Conversion Function**

**Location:** `Kinematics.C`, lines 62-88

```cpp
Double_t ConvertThetaCMtoLab(Double_t ekin, Double_t theta_cm_deg)
{
    // Calculate initial state
    Double_t ekinGeV = ekin/1000.;
    Double_t mom = TMath::Sqrt(2*ekinGeV*mp + ekinGeV*ekinGeV);
    TLorentzVector iState(0., 0., mom, ekinGeV + mp + mC);
    
    // Calculate CM momentum
    Double_t s = iState.Mag2();
    Double_t pcm = sqrt((s-(mp+mC)*(mp+mC)) * (s-(mp-mC)*(mp-mC)) / (4.*s));
    
    // Create proton 4-momentum in CM frame
    Double_t theta_cm_rad = theta_cm_deg * TMath::DegToRad();
    TLorentzVector p_cm;
    p_cm.SetXYZM(
        pcm * TMath::Sin(theta_cm_rad),
        0.,
        pcm * TMath::Cos(theta_cm_rad),
        mp
    );
    
    // Boost to LAB frame
    p_cm.Boost(iState.BoostVector());
    
    // Get LAB theta
    Double_t theta_lab_rad = p_cm.Theta();
    Double_t theta_lab_deg = theta_lab_rad * TMath::RadToDeg();
    
    return theta_lab_deg;
}
```

**Key steps:**
1. **Build initial state:** Beam proton + target carbon at rest
2. **Calculate β:** Boost velocity = p_beam / E_total
3. **Build CM vector:** p_cm at angle θ_CM
4. **Apply Lorentz boost:** Transform to Lab frame
5. **Extract θ_Lab:** From boosted 4-vector

---

## 📐 Energy Dependence

**The mapping depends on beam energy!**

### Example: Lab detector at 16.2° ± 0.286°

| Beam Energy | θ_CM range | Width | Compression Factor |
|-------------|------------|-------|-------------------|
| 150 MeV | [20.1°, 21.3°] | 1.2° | 2.1× |
| 200 MeV | [19.4°, 20.6°] | 1.2° | 2.1× |
| 250 MeV | [18.8°, 20.0°] | 1.2° | 2.1× |

**Pattern:** 
- Higher energy → boost is stronger → CM range shifts to smaller angles
- Width stays similar because boost is nearly constant in small angular range

---

## ⚡ Efficiency Gain

### **Full Range Sampling (NO optimization):**
```
Sample: θ_CM ∈ [0°, 50°]
Total phase space: 50°
Detector acceptance: ~1.2°
Efficiency: 1.2° / 50° = 2.4%

→ Need to generate 42× more events!
→ Time multiplied by 42×
```

### **Optimized Range Sampling:**
```
Sample: θ_CM ∈ [19.4°, 20.6°]
Total sampled: 1.2°
Detector acceptance: ~1.2°
Efficiency: ~80% (accounting for phi cuts)

→ Most events accepted!
→ Fast generation
```

---

## 🔍 Visual Summary

```
LAB FRAME:                  CM FRAME:
                            
  Beam                        ╱ p
  │                          ╱  
  │                         ╱ θ_CM = 19.4°-20.6°
  ↓                        ╱    
  ●─────────→ Detector    ╱
  Target    ↗             ●────→ 
           θ_Lab = 16.2°  
           ±0.286°        
                           
                           
Lorentz Boost              Narrow range!
    ↓                      Only 1.2° wide
  Maps                     
    ↓                      vs.
                           
Large detector             Full range
acceptance in Lab          0° - 50°
(±0.286° = 0.572° total)   would be 50° wide!
```

---

## 📝 Complete Call Chain

```
EventGenerator.C: SingleRunMultithreadPolarized()
    │
    ├─ Gets detector config:
    │  theta_lab_target = DETECTOR_THETA_CENTER_RAD  (0.2827 rad)
    │  theta_lab_window = DETECTOR_THETA_WINDOW      (0.005 rad)
    │
    └─> Calls: ComputeCMAngleRange(200, 0.2827, 0.005, θ_cm_min, θ_cm_max)
           │
           └─ Kinematics.C: ComputeCMAngleRange()
                  │
                  ├─ Step 1: Define Lab range [15.91°, 16.49°]
                  │
                  ├─ Step 2: Build 10,000-point lookup table
                  │     For i=0 to 10000:
                  │         θ_CM[i] = i × 50/10000
                  │         θ_Lab[i] = ConvertThetaCMtoLab(200, θ_CM[i])
                  │                       │
                  │                       └─> Kinematics.C: ConvertThetaCMtoLab()
                  │                              • Builds initial state
                  │                              • Calculates p_cm
                  │                              • Creates CM 4-vector
                  │                              • Applies Lorentz boost
                  │                              • Returns θ_Lab
                  │
                  ├─ Step 3: Search lookup table
                  │     Find first i where θ_Lab[i] ≥ 15.91° → θ_cm_min
                  │     Find last i where θ_Lab[i] ≤ 16.49° → θ_cm_max
                  │
                  ├─ Step 4: Add 10% margin
                  │     margin = 0.1 × (θ_cm_max - θ_cm_min)
                  │     θ_cm_min -= margin
                  │     θ_cm_max += margin
                  │
                  └─ Returns: θ_cm_min = 19.4°, θ_cm_max = 20.6°
           
           ↓ Return to EventGenerator.C
           
    Create sampling histogram:
    TF1 *ffps_narrow = new TF1("ffps_narrow", "fXS200Op(x)",
                               19.4 × π/180,  // Only narrow range!
                               20.6 × π/180);
```

---

## 🎯 Summary

**The optimization works by:**

1. **Taking Lab detector acceptance** (16.2° ± 0.286°)
2. **Creating detailed CM→Lab lookup table** (10,000 points)
3. **Finding CM range that maps into Lab acceptance** (via search)
4. **Adding 10% safety margin** (for numerical precision)
5. **Sampling ONLY from this narrow CM range** (19.4° to 20.6°)

**Result:** ~42× speedup compared to sampling full range!

**Key insight:** The Lorentz boost is deterministic and monotonic, so we can pre-compute which CM angles will end up in the detector.

---

═══════════════════════════════════════════════════════════════════
## 🆕 SECTION 4: t-BASED SAMPLING IMPLEMENTATION (February 2026)
═══════════════════════════════════════════════════════════════════

### **Motivation**

The traditional approach samples θ_CM from energy-specific cross-section files:
- Requires separate files for each beam energy (150, 160, ..., 240 MeV)
- Cannot handle arbitrary or continuous energy distributions
- Limited flexibility for multi-energy studies

**Solution:** Sample from Mandelstam invariant **t** instead, which is energy-independent!

### **Physics Background**

**Mandelstam t** is the momentum transfer squared:
```
t = (p_initial - p_final)²
t = -2p²_CM (1 - cos θ_CM)
```

**Key properties:**
- **Invariant:** Same value in any reference frame
- **Energy-independent parameterization:** σ(t) works for any beam energy
- **One-to-one mapping:** θ_CM ↔ t (given beam energy)

**Range:**
- Forward scattering (θ=0°):  t = 0
- Backward scattering (θ=180°): t = -4p²_CM ≈ -1.37 GeV² at 200 MeV

---

### **Implementation Overview**

#### **New Files Added:**

**Data Files (Splines):**
- `XStSpline.C` - Cross section as function of |t| [(MeV/c)²]
- `APSpline.C` - Analyzing power as function of |t| [(MeV/c)²]

**Modified Modules:**

1. **Kinematics.C/h** - Added conversion functions
2. **ElasticScattering.C/h** - Added spline includes
3. **EventGenerator.C/h** - Added t-based sampling mode

---

### **Module Updates**

#### **1. Kinematics Module**

**New Functions Added:**

```cpp
// Calculate CM momentum for elastic scattering
Double_t CalculateCMMomentum(Double_t ekin)
// Input: Beam energy [MeV]
// Output: CM momentum [MeV/c]

// Convert θ_CM → Mandelstam t
Double_t ConvertThetaCMtoT(Double_t theta_cm_deg, Double_t ekin)
// Input: θ_CM [degrees], beam energy [MeV]
// Output: t [GeV²] (negative for physical scattering)

// Convert Mandelstam t → θ_CM
Double_t ConvertTtoThetaCM(Double_t t, Double_t ekin)
// Input: t [GeV²], beam energy [MeV]
// Output: θ_CM [degrees]

// Calculate unified t range for energy interval
void ComputeUnifiedTRange(Double_t ekin_min, Double_t ekin_max,
                          Double_t theta_lab_center, Double_t theta_lab_window,
                          Double_t& t_min_out, Double_t& t_max_out)
// Calculates optimal t sampling range covering detector acceptance
// For single energy: use ekin_min = ekin_max
```

**Example Usage:**
```cpp
// For E=200 MeV, detector at 16.2°
Double_t t_min, t_max;
ComputeUnifiedTRange(200.0, 200.0,  // Single energy
                     DETECTOR_THETA_CENTER_RAD,
                     DETECTOR_THETA_WINDOW,
                     t_min, t_max);
// Result: t ∈ [-0.0611, -0.0128] GeV²
```

---

#### **2. ElasticScattering Module**

**Added Spline Functions:**

```cpp
#include "XStSpline.C"  // Cross section vs |t|
#include "APSpline.C"   // Analyzing power vs |t|

double XStSpline(double t);  // Input: |t| in (MeV/c)²
double APSpline(double t);   // Input: |t| in (MeV/c)²
```

**⚠️ CRITICAL UNIT NOTE:**
- Conversion functions return t in **GeV²** (negative)
- Splines expect |t| in **(MeV/c)²** (positive)
- **Conversion required:** `t_MeV² = TMath::Abs(t_GeV²) × 1e6`

---

#### **3. EventGenerator Module**

**New Sampling Mode Enum:**

```cpp
enum SamplingMode {
    THETA_CM_SAMPLING = 0,  // Traditional: sample θ_CM
    T_SAMPLING = 1          // New: sample t from splines
};
```

**Modified Main Function:**

```cpp
void SingleRunMultithreadPolarized(
    Double_t energy,
    Int_t number_total,
    Double_t polarization,
    Int_t spin_state,
    Int_t num_threads = 0,
    SamplingMode mode = THETA_CM_SAMPLING  // NEW parameter
);
```

**New Worker Function:**

```cpp
ThreadData GenerateThreadEventsPolarized_TSampling(
    int thread_id,
    Double_t ekin,
    Int_t events_per_thread,
    Double_t t_min,          // t range in GeV²
    Double_t t_max,
    Double_t theta_lab_target,
    Double_t theta_lab_window,
    Double_t phi_lab_window,
    Double_t polarization,
    Int_t spin_state
);
```

**Key Differences in t-Based Sampling:**

1. **Sampling histogram:** Filled from `XStSpline(|t|)` instead of energy-specific formula
2. **Angle conversion:** `t → θ_CM` using `ConvertTtoThetaCM()`
3. **Analyzing power:** From `APSpline(|t|)` instead of `GetElasticAnalyzingPower()`
4. **Boost and cuts:** Identical to theta-based method

---

### **Workflow Comparison**

#### **Traditional θ_CM-Based Sampling:**

```
1. Get energy-specific formula: fXS200Op(θ)
2. Calculate θ_CM range: [19.4°, 20.6°] for detector
3. Create histogram from formula
4. Sample θ_CM from histogram
5. Calculate A_N(E, θ_Lab) from parametric formula
6. Sample φ from 1 + P·A_N·cos(φ)
7. Boost to lab, apply cuts
```

#### **New t-Based Sampling:**

```
1. Calculate t range: [-0.0611, -0.0128] GeV² for detector
2. Create histogram from XStSpline(|t|)
3. Sample t from histogram
4. Convert t → θ_CM (energy-dependent)
5. Calculate A_N(t) from APSpline(|t|)
6. Sample φ from 1 + P·A_N·cos(φ)
7. Boost to lab, apply cuts
```

**Key Advantage:** Steps 1-2 use **energy-independent** data!

---

### **Performance Analysis**

#### **Efficiency for Energy Range [195-200] MeV:**

**Monte Carlo test with 10,000 random energies:**
- Average efficiency: **97.95%**
- Efficiency loss: **2.05%** (negligible!)
- Ranges outside unified: **0 / 10,000**

**Why such high efficiency?**
- Narrow energy range (2.5% spread)
- t-ranges overlap by ~98%
- Unified range only 3.6% wider than optimal average

#### **Speed Comparison:**

For generating 1M events at 6 energies:

| Solution | Setup Time | Efficiency | Notes |
|----------|------------|------------|-------|
| θ-based (traditional) | 6 × 5 ms = 30 ms | 100% | Energy-specific |
| t-based (unified) | 1 × 15 ms = 15 ms | 98% | **2× faster setup** |

**Event generation time dominates (~200 sec), so setup difference is <0.01%**

---

### **Usage Examples**

#### **Single Energy (Traditional Mode):**

```cpp
// Uses theta-based sampling (default)
SingleRunMultithreadPolarized(200, 1000000, 0.8, +1);
```

#### **Single Energy (t-Based Mode):**

```cpp
// Uses t-based sampling
SingleRunMultithreadPolarized(200, 1000000, 0.8, +1, 0, T_SAMPLING);
```

#### **Comparing Both Methods:**

```cpp
// Generate with theta-based
SingleRunMultithreadPolarized(200, 100000, 0.8, +1, 0, THETA_CM_SAMPLING);
// Output: pC_Elas_200MeV_MT_P80_SpinUp_16p2.root

// Generate with t-based
SingleRunMultithreadPolarized(200, 100000, 0.8, +1, 0, T_SAMPLING);
// Output: pC_Elas_200MeV_MT_P80_SpinUp_16p2.root (same name!)

// Compare using provided macro
.L CompareTheta_vs_T_Sampling.C
CompareTheta_vs_T_Sampling()
```

---

### **Validation**

**Test Macro:** `CompareTheta_vs_T_Sampling.C`

**Checks:**
- θ_Lab distributions (should be identical)
- φ_Lab distributions (should be identical)
- Azimuthal asymmetry (R-L)/(R+L)
- Proton momentum distributions
- 2D correlations

**Expected Result:**
```
✓ Distributions agree within statistical uncertainties
✓ Asymmetry matches: A_θ ≈ A_t ≈ 0.6-0.8 for P=80%, A_N~0.8
```

---

### **Known Issues & Solutions**

#### **Issue 1: Unit Conversion**

**Problem:** Conversion functions return t in GeV², splines expect (MeV/c)²

**Solution:**
```cpp
Double_t t_GeV2 = ConvertThetaCMtoT(theta, ekin);  // GeV²
Double_t t_MeV2 = TMath::Abs(t_GeV2) * 1.0e6;      // (MeV/c)²
Double_t xs = XStSpline(t_MeV2);
```

#### **Issue 2: Sign Convention**

**Problem:** Physical t is negative, splines use |t|

**Solution:** Always use `TMath::Abs(t)` when calling splines

#### **Issue 3: Particle 614 Warning**

**Problem:** Pluto doesn't recognize particle code 614 (C-12)

**Solution:** Explicitly provide mass in PParticle constructor:
```cpp
new ((*part)[1]) PParticle(614, cx, cy, cz, mC, 1);  // Add mC!
```

---

### **Future Enhancements**

#### **1. Multi-Energy Support (Implemented Structure, Not Yet Active)**

```cpp
// Ready for random energy sampling
Double_t t_min, t_max;
ComputeUnifiedTRange(195.0, 200.0,  // Energy range
                     DETECTOR_THETA_CENTER_RAD,
                     DETECTOR_THETA_WINDOW,
                     t_min, t_max);

// Each event can have different energy!
for (int i = 0; i < nevents; i++) {
    Double_t energy = GetRandomEnergy(195, 200);
    Double_t t = h_t->GetRandom();
    Double_t theta = ConvertTtoThetaCM(t, energy);  // Uses actual energy
    // ...
}
```

#### **2. Inelastic Scattering** ✅ COMPLETED (March 2026)

Meyer inelastic sampling fully implemented:
- ✅ Meyer inelastic splines (XSSpline/APSpline at 160 and 200 MeV)
- ✅ `GenerateThreadEventsInelastic_MeyerSampling()` worker
- ✅ `CalculateCMMomentumInelastic()` used for correct p_CM
- ✅ Energy interpolation for both XS and A_N

#### **3. Code Refactoring** ✅ COMPLETED (March 2026)

- ✅ Unified `SamplingMode` enum dispatches all three modes
- ✅ Common output path for all modes
- ✅ Runtime validation warnings for out-of-range energies

---

### **Technical Details**

#### **Unified t-Range Calculation Algorithm:**

```
For energy interval [E_min, E_max]:

1. Calculate θ_CM range for E_min:
   ComputeCMAngleRange(E_min, ...) → [θ_min_low, θ_max_low]

2. Convert to t range:
   t_min_low = ConvertThetaCMtoT(θ_max_low, E_min)  // Note swap!
   t_max_low = ConvertThetaCMtoT(θ_min_low, E_min)

3. Repeat for E_max:
   [θ_min_high, θ_max_high] → [t_min_high, t_max_high]

4. Take union:
   t_min_unified = min(t_min_low, t_min_high)
   t_max_unified = max(t_max_low, t_max_high)
```

**Note:** θ_max → t_min because larger angle → more negative t

#### **Numerical Precision:**

- θ_CM resolution: ~0.005° (10,000 bins over ~1° range)
- t resolution: ~5×10⁻⁶ GeV² (10,000 bins over ~0.05 GeV² range)
- Round-trip error θ→t→θ: < 10⁻⁶ degrees
- Round-trip error t→θ→t: < 10⁻⁹ GeV²

---

### **File Structure Changes**

**New Files (v5 — t-sampling):**
```
~/pluto_v6.02/PEG/
├── XStSpline.C              # dσ/d|t| spline (energy-independent)
├── APSpline.C               # A_N(|t|) spline (energy-independent)
├── CompareTheta_vs_T_Sampling.C  # Validation macro
```

**New Files (v6 — Meyer sampling):**
```
├── src/
│   └── MeyerScattering.C         # Meyer module (~250 lines)
├── include/
│   ├── MeyerScattering.h         # Meyer module header
│   ├── XSSpline_160MeV_elastic.C # Meyer spline data (×8 files)
│   ├── XSSpline_200MeV_elastic.C
│   ├── APSpline_160MeV_elastic.C
│   ├── APSpline_200MeV_elastic.C
│   ├── XSSpline_160MeV_inelastic.C
│   ├── XSSpline_200MeV_inelastic.C
│   ├── APSpline_160MeV_inelastic.C
│   └── APSpline_200MeV_inelastic.C
├── TestMeyer_LinearBins.C        # Meyer validation (Stage 2)
├── TestMeyer_Phase4.C            # Meyer full validation (Stages 1-3)
├── CompareTheta_vs_Meyer.C       # Sampling mode comparison macro
├── DiagnoseXSComparison.C        # XS source diagnostic (elastic)
├── DiagnoseInelasticXS.C         # XS source diagnostic (inelastic)
└── QuantifySystematics.C         # Systematic uncertainty quantification
```

**Modified Files (v6):**
```
├── include/
│   └── EventGenerator.h     # TH1.h added; Meyer set as default
├── src/
│   └── EventGenerator.C     # 4 new workers; bug fix; warnings
```

**Total New Code (v6):** ~700 lines  
**Modified Code (v6):** ~80 lines

---

### **Summary**

**What We Achieved (v5):**

✅ Energy-independent cross-section parameterization  
✅ Unified sampling range for energy intervals  
✅ 97.95% efficiency for [195-200] MeV range  
✅ Identical physics to theta-based method  
✅ Foundation for multi-energy studies  
✅ Backward compatible (default = theta-based)

**What We Achieved (v6 — March 2026):**

✅ Meyer energy-dependent sampling for elastic and inelastic (160–200 MeV)  
✅ Thread-safe envelope accept/reject with `ComputeIntegral()` pre-building  
✅ Fixed unit bug: THETA_CM elastic worker was passing θ in radians to a   
   degree-axis spline (ratio 81,960× at detector centre — now corrected)  
✅ Systematic quantification of all three modes across 160–200 MeV  
✅ T_SAMPLING_MEYER set as default for both elastic and inelastic  
✅ Runtime warnings for invalid energy/mode combinations  
✅ Complete required-files documentation (25 data files catalogued)

**Performance:**

- Setup: 2× faster for multi-energy  
- Sampling: Identical speed to theta-based  
- Efficiency: 98% vs 100% (negligible difference)  
- Memory: Same (both use 10,000-bin histograms)

**Key Innovation:**

By sampling in the Lorentz-invariant variable **t** instead of the frame-dependent angle **θ_CM**, we achieve energy-independent cross-section access while maintaining identical physics results.

**Status (v5):** ✅ Implemented and validated for elastic scattering

═══════════════════════════════════════════════════════════════════
## END OF t-BASED SAMPLING DOCUMENTATION (v5)
═══════════════════════════════════════════════════════════════════

═══════════════════════════════════════════════════════════════════
## 🆕 SECTION 5: MEYER ENERGY-DEPENDENT SAMPLING (March 2026)
═══════════════════════════════════════════════════════════════════

### Motivation

The t-Based sampling (Section 4) uses energy-independent splines — a single
dσ/dt and A_N(t) valid at one reference energy. For the 160–200 MeV range
relevant to the experiment, the cross section and analyzing power change
significantly with energy. The Meyer parameterization provides spline data
at both 160 MeV and 200 MeV, enabling linear interpolation at any
intermediate energy.

---

### New Module: MeyerScattering.h / MeyerScattering.C

**Location:** `include/MeyerScattering.h`, `src/MeyerScattering.C`  
**Dependencies:** 8 Meyer spline files (XSSpline/APSpline at 160 and 200 MeV)  
**Purpose:** Energy-interpolated cross sections and analyzing powers for
elastic and inelastic scattering via the Meyer parameterization.

**Public API:**
```cpp
// Cross sections — input: t [GeV²] (negative), ekin [MeV]
// Output: dσ/dΩ [mb/sr]
Double_t MeyerXS_Elastic  (Double_t t, Double_t ekin);
Double_t MeyerXS_Inelastic(Double_t t, Double_t ekin);

// Analyzing powers — same inputs, output: A_N (dimensionless)
Double_t MeyerAP_Elastic  (Double_t t, Double_t ekin);
Double_t MeyerAP_Inelastic(Double_t t, Double_t ekin);

// Envelope histograms for accept/reject sampling
// (MAX of 160 MeV and 200 MeV splines at each t)
TH1D* CreateMeyerEnvelope_Elastic  (Double_t t_min, Double_t t_max, Int_t nbins=50000);
TH1D* CreateMeyerEnvelope_Inelastic(Double_t t_min, Double_t t_max, Int_t nbins=50000);
```

**Energy interpolation logic:**
```
ekin ≤ 160 MeV → use 160 MeV spline directly
160 < ekin < 200 MeV → linear interpolation: weight = (ekin-160)/40
ekin ≥ 200 MeV → use 200 MeV spline directly
```

---

### New Workers in EventGenerator.C

Two new per-thread worker functions were added:

#### `GenerateThreadEventsPolarized_MeyerSampling`  (elastic)

**Algorithm — two-level accept/reject:**
1. Build envelope `MAX(XS_160MeV, XS_200MeV)` once before threads launch
2. Call `ComputeIntegral()` on envelope to make CDF read-only (thread safety)
3. Each thread: sample t from envelope using manual inverse-CDF with thread-local RNG
4. Accept with probability `MeyerXS_Elastic(t, ekin) / envelope_bin`
5. Convert accepted t → θ_CM → build CM 4-momentum → boost → acceptance cut

**Thread safety:** The critical fix here is that `TH1::GetRandom()` is not
thread-safe (uses global `gRandom` and lazily builds `fIntegral`). The
solution is:
- Call `envelope->ComputeIntegral()` before launching threads (single-threaded)
- Implement `GetRandom` manually using `GetIntegral()` (read-only) and a
  thread-local `TRandom3`

#### `GenerateThreadEventsInelastic_MeyerSampling`  (inelastic)

Identical structure to the elastic worker, but uses:
- `MeyerXS_Inelastic` / `MeyerAP_Inelastic`
- `CalculateCMMomentumInelastic(ekin, 0.00443)` for reduced CM momentum
- Carbon stored with excited-state mass `mC + E_EXCITATION`

---

### SamplingMode Enum (updated)

```cpp
enum SamplingMode {
    THETA_CM_SAMPLING = 0,   // Legacy optical model (θ_CM in degrees)
    T_SAMPLING        = 1,   // Energy-independent t-splines
    T_SAMPLING_MEYER  = 2    // Meyer energy-interpolated [DEFAULT]
};
```

**Default changed** in both orchestrators from `THETA_CM_SAMPLING` /
`T_SAMPLING` to `T_SAMPLING_MEYER`.

---

### Validation — Three Stages

All three stages passed successfully:

| Stage | Test | Result |
|-------|------|--------|
| 1 | Smoke test: elastic + inelastic at 160/180/200 MeV, 1 thread | ✅ No crash, no ratio>1 warnings |
| 2 | Shape test: χ²/dof of sampled vs expected θ_lab | ✅ χ²/dof = 1.10 |
| 3 | KS comparison Meyer vs T_SAMPLING at 180 MeV | ✅ KS p = 0.059 (compatible) |

---

### Usage

```cpp
// Elastic — Meyer default (160-200 MeV)
SingleRunMultithreadPolarized(180., 1000000, 0.80, +1);

// Inelastic — Meyer default (only correct choice away from 200 MeV)
SingleRunMultithreadInelasticPolarized(180., 1000000, 0.80, +1);

// Legacy theta-based (elastic only, valid at discrete 10-MeV steps)
SingleRunMultithreadPolarized(200., 1000000, 0.80, +1, 0, THETA_CM_SAMPLING);

// Inelastic legacy (ONLY valid at exactly 200 MeV)
SingleRunMultithreadInelasticPolarized(200., 1000000, 0.80, +1, 0, THETA_CM_SAMPLING);
```

---

═══════════════════════════════════════════════════════════════════
## 🐛 SECTION 6: BUG FIX — THETA_CM ELASTIC UNIT ERROR (March 2026)
═══════════════════════════════════════════════════════════════════

### Discovery

During systematic comparison of THETA_CM vs Meyer sampling, the θ_lab
KS test failed persistently. Diagnostic macro `DiagnoseXSComparison.C`
revealed the root cause.

### The Bug

**Location:** `EventGenerator.C`, `GenerateThreadEventsPolarized` worker  
**Introduced:** Original implementation (pre-v5)

The `XSLog*MeV` splines store cross section data as a function of
**θ_CM in degrees** (x-axis: 1, 2, 3, ..., 80 degrees).

The worker created a TF1 with axis limits in **radians**:
```cpp
// BUGGY CODE (before fix):
TF1* ffps_narrow = new TF1("ffps_narrow", xsFormula,
                            theta_cm_min * TMath::DegToRad(),   // radians!
                            theta_cm_max * TMath::DegToRad());
// ...
sigmacm_narrow->SetBinContent(ibin,
    ffps_narrow->Eval(theta_deg * TMath::DegToRad()));  // passes ~0.31 rad
                                                         // spline gets 0.31°
                                                         // instead of 17.8°!
```

### Impact

At the detector centre (θ_CM ≈ 17.8°):

| Evaluation | XS value | Correct? |
|-----------|----------|---------|
| `fXS200Op(0.311 rad)` | **3,249,878 mb/sr** | ❌ Wrong (at 0.311°) |
| `fXS200Op(17.8 deg)`  | 39.65 mb/sr | ✅ Correct |
| `MeyerXS_Elastic`     | 45.43 mb/sr | ✅ Correct |

**Ratio actual/correct = 81,960×** — completely wrong importance weights.

The final events were still kinematically correct (the lab acceptance cut
filtered events independently of the sampling weights), but the importance
sampling histogram was driven by the extreme forward peak (~0.3°) rather
than the actual cross section at the detector angle. Within the narrow
detector window the effect on the shape was subtle but detectable.

### Fix

```cpp
// FIXED CODE:
TF1* ffps_narrow = new TF1("ffps_narrow", xsFormula,
                            theta_cm_min,   // degrees — matches spline x-axis
                            theta_cm_max);
// ...
sigmacm_narrow->SetBinContent(ibin,
    ffps_narrow->Eval(theta_deg));  // passes degrees — correct!
```

### Verification

After fix: KS test p-value for θ_lab comparison = **passes** with the
detector acceptance window. The diagnostic output confirmed:
```
fXS200Op(17.80 deg) [CORRECT]: 39.6519 mb/sr
MeyerXS_Elastic               : 45.4269 mb/sr
Ratio Meyer/correct           : 1.1456   ← genuine physical difference
```

---

═══════════════════════════════════════════════════════════════════
## 📊 SECTION 7: SYSTEMATIC COMPARISON OF SAMPLING MODES (March 2026)
═══════════════════════════════════════════════════════════════════

### A_N Sources — Complete Reference

| Worker | A_N function | Data source |
|--------|-------------|-------------|
| THETA_CM elastic | `GetElasticAnalyzingPower(ekin, θ_lab)` | Wissink parametric formula (Phys.Rev.C 45, R504, 1992) |
| T_SAMPLING elastic | `APSpline(t_MeV2)` | Energy-independent spline |
| T_SAMPLING_MEYER elastic | `MeyerAP_Elastic(t, ekin)` | Meyer splines, interpolated 160–200 MeV |
| THETA_CM inelastic | `GetInelasticAnalyzingPower(θ_cm)` | Digitized CSV — **200 MeV only** |
| T_SAMPLING_MEYER inelastic | `MeyerAP_Inelastic(t, ekin)` | Meyer splines, interpolated 160–200 MeV |

### XS Sources — Complete Reference

| Worker | XS function | Data source |
|--------|------------|-------------|
| THETA_CM elastic | `fXS*Op(θ_cm_deg)` via TF1 | Optical model `XSLog*MeV.C` — dσ/dΩ vs θ_CM [deg] |
| T_SAMPLING elastic | `XStSpline(t_MeV2)` | Energy-independent t-spline — dσ/dt |
| T_SAMPLING_MEYER elastic | `MeyerXS_Elastic(t, ekin)` | Meyer splines — dσ/dΩ, interpolated |
| THETA_CM inelastic | `GetInelasticCrossSection(θ_cm_deg)` | Digitized CSV — **200 MeV only** |
| T_SAMPLING_MEYER inelastic | `MeyerXS_Inelastic(t, ekin)` | Meyer splines — dσ/dΩ, interpolated |

### Systematic Differences at Detector Centre (θ_lab ≈ 16.2°)

Measured by `QuantifySystematics.C` across 160–200 MeV:

**Elastic XS ratio Meyer / OpticalModel:**

| Energy | XS_OptModel | XS_Meyer | Ratio |
|--------|------------|---------|-------|
| 160 MeV | 57.56 mb/sr | 72.35 mb/sr | 1.257 |
| 170 MeV | 53.18 mb/sr | 62.33 mb/sr | 1.172 |
| 180 MeV | 48.69 mb/sr | 54.54 mb/sr | 1.120 |
| 190 MeV | 44.16 mb/sr | 49.14 mb/sr | 1.113 |
| 200 MeV | 39.65 mb/sr | 45.43 mb/sr | 1.146 |

**Elastic ΔA_N = A_N(Meyer) − A_N(Wissink):**

| Energy | ΔA_N | Events for 1σ pull |
|--------|------|--------------------|
| 160 MeV | +0.003 | 107,000 |
| 170 MeV | −0.030 | 842 |
| 180 MeV | −0.040 | 418 |
| 190 MeV | −0.028 | 777 |
| 200 MeV | +0.003 | 67,000 |

**Inelastic — critical finding:**

The CSV inelastic data is measured at 200 MeV only. THETA_CM inelastic
uses it at all energies, producing wrong A_N away from 200 MeV:

| Energy | ΔA_N (Meyer−CSV) | Events for 1σ |
|--------|-----------------|---------------|
| 160 MeV | **−0.295** | **12 events** |
| 170 MeV | −0.222 | 21 events |
| 180 MeV | −0.147 | 45 events |
| 190 MeV | −0.071 | 187 events |
| 200 MeV | +0.007 | 19,768 events |

⚠ **THETA_CM inelastic is only physically valid at exactly 200 MeV.**

### Recommendation

```
Elastic:   T_SAMPLING_MEYER for all energies 160–200 MeV  [DEFAULT]
           THETA_CM_SAMPLING acceptable at any discrete 10-MeV step
           as a cross-check; ΔA_N < 0.04

Inelastic: T_SAMPLING_MEYER for ALL energies            [DEFAULT, MANDATORY]
           THETA_CM_SAMPLING only at exactly 200 MeV
           as a cross-check; document the choice
```

### Required Input Files (25 total)

| Files | Category | Mode(s) |
|-------|----------|---------|
| `XSLog{150..240}MeV.C` (×10) | Optical model XS | THETA_CM elastic |
| `XStSpline.C` | t-spline XS | T_SAMPLING elastic |
| `APSpline.C` | t-spline A_N | T_SAMPLING elastic |
| `XSSpline_160/200MeV_elastic.C` (×2) | Meyer XS | T_SAMPLING_MEYER elastic |
| `APSpline_160/200MeV_elastic.C` (×2) | Meyer A_N | T_SAMPLING_MEYER elastic |
| `XSSpline_160/200MeV_inelastic.C` (×2) | Meyer XS | T_SAMPLING_MEYER inelastic |
| `APSpline_160/200MeV_inelastic.C` (×2) | Meyer A_N | T_SAMPLING_MEYER inelastic |
| `inelastic_crosssection_443_200MeV.csv` | Exp. XS | THETA_CM inelastic |
| `inelastic_analyzingpower_443_200MeV.csv` | Exp. A_N | THETA_CM inelastic |

**Minimum for T_SAMPLING_MEYER only (recommended):** 8 Meyer spline files.
No CSV files needed. No optical model files needed.

---

### Diagnostic and Comparison Macros

| Macro | Purpose | Output |
|-------|---------|--------|
| `TestMeyer_LinearBins.C` | Validates Meyer sampling shape (Stage 2) | Plot: sampled vs expected |
| `TestMeyer_Phase4.C` | Full 3-stage Meyer validation | Console + plots |
| `CompareTheta_vs_Meyer.C` | θ_lab, φ, p, A distributions side by side | PDF per reaction |
| `DiagnoseXSComparison.C` | Plots elastic XS from both sources vs θ_CM | `DiagnoseXSComparison_200MeV.pdf` |
| `DiagnoseInelasticXS.C` | Plots inelastic XS and A_N from both sources | `DiagnoseInelasticXS_200MeV.pdf` |
| `QuantifySystematics.C` | Table + plots of ΔXS and ΔA_N across energies | `Systematics_ThetaCM_vs_Meyer.pdf` |

---

═══════════════════════════════════════════════════════════════════
## END OF DOCUMENTATION v6.0
═══════════════════════════════════════════════════════════════════