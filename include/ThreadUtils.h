// ThreadUtils.h
// Thread management utilities for parallel event generation

#ifndef THREAD_UTILS_H
#define THREAD_UTILS_H

#include <vector>

// ====================================================================
// ThreadData - Data structure for thread results
// ====================================================================
// Each thread generates events and returns them in this structure
// The main thread then combines all ThreadData objects and writes to file
// ====================================================================
struct ThreadData {
    std::vector<int> event_ids;       // Event numbering (for bookkeeping)
    std::vector<double> px, py, pz;   // Proton momentum components [GeV/c]
    std::vector<double> cx, cy, cz;   // Carbon momentum components [GeV/c]
    int count;                         // Number of events generated
};

#endif // THREAD_UTILS_H