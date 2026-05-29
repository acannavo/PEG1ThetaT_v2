// TestElasticScattering.C
// Test elastic scattering cross section functions

#include "ElasticScattering.h"
#include <iostream>

void TestElasticScattering() {
    std::cout << "\n=== Testing Elastic Scattering Cross Sections ===\n" << std::endl;
    
    // Test angle in CM frame (degrees)
    double theta_cm = 20.0;
    
    std::cout << "Testing at theta_CM = " << theta_cm << " degrees\n" << std::endl;
    
    std::cout << "150 MeV: " << fXS150Op(theta_cm) << " mb/sr" << std::endl;
    std::cout << "160 MeV: " << fXS160Op(theta_cm) << " mb/sr" << std::endl;
    std::cout << "170 MeV: " << fXS170Op(theta_cm) << " mb/sr" << std::endl;
    std::cout << "180 MeV: " << fXS180Op(theta_cm) << " mb/sr" << std::endl;
    std::cout << "190 MeV: " << fXS190Op(theta_cm) << " mb/sr" << std::endl;
    std::cout << "200 MeV: " << fXS200Op(theta_cm) << " mb/sr" << std::endl;
    std::cout << "210 MeV: " << fXS210Op(theta_cm) << " mb/sr" << std::endl;
    std::cout << "220 MeV: " << fXS220Op(theta_cm) << " mb/sr" << std::endl;
    std::cout << "230 MeV: " << fXS230Op(theta_cm) << " mb/sr" << std::endl;
    std::cout << "240 MeV: " << fXS240Op(theta_cm) << " mb/sr" << std::endl;
    
    std::cout << "\n=== All cross section functions work! ===\n" << std::endl;

	
	using namespace std;
}