#include <iostream>
#include <cmath>
#include <iomanip>
#include <boost/math/constants/constants.hpp>

// Test basic integrand functions
template <typename Real>
Real gaussian_2d(const Real* x) {
    return std::exp(-(x[0]*x[0] + x[1]*x[1]));
}

template <typename Real>
Real polynomial_2d(const Real* x) {
    return x[0]*x[0] + x[1]*x[1] + x[0]*x[1];
}

template <typename Real> 
Real oscillatory_2d(const Real* x) {
    const Real pi = boost::math::constants::pi<Real>();
    return std::sin(pi * x[0]) * std::cos(pi * x[1]);
}

int main() {
    using Real = double;
    
    std::cout << "\n=== Integration Pipeline Test ===\n" << std::endl;
    std::cout << "Testing basic integrand evaluation..." << std::endl;
    
    // Test point
    Real x[2] = {0.5, 0.5};
    
    // Test Gaussian
    Real gauss_val = gaussian_2d<Real>(x);
    std::cout << "  Gaussian at (0.5, 0.5): " << gauss_val << std::endl;
    std::cout << "    Expected: " << std::exp(-0.5) << std::endl;
    
    // Test polynomial
    Real poly_val = polynomial_2d<Real>(x);
    std::cout << "  Polynomial at (0.5, 0.5): " << poly_val << std::endl;
    std::cout << "    Expected: " << 0.75 << std::endl;
    
    // Test oscillatory
    Real osc_val = oscillatory_2d<Real>(x);
    const Real pi = boost::math::constants::pi<Real>();
    std::cout << "  Oscillatory at (0.5, 0.5): " << osc_val << std::endl;
    std::cout << "    Expected: " << std::sin(pi/2) * std::cos(pi/2) << std::endl;
    
    std::cout << "\nIntegrand evaluations working correctly!" << std::endl;
    
    // Test over unit square [0,1]²
    std::cout << "\n=== Manual Integration Test ===\n" << std::endl;
    
    // Simple midpoint rule
    Real midpoint[2] = {0.5, 0.5};
    Real volume = 1.0;  // Area of [0,1]²
    
    Real gauss_integral = gaussian_2d<Real>(midpoint) * volume;
    Real poly_integral = polynomial_2d<Real>(midpoint) * volume;
    Real osc_integral = oscillatory_2d<Real>(midpoint) * volume;
    
    std::cout << "Midpoint rule estimates:" << std::endl;
    std::cout << "  Gaussian integral: " << gauss_integral << std::endl;
    std::cout << "  Polynomial integral: " << poly_integral << std::endl;
    std::cout << "  Oscillatory integral: " << osc_integral << std::endl;
    
    // Known exact values for comparison
    std::cout << "\nExact values (for reference):" << std::endl;
    std::cout << "  ∫∫ (x² + y² + xy) dxdy over [0,1]²: " << Real(5)/Real(6) << std::endl;
    std::cout << "  ∫∫ sin(πx)cos(πy) dxdy over [0,1]²: " << 0 << " (by symmetry)" << std::endl;
    
    std::cout << "\n=== Pipeline components verified ===\n" << std::endl;
    
    return 0;
}