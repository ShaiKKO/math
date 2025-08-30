#include <boost/math/cubature/adaptive.hpp>
#include <iostream>
#include <cmath>

using namespace boost::math::cubature;

// Simple 2D Gaussian
double gaussian_2d(const double* x, std::size_t) {
    double r2 = x[0]*x[0] + x[1]*x[1];
    return std::exp(-r2);
}

int main() {
    hypercube<double> box(2);
    box.lower = {-2.0, -2.0};
    box.upper = {2.0, 2.0};
    
    auto result = integrate_adaptive<double>(
        gaussian_2d, box, 1e-6, 1e-6);
    
    std::cout << "Integral: " << result.value << std::endl;
    std::cout << "Error: " << result.error << std::endl;
    std::cout << "Evaluations: " << result.evaluations << std::endl;
    
    // Check against known value (should be close to pi)
    double expected = M_PI;
    double diff = std::abs(result.value - expected);
    
    if (diff < 0.05) {  // Allow 5% error
        std::cout << "✓ Result matches expected value" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Result differs from expected" << std::endl;
        return 1;
    }
}
