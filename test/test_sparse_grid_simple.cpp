// Simple test for sparse grid integration
#include <boost/math/cubature/sparse_grid.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace boost::math::cubature;

int main() {
    std::cout << "Testing Sparse Grid Integration\n\n";
    
    // Test 1: Simple 2D function
    std::cout << "Test 1: 2D product integral\n";
    auto f2d = [](const double* x, std::size_t) -> double {
        return x[0] * x[1];
    };
    
    hypercube<double> box2d(2);
    std::fill(box2d.lower.begin(), box2d.lower.end(), 0.0);
    std::fill(box2d.upper.begin(), box2d.upper.end(), 1.0);
    
    auto result2d = integrate_sparse_grid<double>(f2d, box2d, 3);
    std::cout << "  Result: " << result2d.value << "\n";
    std::cout << "  Expected: 0.25\n";
    std::cout << "  Error: " << std::abs(result2d.value - 0.25) << "\n";
    std::cout << "  Evaluations: " << result2d.evaluations << "\n\n";
    
    // Test 2: Gaussian-weighted integration
    std::cout << "Test 2: Gaussian-weighted integral\n";
    auto f_gauss = [](const double* x, std::size_t) -> double {
        // f(x) = 1 (constant), integral should be (√π)^d for exp(-||x||^2)
        return 1.0;
    };
    
    auto result_gauss = integrate_sparse_grid_gaussian<double>(
        f_gauss, 2, 3, false);
    
    double expected_gauss = M_PI;  // π for 2D
    std::cout << "  Result: " << result_gauss.value << "\n";
    std::cout << "  Expected: " << expected_gauss << "\n";
    std::cout << "  Error: " << std::abs(result_gauss.value - expected_gauss) << "\n";
    std::cout << "  Evaluations: " << result_gauss.evaluations << "\n\n";
    
    std::cout << "All tests completed!\n";
    return 0;
}