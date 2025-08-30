// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Test for dimension-adaptive sparse grids and Gauss-Hermite rules

#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/detail/gauss_hermite_rules.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace boost::math::cubature;
using namespace boost::math::cubature::detail;

// Anisotropic test function - different scales per dimension
template <typename Real>
class anisotropic_function {
public:
    Real operator()(const Real* x, std::size_t) const {
        // f(x,y,z) = exp(-10*x^2 - y^2 - 0.1*z^2)
        // Much sharper in x direction than z
        return std::exp(-Real(10)*x[0]*x[0] - x[1]*x[1] - Real(0.1)*x[2]*x[2]);
    }
};

// Test function for Gaussian weight integration
template <typename Real>
class gaussian_weight_function {
public:
    Real operator()(const Real* x, std::size_t dim) const {
        // f(x) = product of (1 + x_i^2)
        // Will be integrated with exp(-||x||^2) weight
        Real prod = 1;
        for (std::size_t i = 0; i < dim; ++i) {
            prod *= (Real(1) + x[i]*x[i]);
        }
        return prod;
    }
};

void test_adaptive_sparse_grid() {
    std::cout << "\n=== Testing Dimension-Adaptive Sparse Grid ===\n\n";
    
    anisotropic_function<double> f;
    hypercube<double> box(3);
    std::fill(box.lower.begin(), box.lower.end(), -1.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // Test 1: Uniform levels (isotropic)
    std::cout << "Test 1: Isotropic grid (uniform levels)\n";
    std::vector<std::size_t> uniform_levels = {3, 3, 3};
    auto result_iso = integrate_sparse_grid_adaptive<double>(
        f, box, uniform_levels, {}, 5000);
    
    std::cout << "  Result: " << std::scientific << std::setprecision(8) 
              << result_iso.value << "\n";
    std::cout << "  Error estimate: " << result_iso.error << "\n";
    std::cout << "  Evaluations: " << result_iso.evaluations << "\n\n";
    
    // Test 2: Anisotropic initial levels
    std::cout << "Test 2: Anisotropic grid (different initial levels)\n";
    std::vector<std::size_t> aniso_levels = {5, 3, 2};  // More in x, less in z
    auto result_aniso = integrate_sparse_grid_adaptive<double>(
        f, box, aniso_levels, {}, 5000);
    
    std::cout << "  Result: " << result_aniso.value << "\n";
    std::cout << "  Error estimate: " << result_aniso.error << "\n";
    std::cout << "  Evaluations: " << result_aniso.evaluations << "\n\n";
    
    // Test 3: With dimension weights
    std::cout << "Test 3: Weighted dimensions\n";
    std::vector<std::size_t> base_levels = {2, 2, 2};
    std::vector<double> dim_weights = {10.0, 1.0, 0.1};  // Match function anisotropy
    auto result_weighted = integrate_sparse_grid_adaptive<double>(
        f, box, base_levels, dim_weights, 5000);
    
    std::cout << "  Result: " << result_weighted.value << "\n";
    std::cout << "  Error estimate: " << result_weighted.error << "\n";
    std::cout << "  Evaluations: " << result_weighted.evaluations << "\n\n";
    
    // Compare efficiency
    std::cout << "Efficiency comparison:\n";
    std::cout << "  Isotropic:  " << result_iso.evaluations << " evals\n";
    std::cout << "  Anisotropic: " << result_aniso.evaluations << " evals\n";
    std::cout << "  Weighted:    " << result_weighted.evaluations << " evals\n";
    
    if (result_aniso.evaluations < result_iso.evaluations) {
        double savings = 100.0 * (1.0 - double(result_aniso.evaluations) / 
                                        double(result_iso.evaluations));
        std::cout << "  -> Anisotropic saved " << std::fixed << std::setprecision(1) 
                  << savings << "% evaluations\n";
    }
}

void test_gauss_hermite_rules() {
    std::cout << "\n=== Testing Gauss-Hermite Rules ===\n\n";
    
    // Test 1: Standard Gauss-Hermite nodes and weights
    std::cout << "Test 1: Gauss-Hermite quadrature (n=5)\n";
    auto [nodes, weights] = gauss_hermite<double>::get_rule(5);
    
    std::cout << "  Nodes and weights:\n";
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        std::cout << "    x[" << i << "] = " << std::setw(12) << std::fixed 
                  << std::setprecision(8) << nodes[i]
                  << ", w[" << i << "] = " << weights[i] << "\n";
    }
    
    // Test weight sum (should equal sqrt(pi))
    double weight_sum = 0;
    for (auto w : weights) weight_sum += w;
    double expected = std::sqrt(M_PI);
    std::cout << "\n  Weight sum: " << weight_sum 
              << " (expected: " << expected << ")\n";
    std::cout << "  Error: " << std::scientific << std::abs(weight_sum - expected) << "\n";
    
    // Test 2: Genz-Keister nested sequence
    std::cout << "\nTest 2: Genz-Keister nested sequence\n";
    std::cout << "  Level -> Points:\n";
    for (std::size_t level = 0; level <= 5; ++level) {
        std::size_t n_points = genz_keister<double>::num_points(level);
        std::cout << "    " << level << " -> " << n_points << "\n";
    }
    
    // Test 3: Polynomial exactness
    std::cout << "\nTest 3: Polynomial exactness (x^4 * exp(-x^2))\n";
    auto [nodes5, weights5] = gauss_hermite<double>::get_rule(5);
    
    // Integrate x^4 * exp(-x^2) over (-inf, inf)
    // Exact value: 3*sqrt(pi)/4
    double integral = 0;
    for (std::size_t i = 0; i < nodes5.size(); ++i) {
        double x = nodes5[i];
        integral += std::pow(x, 4) * weights5[i];
    }
    
    double exact = 3.0 * std::sqrt(M_PI) / 4.0;
    std::cout << "  Computed: " << integral << "\n";
    std::cout << "  Exact:    " << exact << "\n";
    std::cout << "  Error:    " << std::scientific << std::abs(integral - exact) << "\n";
}

void test_sparse_grid_gaussian() {
    std::cout << "\n=== Testing Sparse Grid with Gaussian Weight ===\n\n";
    
    // Note: This is a placeholder test since the implementation is TODO
    gaussian_weight_function<double> f;
    
    std::cout << "Test: Gaussian-weighted integration (NOT YET IMPLEMENTED)\n";
    auto result = integrate_sparse_grid_gaussian<double>(f, 2, 3, false);
    
    if (result.status == status_code::dimension_error) {
        std::cout << "  Status: Not yet implemented (as expected)\n";
        std::cout << "  This will integrate f(x) * exp(-||x||^2) over R^d\n";
        std::cout << "  using Gauss-Hermite or Genz-Keister nodes\n";
    }
}

int main() {
    std::cout << "Boost.Math Cubature - Advanced Sparse Grid Tests\n";
    std::cout << "================================================\n";
    
    try {
        test_adaptive_sparse_grid();
        test_gauss_hermite_rules();
        test_sparse_grid_gaussian();
        
        std::cout << "\n================================================\n";
        std::cout << "✓ All tests completed successfully!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}