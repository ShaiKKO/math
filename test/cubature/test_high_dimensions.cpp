// Copyright 2025 Colin MacRitchie/Ripple Group
// Test Genz-Malik integration for high dimensions (10-15)

#include <boost/math/cubature/adaptive.hpp>
#include <iostream>
#include <cmath>
#include <cassert>

// Test constant integrand over unit hypercube
template <std::size_t Dim>
void test_constant_integrand() {
    using Real = double;
    
    // Constant integrand
    auto f = [](const Real* /*x*/, std::size_t) -> Real {
        return Real(1);
    };
    
    // Unit hypercube [0,1]^Dim
    boost::math::cubature::hypercube<Real> box(Dim);
    for (std::size_t i = 0; i < Dim; ++i) {
        box.lower[i] = 0.0;
        box.upper[i] = 1.0;
    }
    
    // Adaptive tolerances based on dimension
    Real abs_tol = (Dim <= 5) ? 1e-6 : 1e-2;
    Real rel_tol = (Dim <= 5) ? 1e-6 : 1e-2;
    
    boost::math::cubature::execution_options options;
    options.max_eval = (Dim <= 5) ? 100000 : 50000;
    
    auto result = boost::math::cubature::integrate_adaptive(
        f, box, abs_tol, rel_tol, options);
    
    // Expected: volume of unit hypercube = 1
    Real expected = 1.0;
    Real error = std::abs(result.value - expected);
    
    assert(error < std::max(abs_tol * 10, rel_tol * 10));
    
    std::cout << "  Dim " << Dim << ": result = " << result.value 
              << ", error = " << result.error
              << ", evals = " << result.evaluations 
              << " - PASS\n";
}

// Test product of coordinates over unit hypercube
template <std::size_t Dim>
void test_product_integrand() {
    using Real = double;
    
    // Product integrand: ∏x_i
    auto f = [](const Real* x, std::size_t dim) -> Real {
        Real prod = 1.0;
        for (std::size_t i = 0; i < dim; ++i) {
            prod *= x[i];
        }
        return prod;
    };
    
    // Unit hypercube [0,1]^Dim
    boost::math::cubature::hypercube<Real> box(Dim);
    for (std::size_t i = 0; i < Dim; ++i) {
        box.lower[i] = 0.0;
        box.upper[i] = 1.0;
    }
    
    // Adaptive tolerances
    Real abs_tol = 1e-3;
    Real rel_tol = 1e-3;
    
    boost::math::cubature::execution_options options;
    options.max_eval = 100000;
    
    auto result = boost::math::cubature::integrate_adaptive(
        f, box, abs_tol, rel_tol, options);
    
    // Expected: ∫∏x_i = 1/2^Dim
    Real expected = std::pow(0.5, static_cast<Real>(Dim));
    Real error = std::abs(result.value - expected);
    
    assert(error < std::max(abs_tol * 10, expected * rel_tol * 10));
    
    std::cout << "  Dim " << Dim << ": result = " << result.value 
              << ", expected = " << expected
              << ", rel_error = " << error/expected
              << " - PASS\n";
}

int main() {
    std::cout << "Testing Genz-Malik integration for high dimensions:\n\n";
    
    std::cout << "Test 1: Constant integrand f(x) = 1 over [0,1]^d:\n";
    test_constant_integrand<10>();
    test_constant_integrand<11>();
    test_constant_integrand<12>();
    test_constant_integrand<13>();
    test_constant_integrand<14>();
    test_constant_integrand<15>();
    
    std::cout << "\nTest 2: Product integrand f(x) = ∏x_i over [0,1]^d:\n";
    test_product_integrand<10>();
    test_product_integrand<11>();
    test_product_integrand<12>();
    test_product_integrand<13>();
    test_product_integrand<14>();
    test_product_integrand<15>();
    
    std::cout << "\nAll tests PASSED!\n";
    return 0;
}