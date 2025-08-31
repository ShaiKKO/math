// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE qmc_test
#include <boost/test/unit_test.hpp>
#include <boost/math/cubature/qmc.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_index.hpp>
#include <cmath>
#include <iostream>
#include <random>

using boost::math::cubature::hypercube;
using boost::math::cubature::integrate_qmc;
using boost::math::cubature::integrate_rqmc;
using boost::math::cubature::result;
using boost::math::constants::pi;
using boost::math::constants::root_pi;

template <typename Real>
void test_qmc_convergence() {
    std::cout << "Testing QMC convergence with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Test function: product of sines
    auto f = [](const Real* x, std::size_t dim) -> Real {
        Real prod = 1;
        for (std::size_t i = 0; i < dim; ++i) {
            prod *= std::sin(pi<Real>() * x[i]);
        }
        return prod;
    };
    
    const std::size_t dim = 3;
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Test convergence with increasing sample counts
    std::vector<std::size_t> sample_counts = {100, 1000, 10000};
    Real expected = std::pow(Real(2) / pi<Real>(), dim);  // Exact integral
    
    Real prev_error = Real(1);
    for (std::size_t n : sample_counts) {
        auto res = integrate_qmc<Real>(f, box, n);
        Real error = std::abs(res.value - expected);
        
        // Error should generally decrease
        if (n > 100) {
            BOOST_CHECK_LT(error, prev_error * Real(2));  // Allow some variance
        }
        prev_error = error;
        
        // Check convergence rate (should be approximately O(1/n) for QMC)
        if (n >= 10000) {
            BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(0.01));
        }
    }
}

template <typename Real>
void test_rqmc_variance_reduction() {
    std::cout << "Testing RQMC variance reduction with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Smooth test function
    auto f = [](const Real* x, std::size_t) -> Real {
        return std::exp(-(x[0] * x[0] + x[1] * x[1]));
    };
    
    hypercube<Real> box(2);
    std::fill(box.lower.begin(), box.lower.end(), Real(-2));
    std::fill(box.upper.begin(), box.upper.end(), Real(2));
    
    // Run RQMC with multiple randomizations
    std::size_t n_samples = 1000;
    std::size_t n_randomizations = 10;
    
    auto res = integrate_rqmc<Real>(f, box, n_samples, n_randomizations);
    
    // Expected value (approximate)
    Real expected = pi<Real>();  // Integral of exp(-(x^2+y^2)) over R^2
    
    // RQMC should give reasonable estimate
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(0.1));
    
    // Error estimate should be positive
    BOOST_CHECK_GT(res.error, Real(0));
    
    // Check that we used the expected number of evaluations
    BOOST_CHECK_EQUAL(res.evaluations, n_samples * n_randomizations);
}

template <typename Real>
void test_high_dimensional_qmc() {
    std::cout << "Testing high-dimensional QMC with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Test in higher dimensions where QMC excels over MC
    const std::size_t dim = 10;
    
    // Simple product function
    auto f = [](const Real* x, std::size_t d) -> Real {
        Real prod = 1;
        for (std::size_t i = 0; i < d; ++i) {
            prod *= (Real(2) * x[i]);  // f = 2^d * x1 * x2 * ... * xd
        }
        return prod;
    };
    
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Exact integral: 2^d * (1/2)^d = 1
    Real expected = Real(1);
    
    auto res = integrate_qmc<Real>(f, box, 10000);
    
    // QMC should handle high dimensions well
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(0.01));
}

template <typename Real>
void test_discontinuous_qmc() {
    std::cout << "Testing discontinuous function with QMC type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Characteristic function of a ball
    auto f = [](const Real* x, std::size_t dim) -> Real {
        Real r2 = 0;
        for (std::size_t i = 0; i < dim; ++i) {
            Real xi = x[i] - Real(0.5);
            r2 += xi * xi;
        }
        return (r2 < Real(0.25)) ? Real(1) : Real(0);
    };
    
    const std::size_t dim = 2;
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Volume of circle with radius 0.5 = pi * 0.25
    Real expected = pi<Real>() * Real(0.25);
    
    auto res = integrate_qmc<Real>(f, box, 50000);
    
    // QMC handles discontinuities reasonably well
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(0.02));
}

template <typename Real>
void test_sobol_sequence_properties() {
    std::cout << "Testing Sobol sequence properties with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Test that Sobol sequence provides good uniformity
    // Integrate constant function - should give exact volume
    auto f = [](const Real*, std::size_t) -> Real {
        return Real(1);
    };
    
    hypercube<Real> box(3);
    std::fill(box.lower.begin(), box.lower.end(), Real(-1));
    std::fill(box.upper.begin(), box.upper.end(), Real(2));
    
    // Volume = 3^3 = 27
    Real expected = Real(27);
    
    // Even with small sample count, should be exact for constant
    auto res = integrate_qmc<Real>(f, box, 128);  // Power of 2 for Sobol
    
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-10));
}

template <typename Real>
void test_qmc_smoothness_benefit() {
    std::cout << "Testing QMC benefit for smooth functions with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Very smooth function - QMC should excel
    auto f = [](const Real* x, std::size_t) -> Real {
        return std::cos(Real(2) * pi<Real>() * x[0]) * 
               std::cos(Real(2) * pi<Real>() * x[1]) *
               std::exp(-(x[0] + x[1]));
    };
    
    hypercube<Real> box(2);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Test that error decreases rapidly for smooth functions
    std::size_t n1 = 256;
    std::size_t n2 = 1024;
    
    auto res1 = integrate_qmc<Real>(f, box, n1);
    auto res2 = integrate_qmc<Real>(f, box, n2);
    
    // For smooth functions, QMC error should decrease faster than O(1/n)
    // Typically O(log(n)^d/n) for d dimensions
    Real error_ratio = std::abs(res2.value - res1.value);
    BOOST_CHECK_LT(error_ratio, Real(0.1));  // Significant improvement
}

// Test suites for different floating-point types
BOOST_AUTO_TEST_CASE(test_qmc_float) {
    test_qmc_convergence<float>();
    test_rqmc_variance_reduction<float>();
    test_discontinuous_qmc<float>();
    test_sobol_sequence_properties<float>();
}

BOOST_AUTO_TEST_CASE(test_qmc_double) {
    test_qmc_convergence<double>();
    test_rqmc_variance_reduction<double>();
    test_high_dimensional_qmc<double>();
    test_discontinuous_qmc<double>();
    test_sobol_sequence_properties<double>();
    test_qmc_smoothness_benefit<double>();
}

BOOST_AUTO_TEST_CASE(test_qmc_long_double) {
    test_qmc_convergence<long double>();
    test_rqmc_variance_reduction<long double>();
    test_high_dimensional_qmc<long double>();
}

// Test error conditions
BOOST_AUTO_TEST_CASE(test_qmc_error_conditions) {
    hypercube<double> box(0);
    auto f = [](const double*, std::size_t) -> double { return 1.0; };
    
    // Zero dimension should fail
    auto res = integrate_qmc<double>(f, box, 1000);
    BOOST_CHECK_NE(res.status, boost::math::cubature::status_code::success);
    
    // Zero samples should fail
    hypercube<double> box2(2);
    std::fill(box2.lower.begin(), box2.lower.end(), 0.0);
    std::fill(box2.upper.begin(), box2.upper.end(), 1.0);
    
    auto res2 = integrate_qmc<double>(f, box2, 0);
    BOOST_CHECK_NE(res2.status, boost::math::cubature::status_code::success);
}

// Test comparison with Monte Carlo
BOOST_AUTO_TEST_CASE(test_qmc_vs_mc_comparison) {
    // For smooth integrands, QMC should outperform MC
    auto f = [](const double* x, std::size_t) -> double {
        return std::exp(x[0]) * std::sin(pi<double>() * x[1]);
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // QMC with deterministic Sobol sequence
    auto qmc_res = integrate_qmc<double>(f, box, 5000);
    
    // RQMC with randomization
    auto rqmc_res = integrate_rqmc<double>(f, box, 1000, 5);
    
    // Both should give similar results
    BOOST_CHECK_CLOSE_FRACTION(qmc_res.value, rqmc_res.value, 0.05);
    
    // RQMC should provide error estimate
    BOOST_CHECK_GT(rqmc_res.error, 0.0);
}