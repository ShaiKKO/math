// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE adaptive_integration_test
#include <boost/test/unit_test.hpp>
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_index.hpp>
#include <cmath>
#include <iostream>

using boost::math::cubature::hypercube;
using boost::math::cubature::integrate_adaptive;
using boost::math::cubature::result;
using boost::math::constants::pi;
using boost::math::constants::e;
using boost::math::constants::root_pi;

// Genz test functions for multidimensional integration
// Reference: Alan Genz, "Testing Multidimensional Integration Routines"

template <typename Real>
struct genz_functions {
    // Product peak
    static Real product_peak(const Real* x, std::size_t dim) {
        Real prod = 1;
        const Real c = Real(2);  // Peak center
        const Real w = Real(0.5); // Peak width
        for (std::size_t i = 0; i < dim; ++i) {
            prod *= 1 / (1 / (w * w) + (x[i] - c) * (x[i] - c));
        }
        return prod;
    }
    
    // Oscillatory
    static Real oscillatory(const Real* x, std::size_t dim) {
        Real sum = 0;
        const Real omega = Real(2) * pi<Real>();
        for (std::size_t i = 0; i < dim; ++i) {
            sum += x[i];
        }
        return std::cos(omega * sum);
    }
    
    // Corner peak
    static Real corner_peak(const Real* x, std::size_t dim) {
        Real sum = 0;
        for (std::size_t i = 0; i < dim; ++i) {
            sum += x[i];
        }
        return std::pow(Real(1) + sum, -Real(dim + 1));
    }
    
    // Gaussian
    static Real gaussian(const Real* x, std::size_t dim) {
        Real sum = 0;
        const Real center = Real(0.5);
        const Real sigma = Real(0.2);
        for (std::size_t i = 0; i < dim; ++i) {
            Real diff = x[i] - center;
            sum += diff * diff;
        }
        return std::exp(-sum / (Real(2) * sigma * sigma));
    }
    
    // C0 function (continuous but not differentiable)
    static Real c0_function(const Real* x, std::size_t dim) {
        Real sum = 0;
        for (std::size_t i = 0; i < dim; ++i) {
            sum += std::abs(x[i] - Real(0.5));
        }
        return std::exp(-sum);
    }
    
    // Discontinuous
    static Real discontinuous(const Real* x, std::size_t dim) {
        Real sum = 0;
        for (std::size_t i = 0; i < dim; ++i) {
            sum += x[i];
        }
        return (sum < Real(dim) / Real(2)) ? Real(1) : Real(0);
    }
};

template <typename Real>
void test_2d_polynomial() {
    std::cout << "Testing 2D polynomial with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // f(x,y) = x^2 + y^2 over [0,1]^2
    // Exact integral = 2/3
    auto f = [](const Real* x, std::size_t) -> Real {
        return x[0] * x[0] + x[1] * x[1];
    };
    
    hypercube<Real> box(2);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    Real abs_tol = Real(1e-10);
    Real rel_tol = Real(1e-8);
    
    auto res = integrate_adaptive<Real>(f, box, abs_tol, rel_tol, 10000);
    
    Real expected = Real(2) / Real(3);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-8));
    BOOST_CHECK_LE(res.error, abs_tol + rel_tol * std::abs(expected));
    BOOST_CHECK_GT(res.evaluations, std::size_t(0));
}

template <typename Real>
void test_3d_gaussian() {
    std::cout << "Testing 3D Gaussian with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Gaussian centered at origin
    auto f = [](const Real* x, std::size_t) -> Real {
        Real r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        return std::exp(-r2);
    };
    
    hypercube<Real> box(3);
    std::fill(box.lower.begin(), box.lower.end(), Real(-3));
    std::fill(box.upper.begin(), box.upper.end(), Real(3));
    
    Real abs_tol = Real(1e-6);
    Real rel_tol = Real(1e-4);
    
    auto res = integrate_adaptive<Real>(f, box, abs_tol, rel_tol, 100000);
    
    // Expected: pi^(3/2) ≈ 5.5683279968317
    Real expected = std::pow(root_pi<Real>(), 3);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-3));
}

template <typename Real>
void test_genz_oscillatory() {
    std::cout << "Testing Genz oscillatory function with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    const std::size_t dim = 3;
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    Real abs_tol = Real(1e-6);
    Real rel_tol = Real(1e-4);
    
    auto res = integrate_adaptive<Real>(
        genz_functions<Real>::oscillatory, box, abs_tol, rel_tol, 100000);
    
    // For cos(2π(x+y+z)) over [0,1]³, integral should be close to 0
    BOOST_CHECK_SMALL(res.value, Real(1e-4));
}

template <typename Real>
void test_corner_peak() {
    std::cout << "Testing corner peak function with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    const std::size_t dim = 2;
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    Real abs_tol = Real(1e-8);
    Real rel_tol = Real(1e-6);
    
    auto res = integrate_adaptive<Real>(
        genz_functions<Real>::corner_peak, box, abs_tol, rel_tol, 50000);
    
    // For 1/(1+x+y)³ over [0,1]², exact = 1/6
    Real expected = Real(1) / Real(6);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-4));
}

template <typename Real>
void test_discontinuous() {
    std::cout << "Testing discontinuous function with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    const std::size_t dim = 2;
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    Real abs_tol = Real(1e-4);
    Real rel_tol = Real(1e-3);
    
    auto res = integrate_adaptive<Real>(
        genz_functions<Real>::discontinuous, box, abs_tol, rel_tol, 100000);
    
    // Volume under x + y < 1 in [0,1]² is 1/2
    Real expected = Real(0.5);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-2));
}

template <typename Real>
void test_high_dimension() {
    std::cout << "Testing high dimension (5D) with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Simple product over [0,1]^5: ∏x_i
    auto f = [](const Real* x, std::size_t dim) -> Real {
        Real prod = 1;
        for (std::size_t i = 0; i < dim; ++i) {
            prod *= x[i];
        }
        return prod;
    };
    
    const std::size_t dim = 5;
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    Real abs_tol = Real(1e-8);
    Real rel_tol = Real(1e-6);
    
    auto res = integrate_adaptive<Real>(f, box, abs_tol, rel_tol, 100000);
    
    // Expected: (1/2)^5 = 1/32
    Real expected = std::pow(Real(0.5), dim);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-6));
}

template <typename Real>
void test_reliability_metrics() {
    std::cout << "Testing reliability metrics with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Use a simple function where we know it should converge well
    auto f = [](const Real* x, std::size_t) -> Real {
        return std::exp(-(x[0]*x[0] + x[1]*x[1]));
    };
    
    hypercube<Real> box(2);
    std::fill(box.lower.begin(), box.lower.end(), Real(-2));
    std::fill(box.upper.begin(), box.upper.end(), Real(2));
    
    Real abs_tol = Real(1e-8);
    Real rel_tol = Real(1e-6);
    
    auto res = integrate_adaptive<Real>(f, box, abs_tol, rel_tol, 50000);
    
    // Check that reliability metrics are populated
    BOOST_CHECK_GE(res.reliability.convergence_rate, Real(0));
    BOOST_CHECK_LE(res.reliability.convergence_rate, Real(10));
    BOOST_CHECK_GE(res.reliability.chi_squared, Real(0));
    BOOST_CHECK_GE(res.reliability.efficiency, Real(0));
    BOOST_CHECK_LE(res.reliability.efficiency, Real(1));
}

// Test suite for different floating-point types
BOOST_AUTO_TEST_CASE(test_adaptive_float) {
    test_2d_polynomial<float>();
    test_3d_gaussian<float>();
    test_genz_oscillatory<float>();
    test_corner_peak<float>();
    test_discontinuous<float>();
}

BOOST_AUTO_TEST_CASE(test_adaptive_double) {
    test_2d_polynomial<double>();
    test_3d_gaussian<double>();
    test_genz_oscillatory<double>();
    test_corner_peak<double>();
    test_discontinuous<double>();
    test_high_dimension<double>();
    test_reliability_metrics<double>();
}

BOOST_AUTO_TEST_CASE(test_adaptive_long_double) {
    test_2d_polynomial<long double>();
    test_3d_gaussian<long double>();
    test_genz_oscillatory<long double>();
    test_corner_peak<long double>();
}

// Test error conditions
BOOST_AUTO_TEST_CASE(test_error_conditions) {
    // Test with invalid dimension
    hypercube<double> box(0);
    auto f = [](const double*, std::size_t) -> double { return 1.0; };
    
    auto res = integrate_adaptive<double>(f, box, 1e-8, 1e-6, 1000);
    BOOST_CHECK_NE(res.status, boost::math::cubature::status_code::success);
    
    // Test with very high dimension (should fail gracefully)
    hypercube<double> box2(20);
    std::fill(box2.lower.begin(), box2.lower.end(), 0.0);
    std::fill(box2.upper.begin(), box2.upper.end(), 1.0);
    
    auto res2 = integrate_adaptive<double>(f, box2, 1e-8, 1e-6, 1000);
    BOOST_CHECK_NE(res2.status, boost::math::cubature::status_code::success);
}

// Test convergence behavior
BOOST_AUTO_TEST_CASE(test_convergence) {
    // Use a smooth function that should converge quickly
    auto f = [](const double* x, std::size_t) -> double {
        return std::exp(x[0]) * std::sin(x[1]);
    };
    
    hypercube<double> box(2);
    box.lower = {0.0, 0.0};
    box.upper = {1.0, pi<double>()};
    
    // Exact: (e - 1) * 2
    double expected = 2.0 * (e<double>() - 1.0);
    
    // Test with different tolerances
    std::vector<double> tolerances = {1e-3, 1e-5, 1e-7};
    std::size_t prev_evals = 0;
    
    for (double tol : tolerances) {
        auto res = integrate_adaptive<double>(f, box, tol, tol, 100000);
        BOOST_CHECK_CLOSE_FRACTION(res.value, expected, tol * 10);
        BOOST_CHECK_GT(res.evaluations, prev_evals);
        prev_evals = res.evaluations;
    }
}