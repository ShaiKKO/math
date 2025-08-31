// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE simplex_test
#include <boost/test/unit_test.hpp>
#include <boost/math/cubature/simplex.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_index.hpp>
#include <cmath>
#include <iostream>
#include <numeric>

using boost::math::cubature::simplex;
using boost::math::cubature::integrate_simplex;
using boost::math::cubature::result;
using boost::math::constants::pi;

// Test basic simplex integration
template <typename Real>
void test_unit_simplex() {
    std::cout << "Testing unit simplex with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Create unit simplex in 2D: vertices at (0,0), (1,0), (0,1)
    simplex<Real> simp(2);
    simp.vertices = {
        Real(0), Real(0),  // vertex 0
        Real(1), Real(0),  // vertex 1
        Real(0), Real(1)   // vertex 2
    };
    
    // Integrate constant function = 1
    auto f = [](const Real*, std::size_t) -> Real {
        return Real(1);
    };
    
    auto res = integrate_simplex<Real>(f, simp, 1e-8, 1e-6, 10000);
    
    // Volume of unit simplex in 2D = 1/2
    Real expected = Real(0.5);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-8));
}

// Test polynomial exactness
template <typename Real>
void test_polynomial_exactness() {
    std::cout << "Testing polynomial exactness with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // 2D simplex
    simplex<Real> simp(2);
    simp.vertices = {
        Real(0), Real(0),
        Real(1), Real(0),
        Real(0), Real(1)
    };
    
    // Test linear function: f(x,y) = x + y
    auto f = [](const Real* x, std::size_t) -> Real {
        return x[0] + x[1];
    };
    
    auto res = integrate_simplex<Real>(f, simp, 1e-10, 1e-8, 10000);
    
    // Exact integral: ∫∫_S (x + y) dA where S is unit simplex
    // Using barycentric coordinates: integral = 1/2 * 2/3 = 1/3
    Real expected = Real(1) / Real(3);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-8));
}

// Test simplex in higher dimensions
template <typename Real>
void test_high_dimension_simplex() {
    std::cout << "Testing high-dimensional simplex with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // 3D simplex (tetrahedron)
    simplex<Real> simp(3);
    simp.vertices = {
        Real(0), Real(0), Real(0),  // vertex 0
        Real(1), Real(0), Real(0),  // vertex 1
        Real(0), Real(1), Real(0),  // vertex 2
        Real(0), Real(0), Real(1)   // vertex 3
    };
    
    // Constant function
    auto f = [](const Real*, std::size_t) -> Real {
        return Real(1);
    };
    
    auto res = integrate_simplex<Real>(f, simp, 1e-8, 1e-6, 10000);
    
    // Volume of unit simplex in 3D = 1/6
    Real expected = Real(1) / Real(6);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-6));
}

// Test Duffy transformation
template <typename Real>
void test_duffy_transform() {
    std::cout << "Testing Duffy transformation with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // 2D simplex with Duffy transformation
    simplex<Real> simp(2);
    simp.vertices = {
        Real(0), Real(0),
        Real(2), Real(0),  // Scaled simplex
        Real(0), Real(3)
    };
    
    // Test quadratic function
    auto f = [](const Real* x, std::size_t) -> Real {
        return x[0] * x[0] + x[1] * x[1];
    };
    
    auto res = integrate_simplex<Real>(f, simp, 1e-8, 1e-6, 10000);
    
    // Check that result is positive (function is always positive)
    BOOST_CHECK_GT(res.value, Real(0));
    
    // Check that it converged
    BOOST_CHECK_EQUAL(res.status, boost::math::cubature::status_code::success);
}

// Test barycentric coordinates
template <typename Real>
void test_barycentric_integration() {
    std::cout << "Testing barycentric integration with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Unit simplex
    simplex<Real> simp(2);
    simp.vertices = {
        Real(0), Real(0),
        Real(1), Real(0),
        Real(0), Real(1)
    };
    
    // Function that depends on barycentric coordinates
    // f(x,y) = x * (1 - x - y)
    auto f = [](const Real* x, std::size_t) -> Real {
        Real lambda1 = x[0];
        Real lambda2 = x[1];
        Real lambda0 = Real(1) - lambda1 - lambda2;
        return lambda1 * lambda0;
    };
    
    auto res = integrate_simplex<Real>(f, simp, 1e-10, 1e-8, 10000);
    
    // Exact integral using barycentric formula
    // ∫∫_S λ₁ * λ₀ dA = 1! * 1! / (2 + 2)! * 2! * volume
    //                 = 1/12 * 1/2 = 1/24
    Real expected = Real(1) / Real(24);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-6));
}

// Test general simplex (not origin-based)
template <typename Real>
void test_translated_simplex() {
    std::cout << "Testing translated simplex with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Simplex translated from origin
    simplex<Real> simp(2);
    simp.vertices = {
        Real(1), Real(1),    // vertex 0 at (1,1)
        Real(3), Real(1),    // vertex 1 at (3,1)
        Real(1), Real(4)     // vertex 2 at (1,4)
    };
    
    // Constant function
    auto f = [](const Real*, std::size_t) -> Real {
        return Real(1);
    };
    
    auto res = integrate_simplex<Real>(f, simp, 1e-10, 1e-8, 10000);
    
    // Volume = 1/2 * base * height = 1/2 * 2 * 3 = 3
    Real expected = Real(3);
    BOOST_CHECK_CLOSE_FRACTION(res.value, expected, Real(1e-8));
}

// Test degenerate simplex detection
template <typename Real>
void test_degenerate_simplex() {
    std::cout << "Testing degenerate simplex with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Degenerate simplex (all vertices on a line)
    simplex<Real> simp(2);
    simp.vertices = {
        Real(0), Real(0),
        Real(1), Real(0),
        Real(2), Real(0)  // Colinear!
    };
    
    auto f = [](const Real*, std::size_t) -> Real {
        return Real(1);
    };
    
    auto res = integrate_simplex<Real>(f, simp, 1e-8, 1e-6, 10000);
    
    // Should either fail or return zero (degenerate simplex has zero volume)
    BOOST_CHECK(res.status != boost::math::cubature::status_code::success ||
                std::abs(res.value) < Real(1e-10));
}

// Test with oscillatory function
template <typename Real>
void test_oscillatory_simplex() {
    std::cout << "Testing oscillatory function on simplex with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    simplex<Real> simp(2);
    simp.vertices = {
        Real(0), Real(0),
        Real(1), Real(0),
        Real(0), Real(1)
    };
    
    // Oscillatory function
    auto f = [](const Real* x, std::size_t) -> Real {
        return std::sin(Real(10) * pi<Real>() * x[0]) * 
               std::cos(Real(10) * pi<Real>() * x[1]);
    };
    
    auto res = integrate_simplex<Real>(f, simp, 1e-4, 1e-3, 100000);
    
    // Result should be small due to oscillations
    BOOST_CHECK_LT(std::abs(res.value), Real(0.1));
    
    // Should converge (though may need many evaluations)
    BOOST_CHECK_EQUAL(res.status, boost::math::cubature::status_code::success);
}

// Test suites for different floating-point types
BOOST_AUTO_TEST_CASE(test_simplex_float) {
    test_unit_simplex<float>();
    test_polynomial_exactness<float>();
    test_duffy_transform<float>();
}

BOOST_AUTO_TEST_CASE(test_simplex_double) {
    test_unit_simplex<double>();
    test_polynomial_exactness<double>();
    test_high_dimension_simplex<double>();
    test_duffy_transform<double>();
    test_barycentric_integration<double>();
    test_translated_simplex<double>();
    test_degenerate_simplex<double>();
    test_oscillatory_simplex<double>();
}

BOOST_AUTO_TEST_CASE(test_simplex_long_double) {
    test_unit_simplex<long double>();
    test_polynomial_exactness<long double>();
}

// Test error conditions
BOOST_AUTO_TEST_CASE(test_simplex_errors) {
    // Test with wrong number of vertices
    simplex<double> simp(2);
    simp.vertices = {0.0, 0.0, 1.0};  // Only 1.5 vertices for 2D!
    
    auto f = [](const double*, std::size_t) -> double {
        return 1.0;
    };
    
    auto res = integrate_simplex<double>(f, simp, 1e-8, 1e-6, 1000);
    BOOST_CHECK_NE(res.status, boost::math::cubature::status_code::success);
    
    // Test with zero dimension
    simplex<double> simp0(0);
    auto res0 = integrate_simplex<double>(f, simp0, 1e-8, 1e-6, 1000);
    BOOST_CHECK_NE(res0.status, boost::math::cubature::status_code::success);
}