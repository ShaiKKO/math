// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE transforms_test
#include <boost/test/unit_test.hpp>
#include <boost/math/cubature/transforms.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_index.hpp>
#include <cmath>
#include <limits>
#include <iostream>

using boost::math::cubature::rational_transform;
using boost::math::cubature::tangent_transform;
using boost::math::cubature::exponential_transform;
using boost::math::cubature::algebraic_transform;
using boost::math::cubature::duffy_transform;
using boost::math::cubature::tent;
using boost::math::constants::pi;

template <typename Real>
void test_rational_transform() {
    std::cout << "Testing rational transform with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    rational_transform<Real> transform;
    
    // Test forward transform
    // u = 0 -> x = 0
    auto result = transform.forward(Real(0));
    BOOST_CHECK_SMALL(result.first, Real(1e-10));
    
    // u = 0.5 -> x = 1
    result = transform.forward(Real(0.5));
    BOOST_CHECK_CLOSE_FRACTION(result.first, Real(1), Real(1e-10));
    BOOST_CHECK_CLOSE_FRACTION(result.second, Real(4), Real(1e-10)); // Jacobian
    
    // u = 0.75 -> x = 3
    result = transform.forward(Real(0.75));
    BOOST_CHECK_CLOSE_FRACTION(result.first, Real(3), Real(1e-10));
    
    // Test inverse transform
    Real u = transform.inverse(Real(0));
    BOOST_CHECK_SMALL(u, Real(1e-10));
    
    u = transform.inverse(Real(1));
    BOOST_CHECK_CLOSE_FRACTION(u, Real(0.5), Real(1e-10));
    
    u = transform.inverse(Real(3));
    BOOST_CHECK_CLOSE_FRACTION(u, Real(0.75), Real(1e-10));
    
    // Test round-trip consistency
    for (Real u_test : {Real(0.1), Real(0.3), Real(0.5), Real(0.7), Real(0.9)}) {
        auto [x, jac] = transform.forward(u_test);
        Real u_back = transform.inverse(x);
        BOOST_CHECK_CLOSE_FRACTION(u_back, u_test, Real(1e-10));
    }
    
    // Test boundary behavior
    result = transform.forward(Real(1) - std::numeric_limits<Real>::epsilon());
    BOOST_CHECK_GT(result.first, Real(1000));  // Should be very large
    
    result = transform.forward(Real(1));
    BOOST_CHECK(std::isinf(result.first));
    BOOST_CHECK(std::isinf(result.second));
}

template <typename Real>
void test_tangent_transform() {
    std::cout << "Testing tangent transform with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    tangent_transform<Real> transform;
    
    // Test forward transform
    // u = 0.5 -> x = 0 (center)
    auto result = transform.forward(Real(0.5));
    BOOST_CHECK_SMALL(result.first, Real(1e-10));
    BOOST_CHECK_CLOSE_FRACTION(result.second, pi<Real>(), Real(1e-10));
    
    // u = 0.25 -> x = -1
    result = transform.forward(Real(0.25));
    BOOST_CHECK_CLOSE_FRACTION(result.first, Real(-1), Real(1e-8));
    
    // u = 0.75 -> x = 1
    result = transform.forward(Real(0.75));
    BOOST_CHECK_CLOSE_FRACTION(result.first, Real(1), Real(1e-8));
    
    // Test inverse transform
    Real u = transform.inverse(Real(0));
    BOOST_CHECK_CLOSE_FRACTION(u, Real(0.5), Real(1e-10));
    
    u = transform.inverse(Real(-1));
    BOOST_CHECK_CLOSE_FRACTION(u, Real(0.25), Real(1e-8));
    
    u = transform.inverse(Real(1));
    BOOST_CHECK_CLOSE_FRACTION(u, Real(0.75), Real(1e-8));
    
    // Test round-trip consistency
    for (Real u_test : {Real(0.1), Real(0.3), Real(0.5), Real(0.7), Real(0.9)}) {
        auto [x, jac] = transform.forward(u_test);
        Real u_back = transform.inverse(x);
        BOOST_CHECK_CLOSE_FRACTION(u_back, u_test, Real(1e-8));
    }
}

template <typename Real>
void test_exponential_transform() {
    std::cout << "Testing exponential transform with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    exponential_transform<Real> transform;
    
    // Test forward transform
    // u = 0 -> x = 0
    auto result = transform.forward(Real(0));
    BOOST_CHECK_SMALL(result.first, Real(1e-10));
    BOOST_CHECK_CLOSE_FRACTION(result.second, Real(1), Real(1e-10));
    
    // u = 1 - 1/e -> x = 1
    Real u_for_1 = Real(1) - std::exp(-Real(1));
    result = transform.forward(u_for_1);
    BOOST_CHECK_CLOSE_FRACTION(result.first, Real(1), Real(1e-10));
    
    // Test inverse transform
    Real u = transform.inverse(Real(0));
    BOOST_CHECK_SMALL(u, Real(1e-10));
    
    u = transform.inverse(Real(1));
    BOOST_CHECK_CLOSE_FRACTION(u, u_for_1, Real(1e-10));
    
    u = transform.inverse(Real(2));
    BOOST_CHECK_CLOSE_FRACTION(u, Real(1) - std::exp(-Real(2)), Real(1e-10));
    
    // Test round-trip consistency
    for (Real u_test : {Real(0.1), Real(0.3), Real(0.5), Real(0.7), Real(0.9)}) {
        auto [x, jac] = transform.forward(u_test);
        Real u_back = transform.inverse(x);
        BOOST_CHECK_CLOSE_FRACTION(u_back, u_test, Real(1e-10));
    }
}

template <typename Real>
void test_algebraic_transform() {
    std::cout << "Testing algebraic transform with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Test with α = β = 1 (reduces to rational transform)
    algebraic_transform<Real> transform1(Real(1), Real(1));
    
    auto result = transform1.forward(Real(0.5));
    BOOST_CHECK_CLOSE_FRACTION(result.first, Real(1), Real(1e-10));
    
    Real u = transform1.inverse(Real(1));
    BOOST_CHECK_CLOSE_FRACTION(u, Real(0.5), Real(1e-10));
    
    // Test with α = 2, β = 1
    algebraic_transform<Real> transform2(Real(2), Real(1));
    
    result = transform2.forward(Real(0.5));
    // x = (0.5)^2 / 0.5 = 0.5
    BOOST_CHECK_CLOSE_FRACTION(result.first, Real(0.5), Real(1e-10));
    
    // Test boundary behavior
    result = transform2.forward(Real(0));
    BOOST_CHECK_SMALL(result.first, Real(1e-10));
    BOOST_CHECK_SMALL(result.second, Real(1e-10));
    
    result = transform2.forward(Real(1));
    BOOST_CHECK(std::isinf(result.first));
    BOOST_CHECK(std::isinf(result.second));
    
    // Test Newton iteration for inverse (α ≠ 1 or β ≠ 1)
    algebraic_transform<Real> transform3(Real(1.5), Real(2));
    
    // Test round-trip for general case
    for (Real u_test : {Real(0.2), Real(0.4), Real(0.6), Real(0.8)}) {
        auto [x, jac] = transform3.forward(u_test);
        Real u_back = transform3.inverse(x);
        // Newton iteration may have some error
        BOOST_CHECK_CLOSE_FRACTION(u_back, u_test, Real(1e-6));
    }
}

template <typename Real>
void test_duffy_transform() {
    std::cout << "Testing Duffy transform with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Test 1D (identity)
    {
        Real u[1] = {Real(0.5)};
        Real x[1];
        Real jac = duffy_transform<Real>::apply(u, x, 1);
        BOOST_CHECK_CLOSE_FRACTION(x[0], Real(0.5), Real(1e-10));
        BOOST_CHECK_CLOSE_FRACTION(jac, Real(1), Real(1e-10));
    }
    
    // Test 2D
    {
        Real u[2] = {Real(0.5), Real(0.5)};
        Real x[2];
        Real jac = duffy_transform<Real>::apply(u, x, 2);
        
        // x[0] = u[0] * (1 - u[1]) = 0.5 * 0.5 = 0.25
        // x[1] = u[0] * u[1] = 0.5 * 0.5 = 0.25
        BOOST_CHECK_CLOSE_FRACTION(x[0], Real(0.25), Real(1e-10));
        BOOST_CHECK_CLOSE_FRACTION(x[1], Real(0.25), Real(1e-10));
        BOOST_CHECK_CLOSE_FRACTION(jac, Real(0.5), Real(1e-10)); // u[0]
        
        // Test that points map to simplex
        BOOST_CHECK_GE(x[0], Real(0));
        BOOST_CHECK_GE(x[1], Real(0));
        BOOST_CHECK_LE(x[0] + x[1], Real(1) + Real(1e-10));
    }
    
    // Test 3D
    {
        Real u[3] = {Real(0.5), Real(0.5), Real(0.5)};
        Real x[3];
        Real jac = duffy_transform<Real>::apply(u, x, 3);
        
        // Check Jacobian
        BOOST_CHECK_CLOSE_FRACTION(jac, Real(0.125), Real(1e-10)); // u[0]^2 * u[1]
        
        // Test that points map to simplex
        Real sum = x[0] + x[1] + x[2];
        BOOST_CHECK_LE(sum, Real(1) + Real(1e-10));
    }
    
    // Test inverse
    {
        Real x[2] = {Real(0.3), Real(0.2)};
        Real u[2];
        duffy_transform<Real>::inverse(x, u, 2);
        
        // u[0] = x[0] + x[1] = 0.5
        // u[1] = x[1] / x[0] = 0.2 / 0.3 = 2/3
        BOOST_CHECK_CLOSE_FRACTION(u[0], Real(0.5), Real(1e-10));
        BOOST_CHECK_CLOSE_FRACTION(u[1], Real(2.0/3.0), Real(1e-10));
        
        // Verify round-trip
        Real x_back[2];
        duffy_transform<Real>::apply(u, x_back, 2);
        BOOST_CHECK_CLOSE_FRACTION(x_back[0], x[0], Real(1e-10));
        BOOST_CHECK_CLOSE_FRACTION(x_back[1], x[1], Real(1e-10));
    }
}

template <typename Real>
void test_tent_transform() {
    std::cout << "Testing tent transform with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Test tent transform (Baker's transformation)
    BOOST_CHECK_CLOSE_FRACTION(tent(Real(0)), Real(1), Real(1e-10));
    BOOST_CHECK_CLOSE_FRACTION(tent(Real(0.25)), Real(0.5), Real(1e-10));
    BOOST_CHECK_CLOSE_FRACTION(tent(Real(0.5)), Real(0), Real(1e-10));
    BOOST_CHECK_CLOSE_FRACTION(tent(Real(0.75)), Real(0.5), Real(1e-10));
    BOOST_CHECK_CLOSE_FRACTION(tent(Real(1)), Real(1), Real(1e-10));
    
    // Test symmetry
    for (Real u : {Real(0.1), Real(0.2), Real(0.3), Real(0.4)}) {
        Real t1 = tent(u);
        Real t2 = tent(Real(1) - u);
        BOOST_CHECK_CLOSE_FRACTION(t1, t2, Real(1e-10));
    }
    
    // Test that it's always in [0,1]
    for (Real u : {Real(0), Real(0.1), Real(0.5), Real(0.9), Real(1)}) {
        Real t = tent(u);
        BOOST_CHECK_GE(t, Real(0));
        BOOST_CHECK_LE(t, Real(1));
    }
}

// Main test suites
BOOST_AUTO_TEST_CASE(test_transforms_float) {
    test_rational_transform<float>();
    test_tangent_transform<float>();
    test_exponential_transform<float>();
    test_algebraic_transform<float>();
    test_duffy_transform<float>();
    test_tent_transform<float>();
}

BOOST_AUTO_TEST_CASE(test_transforms_double) {
    test_rational_transform<double>();
    test_tangent_transform<double>();
    test_exponential_transform<double>();
    test_algebraic_transform<double>();
    test_duffy_transform<double>();
    test_tent_transform<double>();
}

BOOST_AUTO_TEST_CASE(test_transforms_long_double) {
    test_rational_transform<long double>();
    test_tangent_transform<long double>();
    test_exponential_transform<long double>();
    test_algebraic_transform<long double>();
}

// Test transformed integrand wrapper
BOOST_AUTO_TEST_CASE(test_transformed_integrand) {
    using boost::math::cubature::make_transformed_integrand;
    
    // Integrate exp(-x) from 0 to infinity using rational transform
    auto f = [](const double* x, std::size_t) -> double {
        return std::exp(-x[0]);
    };
    
    rational_transform<double> transform;
    auto f_transformed = make_transformed_integrand<double>(
        f, transform, 1, std::vector<std::size_t>{0});
    
    // Integrate over [0,1] in transformed space
    double result = 0.0;
    std::size_t n = 100;
    for (std::size_t i = 0; i < n; ++i) {
        double u = (i + 0.5) / n;
        result += f_transformed(&u) / n;
    }
    
    // Should approximate integral of exp(-x) from 0 to infinity = 1
    BOOST_CHECK_CLOSE_FRACTION(result, 1.0, 0.01);
}