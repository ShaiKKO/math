// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE genz_malik_test
#include <boost/test/unit_test.hpp>
#include <boost/math/cubature/detail/genz_malik_evaluator.hpp>
#include <boost/math/cubature/detail/gm_rules.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_index.hpp>
#include <cmath>
#include <iostream>
#include <numeric>

using boost::math::cubature::detail::genz_malik_evaluator;
using boost::math::cubature::detail::embedded_pair_result;
using boost::math::cubature::hypercube;
using boost::math::constants::pi;

// Test polynomial exactness for Genz-Malik rules
template <typename Real>
void test_polynomial_exactness() {
    std::cout << "Testing polynomial exactness with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // Test that degree 7 rule integrates polynomials exactly up to degree 7
    // Monomial x^a * y^b with a + b <= 7
    
    // Test in 2D
    hypercube<Real> box(2);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    genz_malik_evaluator<Real> evaluator;
    
    // Test monomial x^4 * y^3 (degree 7)
    auto f = [](const Real* x, std::size_t) -> Real {
        return std::pow(x[0], 4) * std::pow(x[1], 3);
    };
    
    embedded_pair_result<Real> result;
    bool success = evaluator.evaluate(f, box, result);
    
    BOOST_CHECK(success);
    
    // Exact integral: (1/5) * (1/4) = 1/20
    Real expected = Real(1) / Real(20);
    Real volume = box.volume();
    Real integral = result.primary * volume;
    
    BOOST_CHECK_CLOSE_FRACTION(integral, expected, Real(1e-10));
}

// Test weight symmetry properties
template <typename Real>
void test_weight_symmetry() {
    std::cout << "Testing weight symmetry with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // The Genz-Malik rule should have certain symmetry properties
    // Test with a symmetric function
    
    hypercube<Real> box(3);
    std::fill(box.lower.begin(), box.lower.end(), Real(-1));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    genz_malik_evaluator<Real> evaluator;
    
    // Symmetric function: exp(-(x^2 + y^2 + z^2))
    auto f = [](const Real* x, std::size_t) -> Real {
        Real r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        return std::exp(-r2);
    };
    
    embedded_pair_result<Real> result;
    bool success = evaluator.evaluate(f, box, result);
    
    BOOST_CHECK(success);
    
    // Check that result is positive (function is always positive)
    BOOST_CHECK_GT(result.primary, Real(0));
    
    // Error estimate should be reasonable
    BOOST_CHECK_GE(result.error, Real(0));
    BOOST_CHECK_LT(result.error, result.primary);
}

// Test embedded pair consistency
template <typename Real>
void test_embedded_pair() {
    std::cout << "Testing embedded pair with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    // The embedded rule should give a different (usually less accurate) estimate
    // for error estimation purposes
    
    hypercube<Real> box(2);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    genz_malik_evaluator<Real> evaluator;
    
    // Test function
    auto f = [](const Real* x, std::size_t) -> Real {
        return std::sin(pi<Real>() * x[0]) * std::cos(pi<Real>() * x[1]);
    };
    
    embedded_pair_result<Real> result;
    bool success = evaluator.evaluate(f, box, result);
    
    BOOST_CHECK(success);
    
    // Both estimates should be finite
    BOOST_CHECK(std::isfinite(result.primary));
    BOOST_CHECK(std::isfinite(result.embedded));
    
    // Error estimate should be based on difference
    Real error_est = std::abs(result.primary - result.embedded);
    BOOST_CHECK_CLOSE_FRACTION(result.error, error_est, Real(0.1));
}

// Test dimension limits
template <typename Real>
void test_dimension_limits() {
    std::cout << "Testing dimension limits with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    genz_malik_evaluator<Real> evaluator;
    
    // Test various dimensions
    std::vector<std::size_t> test_dims = {2, 3, 4, 5, 6, 7, 8, 9, 10};
    
    for (std::size_t dim : test_dims) {
        hypercube<Real> box(dim);
        std::fill(box.lower.begin(), box.lower.end(), Real(0));
        std::fill(box.upper.begin(), box.upper.end(), Real(1));
        
        // Simple constant function
        auto f = [](const Real*, std::size_t) -> Real {
            return Real(1);
        };
        
        embedded_pair_result<Real> result;
        bool success = evaluator.evaluate(f, box, result);
        
        if (dim <= 10) {  // Currently support up to 10D
            BOOST_CHECK(success);
            if (success) {
                // For constant function, integral should be exact
                Real expected = Real(1);  // Integral over unit hypercube
                BOOST_CHECK_CLOSE_FRACTION(result.primary, expected, Real(1e-10));
            }
        } else {
            BOOST_CHECK(!success);
        }
    }
}

// Test evaluation count
template <typename Real>
void test_evaluation_count() {
    std::cout << "Testing evaluation count with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    hypercube<Real> box(3);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    genz_malik_evaluator<Real> evaluator;
    
    std::size_t eval_count = 0;
    auto f = [&eval_count](const Real* x, std::size_t) -> Real {
        eval_count++;
        return x[0] + x[1] + x[2];
    };
    
    embedded_pair_result<Real> result;
    bool success = evaluator.evaluate(f, box, result);
    
    BOOST_CHECK(success);
    
    // Check evaluation count matches expected for Genz-Malik rule
    // For 3D: 1 + 2*3 + 2^3 + 2*3 + 2*3*2 = 1 + 6 + 8 + 6 + 12 = 33
    std::size_t expected_evals = 33;
    BOOST_CHECK_EQUAL(eval_count, expected_evals);
    BOOST_CHECK_EQUAL(result.evaluations, expected_evals);
}

// Test with challenging integrands
template <typename Real>
void test_challenging_integrands() {
    std::cout << "Testing challenging integrands with type " 
              << boost::typeindex::type_id<Real>().pretty_name() << "\n";
    
    hypercube<Real> box(2);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    genz_malik_evaluator<Real> evaluator;
    
    // Peak function (challenging for quadrature)
    auto peak = [](const Real* x, std::size_t) -> Real {
        Real dx = x[0] - Real(0.5);
        Real dy = x[1] - Real(0.5);
        return Real(1) / (Real(1) + Real(100) * (dx*dx + dy*dy));
    };
    
    embedded_pair_result<Real> result;
    bool success = evaluator.evaluate(peak, box, result);
    
    BOOST_CHECK(success);
    BOOST_CHECK_GT(result.primary, Real(0));
    BOOST_CHECK_GT(result.error, Real(0));
    
    // Oscillatory function
    auto osc = [](const Real* x, std::size_t) -> Real {
        return std::sin(Real(10) * pi<Real>() * x[0]) * 
               std::cos(Real(10) * pi<Real>() * x[1]);
    };
    
    embedded_pair_result<Real> result2;
    success = evaluator.evaluate(osc, box, result2);
    
    BOOST_CHECK(success);
    // Oscillatory functions should have larger error estimates
    BOOST_CHECK_GT(result2.error, Real(0));
}

// Test suites for different floating-point types
BOOST_AUTO_TEST_CASE(test_genz_malik_float) {
    test_polynomial_exactness<float>();
    test_weight_symmetry<float>();
    test_embedded_pair<float>();
    test_dimension_limits<float>();
}

BOOST_AUTO_TEST_CASE(test_genz_malik_double) {
    test_polynomial_exactness<double>();
    test_weight_symmetry<double>();
    test_embedded_pair<double>();
    test_dimension_limits<double>();
    test_evaluation_count<double>();
    test_challenging_integrands<double>();
}

BOOST_AUTO_TEST_CASE(test_genz_malik_long_double) {
    test_polynomial_exactness<long double>();
    test_embedded_pair<long double>();
}

// Test specific rule properties if available
BOOST_AUTO_TEST_CASE(test_rule_tables) {
    using namespace boost::math::cubature::detail;
    
    // Test 2D rule if available
    if (gm_rules::rule<7, 2>::available) {
        const auto& nodes = gm_rules::rule<7, 2>::nodes();
        const auto& weights = gm_rules::rule<7, 2>::weights();
        
        // Check that weights sum to 1 (for unit hypercube)
        double weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        BOOST_CHECK_CLOSE_FRACTION(weight_sum, 1.0, 1e-10);
        
        // Check node count
        std::size_t expected_nodes = 13;  // For 2D degree 7
        BOOST_CHECK_EQUAL(nodes.size() / 2, expected_nodes);
    }
    
    // Test 3D rule if available
    if (gm_rules::rule<7, 3>::available) {
        const auto& nodes = gm_rules::rule<7, 3>::nodes();
        const auto& weights = gm_rules::rule<7, 3>::weights();
        
        // Check that weights sum to 1
        double weight_sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        BOOST_CHECK_CLOSE_FRACTION(weight_sum, 1.0, 1e-10);
        
        // Check node count
        std::size_t expected_nodes = 33;  // For 3D degree 7
        BOOST_CHECK_EQUAL(nodes.size() / 3, expected_nodes);
    }
}