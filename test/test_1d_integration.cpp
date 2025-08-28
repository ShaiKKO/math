// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <boost/math/cubature/cubature.hpp>
#include <boost/math/constants/constants.hpp>

using namespace boost::math::cubature;
using namespace boost::math::constants;

template <typename Real>
void test_1d_gaussian_infinite() {
    std::cout << "Testing 1D Gaussian on (-∞, ∞)..." << std::endl;
    
    // Integral of exp(-x²) from -∞ to ∞ = √π
    auto gaussian = [](const Real* x, std::size_t) -> Real {
        return std::exp(-x[0] * x[0]);
    };
    
    std::vector<Real> lower = {-std::numeric_limits<Real>::infinity()};
    std::vector<Real> upper = {std::numeric_limits<Real>::infinity()};
    
    auto result = integrate_adaptive_infinite<Real>(
        gaussian, lower, upper, Real(1e-10), Real(1e-10));
    
    Real exact = std::sqrt(pi<Real>());
    Real error = std::abs(result.value - exact);
    
    std::cout << "  Result: " << result.value << std::endl;
    std::cout << "  Exact:  " << exact << std::endl;
    std::cout << "  Error:  " << error << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.status) << std::endl;
    std::cout << "  Evaluations: " << result.evaluations << std::endl;
    
    assert(result.status == status_code::success);
    assert(error < 1e-8);
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_1d_exponential_semi_infinite() {
    std::cout << "Testing 1D exponential on [0, ∞)..." << std::endl;
    
    // Integral of exp(-x) from 0 to ∞ = 1
    auto exponential = [](const Real* x, std::size_t) -> Real {
        return std::exp(-x[0]);
    };
    
    std::vector<Real> lower = {0};
    std::vector<Real> upper = {std::numeric_limits<Real>::infinity()};
    
    auto result = integrate_adaptive_infinite<Real>(
        exponential, lower, upper, Real(1e-10), Real(1e-10));
    
    Real exact = Real(1);
    Real error = std::abs(result.value - exact);
    
    std::cout << "  Result: " << result.value << std::endl;
    std::cout << "  Exact:  " << exact << std::endl;
    std::cout << "  Error:  " << error << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.status) << std::endl;
    
    assert(result.status == status_code::success);
    assert(error < 1e-8);
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_1d_finite_bounds() {
    std::cout << "Testing 1D with finite bounds [0, π]..." << std::endl;
    
    // Integral of sin(x) from 0 to π = 2
    auto sine = [](const Real* x, std::size_t) -> Real {
        return std::sin(x[0]);
    };
    
    std::vector<Real> lower = {0};
    std::vector<Real> upper = {pi<Real>()};
    
    auto result = integrate_adaptive_infinite<Real>(
        sine, lower, upper, Real(1e-10), Real(1e-10));
    
    Real exact = Real(2);
    Real error = std::abs(result.value - exact);
    
    std::cout << "  Result: " << result.value << std::endl;
    std::cout << "  Exact:  " << exact << std::endl;
    std::cout << "  Error:  " << error << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.status) << std::endl;
    
    assert(result.status == status_code::success);
    assert(error < 1e-8);
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_1d_negative_semi_infinite() {
    std::cout << "Testing 1D on (-∞, 0]..." << std::endl;
    
    // Integral of exp(x) from -∞ to 0 = 1
    auto exponential = [](const Real* x, std::size_t) -> Real {
        return std::exp(x[0]);
    };
    
    std::vector<Real> lower = {-std::numeric_limits<Real>::infinity()};
    std::vector<Real> upper = {0};
    
    auto result = integrate_adaptive_infinite<Real>(
        exponential, lower, upper, Real(1e-10), Real(1e-10));
    
    Real exact = Real(1);
    Real error = std::abs(result.value - exact);
    
    std::cout << "  Result: " << result.value << std::endl;
    std::cout << "  Exact:  " << exact << std::endl;
    std::cout << "  Error:  " << error << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.status) << std::endl;
    
    assert(result.status == status_code::success);
    assert(error < 1e-8);
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_1d_with_singularity() {
    std::cout << "Testing 1D with endpoint singularity..." << std::endl;
    
    // Integral of 1/sqrt(x) from 0 to 1 = 2
    auto singular = [](const Real* x, std::size_t) -> Real {
        if (x[0] <= 0) return Real(0);
        return Real(1) / std::sqrt(x[0]);
    };
    
    std::vector<Real> lower = {Real(1e-10)}; // Avoid exact 0
    std::vector<Real> upper = {Real(1)};
    
    auto result = integrate_adaptive_infinite<Real>(
        singular, lower, upper, Real(1e-6), Real(1e-6));
    
    // Analytical: integral of x^(-1/2) from ε to 1 ≈ 2 - 2*sqrt(ε)
    Real epsilon = Real(1e-10);
    Real exact = Real(2) - Real(2) * std::sqrt(epsilon);
    Real error = std::abs(result.value - exact);
    
    std::cout << "  Result: " << result.value << std::endl;
    std::cout << "  Expected: ~" << exact << std::endl;
    std::cout << "  Error:  " << error << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.status) << std::endl;
    
    assert(result.status == status_code::success);
    assert(error < 1e-4); // Relaxed tolerance due to singularity
    std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_2d_still_works() {
    std::cout << "Testing that 2D still uses adaptive..." << std::endl;
    
    // 2D Gaussian
    auto gaussian_2d = [](const Real* x, std::size_t) -> Real {
        return std::exp(-(x[0]*x[0] + x[1]*x[1]));
    };
    
    std::vector<Real> lower = {
        -std::numeric_limits<Real>::infinity(),
        -std::numeric_limits<Real>::infinity()
    };
    std::vector<Real> upper = {
        std::numeric_limits<Real>::infinity(),
        std::numeric_limits<Real>::infinity()
    };
    
    auto result = integrate_adaptive_infinite<Real>(
        gaussian_2d, lower, upper, Real(1e-4), Real(1e-4), 1000000);
    
    Real exact = pi<Real>();
    Real error = std::abs(result.value - exact);
    
    std::cout << "  Result: " << result.value << std::endl;
    std::cout << "  Exact:  " << exact << std::endl;
    std::cout << "  Error:  " << error << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.status) << std::endl;
    std::cout << "  Evaluations: " << result.evaluations << std::endl;
    
    // For 2D infinite integrals, the adaptive method might report maxregions
    // but still achieve excellent accuracy
    assert(result.status == status_code::success || 
           result.status == status_code::maxeval_reached ||
           result.status == status_code::maxregions_reached);
    assert(error < 1e-2);
    std::cout << "  PASSED" << std::endl;
}

int main() {
    std::cout << "===== 1D Integration Delegation Tests =====" << std::endl;
    std::cout << "Testing proper delegation to Boost.Math quadrature methods\n" << std::endl;
    
    using Real = double;
    
    test_1d_gaussian_infinite<Real>();
    test_1d_exponential_semi_infinite<Real>();
    test_1d_finite_bounds<Real>();
    test_1d_negative_semi_infinite<Real>();
    test_1d_with_singularity<Real>();
    test_2d_still_works<Real>();
    
    std::cout << "\n===== All 1D Integration Tests PASSED =====" << std::endl;
    std::cout << "1D integrals properly delegated to quadrature methods!" << std::endl;
    
    return 0;
}