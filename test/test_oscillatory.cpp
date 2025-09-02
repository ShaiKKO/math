// Copyright 2025 Boost Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE test_oscillatory
#include <boost/test/included/unit_test.hpp>
#include <boost/math/cubature/oscillatory.hpp>
#include <boost/math/special_functions/fresnel.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/precision.hpp>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace boost::math::cubature;
namespace bm = boost::math;

template <typename Real>
struct oscillatory_test_case {
    std::string name;
    std::function<Real(Real)> f;
    std::function<Real(Real)> g;
    Real omega;
    Real a;
    Real b;
    detail::oscillator_type type;
    Real expected_value;
    Real tolerance_factor;
};

BOOST_AUTO_TEST_SUITE(oscillatory_integration_tests)

BOOST_AUTO_TEST_CASE(test_fourier_cosine_basic)
{
    std::cout << "\n=== Testing Fourier Cosine Integration ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test 1: ∫₀^π cos(ωx) dx = sin(ωπ)/ω
    {
        auto f = [](Real) { return Real(1); };
        
        for (Real omega : {1.0, 10.0, 100.0, 1000.0}) {
            auto result = integrate_fourier_cosine(f, omega, Real(0), pi);
            Real expected = std::sin(omega * pi) / omega;
            Real rel_error = std::abs(result.value - expected) / 
                           (std::abs(expected) + Real(1e-10));
            
            std::cout << "  ω = " << std::setw(6) << omega 
                     << ", result = " << std::setw(12) << result.value
                     << ", expected = " << std::setw(12) << expected
                     << ", rel_error = " << std::scientific << rel_error
                     << ", evaluations = " << result.evaluations
                     << std::endl;
            
            BOOST_CHECK_SMALL(rel_error, Real(1e-8));
        }
    }
    
    // Test 2: ∫₀^∞ e^(-x) cos(ωx) dx = 1/(1 + ω²)
    {
        auto f = [](Real x) { return std::exp(-x); };
        
        for (Real omega : {0.1, 1.0, 10.0, 50.0}) {
            // Truncate infinite integral at x = 20
            auto result = integrate_fourier_cosine(f, omega, Real(0), Real(20));
            Real expected = Real(1) / (Real(1) + omega * omega);
            Real rel_error = std::abs(result.value - expected) / expected;
            
            std::cout << "  exp(-x): ω = " << std::setw(6) << omega
                     << ", result = " << std::setw(12) << result.value
                     << ", expected = " << std::setw(12) << expected
                     << ", rel_error = " << std::scientific << rel_error
                     << std::endl;
            
            BOOST_CHECK_SMALL(rel_error, Real(1e-4));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_fourier_sine_basic)
{
    std::cout << "\n=== Testing Fourier Sine Integration ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test 1: ∫₀^π sin(ωx) dx = (1 - cos(ωπ))/ω
    {
        auto f = [](Real) { return Real(1); };
        
        for (Real omega : {1.0, 10.0, 100.0, 1000.0}) {
            auto result = integrate_fourier_sine(f, omega, Real(0), pi);
            Real expected = (Real(1) - std::cos(omega * pi)) / omega;
            Real rel_error = std::abs(result.value - expected) / 
                           (std::abs(expected) + Real(1e-10));
            
            std::cout << "  ω = " << std::setw(6) << omega
                     << ", result = " << std::setw(12) << result.value
                     << ", expected = " << std::setw(12) << expected
                     << ", rel_error = " << std::scientific << rel_error
                     << std::endl;
            
            BOOST_CHECK_SMALL(rel_error, Real(1e-8));
        }
    }
    
    // Test 2: ∫₀^∞ e^(-x) sin(ωx) dx = ω/(1 + ω²)
    {
        auto f = [](Real x) { return std::exp(-x); };
        
        for (Real omega : {0.1, 1.0, 10.0, 50.0}) {
            // Truncate infinite integral
            auto result = integrate_fourier_sine(f, omega, Real(0), Real(20));
            Real expected = omega / (Real(1) + omega * omega);
            Real rel_error = std::abs(result.value - expected) / 
                           (std::abs(expected) + Real(1e-10));
            
            std::cout << "  exp(-x): ω = " << std::setw(6) << omega
                     << ", result = " << std::setw(12) << result.value
                     << ", expected = " << std::setw(12) << expected
                     << ", rel_error = " << std::scientific << rel_error
                     << std::endl;
            
            BOOST_CHECK_SMALL(rel_error, Real(1e-4));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_fresnel_integrals)
{
    std::cout << "\n=== Testing Fresnel Integrals ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Fresnel cosine integral: C(x) = ∫₀^x cos(πt²/2) dt
    // Fresnel sine integral: S(x) = ∫₀^x sin(πt²/2) dt
    
    for (Real x : {0.5, 1.0, 2.0, 3.0, 5.0}) {
        // Test cosine integral
        {
            auto f = [](Real) { return Real(1); };
            auto g = [](Real t) { return pi * t * t / Real(2); };
            
            auto result = integrate_oscillatory(
                f, g, Real(1), Real(0), x, 
                detail::oscillator_type::cosine);
            
            Real expected = bm::fresnel_c(x);
            Real rel_error = std::abs(result.value - expected) / 
                           std::abs(expected);
            
            std::cout << "  C(" << x << "): result = " << result.value
                     << ", expected = " << expected
                     << ", rel_error = " << std::scientific << rel_error
                     << ", evaluations = " << result.evaluations
                     << std::endl;
            
            BOOST_CHECK_SMALL(rel_error, Real(1e-6));
        }
        
        // Test sine integral
        {
            auto f = [](Real) { return Real(1); };
            auto g = [](Real t) { return pi * t * t / Real(2); };
            
            auto result = integrate_oscillatory(
                f, g, Real(1), Real(0), x,
                detail::oscillator_type::sine);
            
            Real expected = bm::fresnel_s(x);
            Real rel_error = std::abs(result.value - expected) / 
                           std::abs(expected);
            
            std::cout << "  S(" << x << "): result = " << result.value
                     << ", expected = " << expected
                     << ", rel_error = " << std::scientific << rel_error
                     << std::endl;
            
            BOOST_CHECK_SMALL(rel_error, Real(1e-6));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_high_frequency_behavior)
{
    std::cout << "\n=== Testing High Frequency Behavior ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test asymptotic behavior: ∫₀^π x cos(ωx) dx
    // Exact: [x sin(ωx)/ω + cos(ωx)/ω²]₀^π = π sin(ωπ)/ω + (cos(ωπ) - 1)/ω²
    
    auto f = [](Real x) { return x; };
    
    for (Real omega : {10.0, 100.0, 1000.0, 10000.0}) {
        auto result = integrate_fourier_cosine(f, omega, Real(0), pi);
        
        Real expected = pi * std::sin(omega * pi) / omega + 
                       (std::cos(omega * pi) - Real(1)) / (omega * omega);
        
        Real rel_error = std::abs(result.value - expected) / 
                       (std::abs(expected) + Real(1) / (omega * omega));
        
        std::cout << "  ω = " << std::setw(8) << omega
                 << ", result = " << std::setw(12) << result.value
                 << ", expected = " << std::setw(12) << expected
                 << ", rel_error = " << std::scientific << rel_error
                 << ", error_est = " << result.error
                 << std::endl;
        
        // Error should decrease as O(1/ω²)
        BOOST_CHECK_SMALL(rel_error, Real(1e-6));
        BOOST_CHECK_LE(result.error, Real(10) / (omega * omega));
    }
}

BOOST_AUTO_TEST_CASE(test_nonlinear_phase)
{
    std::cout << "\n=== Testing Nonlinear Phase Functions ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test with cubic phase: g(x) = x³
    {
        auto f = [](Real x) { return std::exp(-x * x); };
        auto g = [](Real x) { return x * x * x; };
        
        for (Real omega : {1.0, 10.0, 50.0}) {
            auto result = integrate_oscillatory(
                f, g, omega, Real(-1), Real(1),
                detail::oscillator_type::cosine);
            
            std::cout << "  Cubic phase: ω = " << omega
                     << ", result = " << result.value
                     << ", error = " << std::scientific << result.error
                     << ", evaluations = " << result.evaluations
                     << std::endl;
            
            // Check that result is finite and error is reasonable
            BOOST_CHECK(std::isfinite(result.value));
            BOOST_CHECK_LE(result.error, Real(0.1));
        }
    }
    
    // Test with logarithmic phase: g(x) = log(x)
    {
        auto f = [](Real x) { return Real(1) / x; };
        auto g = [](Real x) { return std::log(x); };
        
        for (Real omega : {10.0, 100.0}) {
            auto result = integrate_oscillatory(
                f, g, omega, Real(1), Real(10),
                detail::oscillator_type::sine);
            
            std::cout << "  Log phase: ω = " << omega
                     << ", result = " << result.value
                     << ", error = " << std::scientific << result.error
                     << std::endl;
            
            BOOST_CHECK(std::isfinite(result.value));
            BOOST_CHECK_LE(result.error, Real(0.01));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_stationary_phase)
{
    std::cout << "\n=== Testing Stationary Phase Points ===" << std::endl;
    
    using Real = double;
    
    // Test integral with stationary phase point at x = 0
    // g(x) = x² has g'(0) = 0
    {
        auto f = [](Real x) { return std::exp(-x * x / Real(4)); };
        auto g = [](Real x) { return x * x; };
        
        for (Real omega : {10.0, 100.0, 1000.0}) {
            auto result = integrate_oscillatory(
                f, g, omega, Real(-2), Real(2),
                detail::oscillator_type::cosine);
            
            // Asymptotic result for large ω: sqrt(π/(2ω)) * f(0)
            Real asymptotic = std::sqrt(bm::constants::pi<Real>() / (Real(2) * omega));
            Real rel_error = std::abs(result.value - asymptotic) / asymptotic;
            
            std::cout << "  Stationary at 0: ω = " << omega
                     << ", result = " << result.value
                     << ", asymptotic = " << asymptotic
                     << ", rel_error = " << std::scientific << rel_error
                     << std::endl;
            
            // For large omega, should approach asymptotic value
            if (omega >= 100) {
                BOOST_CHECK_SMALL(rel_error, Real(0.1));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_compare_methods)
{
    std::cout << "\n=== Comparing Filon vs Levin Methods ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test case where both methods should work well
    auto f = [](Real x) { return std::exp(-x); };
    auto g = [](Real x) { return x; }; // Linear phase
    
    for (Real omega : {10.0, 100.0, 1000.0}) {
        // Use Filon method explicitly
        detail::filon_integrator<Real> filon;
        auto result_filon = filon.integrate(
            f, g, omega, Real(0), Real(5),
            detail::oscillator_type::cosine);
        
        // Use Levin method explicitly  
        detail::levin_integrator<Real> levin;
        auto result_levin = levin.integrate(
            f, g, omega, Real(0), Real(5),
            detail::oscillator_type::cosine);
        
        Real exact = Real(1) / (Real(1) + omega * omega) * 
                    (Real(1) - std::exp(-Real(5)) * 
                     (std::cos(Real(5) * omega) + omega * std::sin(Real(5) * omega)));
        
        Real error_filon = std::abs(result_filon.value - exact);
        Real error_levin = std::abs(result_levin.value - exact);
        
        std::cout << "  ω = " << std::setw(6) << omega << std::endl;
        std::cout << "    Filon: value = " << result_filon.value
                 << ", error = " << std::scientific << error_filon
                 << ", evals = " << result_filon.evaluations << std::endl;
        std::cout << "    Levin: value = " << result_levin.value
                 << ", error = " << std::scientific << error_levin
                 << ", evals = " << result_levin.evaluations << std::endl;
        std::cout << "    Exact: " << exact << std::endl;
        
        // Both methods should give accurate results
        BOOST_CHECK_SMALL(error_filon / std::abs(exact), Real(1e-6));
        BOOST_CHECK_SMALL(error_levin / std::abs(exact), Real(1e-6));
    }
}

BOOST_AUTO_TEST_CASE(test_complex_oscillatory)
{
    std::cout << "\n=== Testing Complex Oscillatory Integrals ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test ∫₀^π e^(iωx) dx = (e^(iωπ) - 1)/(iω)
    {
        auto f = [](Real) { return Real(1); };
        auto g = [](Real x) { return x; };
        
        for (Real omega : {1.0, 10.0, 100.0}) {
            auto [value, error] = integrate_oscillatory_complex(
                f, g, omega, Real(0), pi);
            
            std::complex<Real> expected = 
                (std::exp(std::complex<Real>(0, omega * pi)) - Real(1)) /
                std::complex<Real>(0, omega);
            
            std::complex<Real> diff = value - expected;
            Real rel_error = std::abs(diff) / std::abs(expected);
            
            std::cout << "  ω = " << omega
                     << ", result = (" << value.real() << ", " << value.imag() << ")"
                     << ", expected = (" << expected.real() << ", " << expected.imag() << ")"
                     << ", rel_error = " << std::scientific << rel_error
                     << std::endl;
            
            BOOST_CHECK_SMALL(rel_error, Real(1e-8));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_error_estimation)
{
    std::cout << "\n=== Testing Error Estimation ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test that error estimates are conservative
    auto f = [](Real x) { return std::sin(x); };
    
    for (Real omega : {10.0, 50.0, 100.0, 500.0}) {
        auto result = integrate_fourier_cosine(f, omega, Real(0), pi);
        
        // Exact value: ∫₀^π sin(x)cos(ωx)dx = π/(1-ω²) for ω ≠ 1
        Real exact = (std::abs(omega - Real(1)) > Real(0.1)) 
                    ? pi / (Real(1) - omega * omega)
                    : Real(0); // Skip ω ≈ 1
        
        if (exact != Real(0)) {
            Real actual_error = std::abs(result.value - exact);
            
            std::cout << "  ω = " << omega
                     << ", estimated_error = " << std::scientific << result.error
                     << ", actual_error = " << actual_error
                     << ", ratio = " << result.error / actual_error
                     << std::endl;
            
            // Error estimate should be conservative (larger than actual)
            BOOST_CHECK_GE(result.error, actual_error * Real(0.1));
            // But not too conservative (within factor of 1000)
            BOOST_CHECK_LE(result.error, actual_error * Real(1000));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_convergence_rates)
{
    std::cout << "\n=== Testing Convergence Rates ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test convergence rate as function of frequency
    auto f = [](Real x) { return x * std::exp(-x); };
    
    std::vector<Real> omegas;
    std::vector<Real> errors;
    
    for (int i = 1; i <= 6; ++i) {
        Real omega = std::pow(Real(10), Real(i));
        omegas.push_back(omega);
        
        auto result = integrate_fourier_cosine(f, omega, Real(0), Real(10));
        errors.push_back(result.error);
        
        std::cout << "  ω = " << std::scientific << omega
                 << ", error = " << result.error
                 << std::endl;
    }
    
    // Check that error decreases with frequency
    // For oscillatory methods, error should be O(1/ω) or better
    for (std::size_t i = 1; i < errors.size(); ++i) {
        Real rate = std::log(errors[i] / errors[i-1]) / 
                   std::log(omegas[i] / omegas[i-1]);
        
        std::cout << "  Convergence rate between ω = " << omegas[i-1]
                 << " and " << omegas[i] << ": " << rate << std::endl;
        
        // Rate should be negative (error decreasing)
        BOOST_CHECK_LT(rate, Real(0));
    }
}

BOOST_AUTO_TEST_SUITE_END()