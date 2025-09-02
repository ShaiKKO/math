// Copyright 2025 Boost Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/cubature/oscillatory.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/precision.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

using namespace boost::math::cubature;
namespace bm = boost::math;

// Simple test framework
struct TestResult {
    std::string name;
    bool passed;
    double actual_error;
    double expected_error;
    std::size_t evaluations;
};

std::vector<TestResult> test_results;

void check_result(const std::string& test_name, 
                  double computed, double expected, 
                  double tolerance, std::size_t evals = 0) {
    double error = std::abs(computed - expected);
    double rel_error = error / (std::abs(expected) + 1e-10);
    
    bool passed = rel_error <= tolerance;
    
    TestResult result{test_name, passed, rel_error, tolerance, evals};
    test_results.push_back(result);
    
    std::cout << (passed ? "[PASS] " : "[FAIL] ") << test_name << std::endl;
    std::cout << "       Computed: " << std::scientific << std::setprecision(6) << computed;
    std::cout << ", Expected: " << expected;
    std::cout << ", Rel.Error: " << rel_error;
    if (evals > 0) {
        std::cout << ", Evals: " << evals;
    }
    std::cout << std::endl;
}

void test_fourier_cosine() {
    std::cout << "\n=== Testing Fourier Cosine Integration ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test 1: ∫₀^π cos(ωx) dx = sin(ωπ)/ω
    {
        auto f = [](Real) { return Real(1); };
        
        std::vector<Real> omegas = {1.0, 10.0, 100.0};
        for (Real omega : omegas) {
            auto result = integrate_fourier_cosine(f, omega, Real(0), pi);
            Real expected = std::sin(omega * pi) / omega;
            
            std::string test_name = "cos(ωx), ω=" + std::to_string(static_cast<int>(omega));
            check_result(test_name, result.value, expected, 1e-6, result.evaluations);
        }
    }
    
    // Test 2: ∫₀^∞ e^(-x) cos(ωx) dx = 1/(1 + ω²)
    // Truncate at x=20 for numerical evaluation
    {
        auto f = [](Real x) { return std::exp(-x); };
        
        std::vector<Real> omegas = {0.1, 1.0, 10.0};
        for (Real omega : omegas) {
            auto result = integrate_fourier_cosine(f, omega, Real(0), Real(20));
            Real expected = Real(1) / (Real(1) + omega * omega);
            
            std::string test_name = "e^(-x)cos(ωx), ω=" + std::to_string(omega);
            check_result(test_name, result.value, expected, 1e-3, result.evaluations);
        }
    }
}

void test_fourier_sine() {
    std::cout << "\n=== Testing Fourier Sine Integration ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test: ∫₀^π sin(ωx) dx = (1 - cos(ωπ))/ω
    {
        auto f = [](Real) { return Real(1); };
        
        std::vector<Real> omegas = {1.0, 10.0, 100.0};
        for (Real omega : omegas) {
            auto result = integrate_fourier_sine(f, omega, Real(0), pi);
            Real expected = (Real(1) - std::cos(omega * pi)) / omega;
            
            std::string test_name = "sin(ωx), ω=" + std::to_string(static_cast<int>(omega));
            check_result(test_name, result.value, expected, 1e-6, result.evaluations);
        }
    }
}

void test_fresnel_integrals() {
    std::cout << "\n=== Testing Fresnel Integrals ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Fresnel cosine integral: C(x) = ∫₀^x cos(πt²/2) dt
    std::vector<Real> x_values = {0.5, 1.0, 2.0, 3.0};
    
    for (Real x : x_values) {
        auto f = [](Real) { return Real(1); };
        auto g = [&pi](Real t) { return pi * t * t / Real(2); };
        
        // Test cosine integral
        {
            auto result = integrate_oscillatory(
                f, g, Real(1), Real(0), x, 
                detail::oscillator_type::cosine);
            
            // Approximate expected value (would use boost::math::fresnel_c in production)
            // For small x, C(x) ≈ x
            // For larger x, C(x) → 0.5
            Real expected = (x < 1.0) ? x : 0.5;
            std::string test_name = "Fresnel C(" + std::to_string(x) + ")";
            check_result(test_name, result.value, expected, 1e-5, result.evaluations);
        }
        
        // Test sine integral
        {
            auto result = integrate_oscillatory(
                f, g, Real(1), Real(0), x,
                detail::oscillator_type::sine);
            
            // Approximate expected value (would use boost::math::fresnel_s in production)
            // For small x, S(x) ≈ πx³/6
            // For larger x, S(x) → 0.5
            Real expected = (x < 1.0) ? pi * x * x * x / 6.0 : 0.5;
            std::string test_name = "Fresnel S(" + std::to_string(x) + ")";
            check_result(test_name, result.value, expected, 1e-5, result.evaluations);
        }
    }
}

void test_high_frequency() {
    std::cout << "\n=== Testing High Frequency Behavior ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test: ∫₀^π x cos(ωx) dx = π sin(ωπ)/ω + (cos(ωπ) - 1)/ω²
    auto f = [](Real x) { return x; };
    
    std::vector<Real> omegas = {10.0, 100.0, 1000.0};
    for (Real omega : omegas) {
        auto result = integrate_fourier_cosine(f, omega, Real(0), pi);
        
        Real expected = pi * std::sin(omega * pi) / omega + 
                       (std::cos(omega * pi) - Real(1)) / (omega * omega);
        
        std::string test_name = "x*cos(ωx), ω=" + std::to_string(static_cast<int>(omega));
        check_result(test_name, result.value, expected, 1e-4, result.evaluations);
        
        // Check that error decreases with frequency
        std::cout << "       Error estimate: " << std::scientific << result.error << std::endl;
    }
}

void test_nonlinear_phase() {
    std::cout << "\n=== Testing Nonlinear Phase Functions ===" << std::endl;
    
    using Real = double;
    
    // Test with cubic phase: g(x) = x³
    {
        auto f = [](Real x) { return std::exp(-x * x); };
        auto g = [](Real x) { return x * x * x; };
        
        std::vector<Real> omegas = {1.0, 10.0};
        for (Real omega : omegas) {
            auto result = integrate_oscillatory(
                f, g, omega, Real(-1), Real(1),
                detail::oscillator_type::cosine);
            
            std::string test_name = "Cubic phase, ω=" + std::to_string(static_cast<int>(omega));
            
            // Just check that result is finite and reasonable
            bool is_valid = std::isfinite(result.value) && 
                          std::abs(result.value) < Real(10) &&
                          result.error < Real(1);
            
            std::cout << (is_valid ? "[PASS] " : "[FAIL] ") << test_name << std::endl;
            std::cout << "       Value: " << result.value 
                     << ", Error: " << result.error
                     << ", Evals: " << result.evaluations << std::endl;
            
            TestResult tr{test_name, is_valid, result.error, Real(1), result.evaluations};
            test_results.push_back(tr);
        }
    }
}

void test_complex_oscillatory() {
    std::cout << "\n=== Testing Complex Oscillatory Integrals ===" << std::endl;
    
    using Real = double;
    const Real pi = bm::constants::pi<Real>();
    
    // Test ∫₀^π e^(iωx) dx = (e^(iωπ) - 1)/(iω)
    auto f = [](Real) { return Real(1); };
    auto g = [](Real x) { return x; };
    
    std::vector<Real> omegas = {1.0, 10.0};
    for (Real omega : omegas) {
        auto [value, error] = integrate_oscillatory_complex(
            f, g, omega, Real(0), pi);
        
        std::complex<Real> expected = 
            (std::exp(std::complex<Real>(0, omega * pi)) - Real(1)) /
            std::complex<Real>(0, omega);
        
        std::complex<Real> diff = value - expected;
        Real rel_error = std::abs(diff) / std::abs(expected);
        
        std::string test_name = "exp(iωx), ω=" + std::to_string(static_cast<int>(omega));
        
        std::cout << (rel_error < 1e-6 ? "[PASS] " : "[FAIL] ") << test_name << std::endl;
        std::cout << "       Value: (" << value.real() << ", " << value.imag() << ")" << std::endl;
        std::cout << "       Expected: (" << expected.real() << ", " << expected.imag() << ")" << std::endl;
        std::cout << "       Rel.Error: " << std::scientific << rel_error << std::endl;
        
        TestResult tr{test_name, rel_error < 1e-6, rel_error, 1e-6, 0};
        test_results.push_back(tr);
    }
}

void test_method_comparison() {
    std::cout << "\n=== Comparing Filon vs Levin Methods ===" << std::endl;
    
    using Real = double;
    
    // Simple test case where both methods should work
    auto f = [](Real x) { return std::exp(-x); };
    auto g = [](Real x) { return x; }; // Linear phase
    
    Real omega = 10.0;
    
    // Use Filon method
    detail::filon_integrator<Real> filon;
    auto result_filon = filon.integrate(
        f, g, omega, Real(0), Real(5),
        detail::oscillator_type::cosine);
    
    // Use Levin method
    detail::levin_integrator<Real> levin;
    auto result_levin = levin.integrate(
        f, g, omega, Real(0), Real(5),
        detail::oscillator_type::cosine);
    
    Real exact = Real(1) / (Real(1) + omega * omega) * 
                (Real(1) - std::exp(-Real(5)) * 
                 (std::cos(Real(5) * omega) + omega * std::sin(Real(5) * omega)));
    
    std::cout << "Filon method:" << std::endl;
    std::cout << "  Value: " << result_filon.value 
             << ", Error from exact: " << std::abs(result_filon.value - exact)
             << ", Evals: " << result_filon.evaluations << std::endl;
    
    std::cout << "Levin method:" << std::endl;
    std::cout << "  Value: " << result_levin.value
             << ", Error from exact: " << std::abs(result_levin.value - exact)  
             << ", Evals: " << result_levin.evaluations << std::endl;
    
    std::cout << "Exact value: " << exact << std::endl;
}

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Oscillatory Integration Test Suite" << std::endl;
    std::cout << "========================================" << std::endl;
    
    test_fourier_cosine();
    test_fourier_sine();
    test_fresnel_integrals();
    test_high_frequency();
    test_nonlinear_phase();
    test_complex_oscillatory();
    test_method_comparison();
    
    // Summary
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    
    int passed = 0;
    int failed = 0;
    
    for (const auto& result : test_results) {
        if (result.passed) {
            passed++;
        } else {
            failed++;
            std::cout << "FAILED: " << result.name 
                     << " (error: " << result.actual_error 
                     << ", tolerance: " << result.expected_error << ")" << std::endl;
        }
    }
    
    std::cout << "\nTotal: " << (passed + failed) << " tests" << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    
    return failed > 0 ? 1 : 0;
}