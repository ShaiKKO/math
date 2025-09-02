// Copyright 2025 Boost.Math Authors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Standalone test without Boost.Test dependency

#include <boost/math/cubature/infinite.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <string>
#include <vector>
#include <chrono>
#include <cassert>

using namespace boost::math::cubature;
namespace constants = boost::math::constants;

// Simple test framework
static int test_count = 0;
static int pass_count = 0;
static int fail_count = 0;

#define CHECK_CLOSE_FRACTION(value, expected, tolerance) \
    do { \
        ++test_count; \
        double rel_error = std::abs((value) - (expected)) / std::abs(expected); \
        if (rel_error <= (tolerance)) { \
            ++pass_count; \
        } else { \
            ++fail_count; \
            std::cerr << "Test failed at " << __FILE__ << ":" << __LINE__ << "\n" \
                      << "  Expected: " << (expected) << "\n" \
                      << "  Got: " << (value) << "\n" \
                      << "  Relative error: " << rel_error << " > " << (tolerance) << "\n"; \
        } \
    } while(0)

#define CHECK(condition) \
    do { \
        ++test_count; \
        if (condition) { \
            ++pass_count; \
        } else { \
            ++fail_count; \
            std::cerr << "Test failed at " << __FILE__ << ":" << __LINE__ << "\n" \
                      << "  Condition: " << #condition << "\n"; \
        } \
    } while(0)

#define CHECK_SMALL(value, tolerance) \
    do { \
        ++test_count; \
        if (std::abs(value) <= (tolerance)) { \
            ++pass_count; \
        } else { \
            ++fail_count; \
            std::cerr << "Test failed at " << __FILE__ << ":" << __LINE__ << "\n" \
                      << "  Value " << (value) << " not small enough (> " << (tolerance) << ")\n"; \
        } \
    } while(0)

template <typename Real>
struct test_case {
    std::string name;
    std::function<Real(Real)> f;
    Real lower;
    Real upper;
    Real exact;
    Real tolerance;
};

/// \brief Test known integrals with semi-infinite domains
void test_semi_infinite_integrals()
{
    using Real = double;
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real tol = Real(1e-9);
    
    std::vector<test_case<Real>> test_cases = {
        // ∫₀^∞ e^(-x) dx = 1
        {
            "exp(-x) on [0,∞)",
            [](Real x) { return std::exp(-x); },
            Real(0), inf,
            Real(1),
            tol
        },
        
        // ∫₀^∞ e^(-x²) dx = √π/2
        {
            "exp(-x²) on [0,∞)",
            [](Real x) { return std::exp(-x * x); },
            Real(0), inf,
            constants::root_pi<Real>() / Real(2),
            tol
        },
        
        // ∫₁^∞ 1/x² dx = 1
        {
            "1/x² on [1,∞)",
            [](Real x) { return Real(1) / (x * x); },
            Real(1), inf,
            Real(1),
            tol
        },
        
        // ∫₀^∞ x e^(-x) dx = 1
        {
            "x*exp(-x) on [0,∞)",
            [](Real x) { return x * std::exp(-x); },
            Real(0), inf,
            Real(1),
            tol
        },
        
        // ∫₀^∞ 1/(1+x²) dx = π/2
        {
            "1/(1+x²) on [0,∞)",
            [](Real x) { return Real(1) / (Real(1) + x * x); },
            Real(0), inf,
            constants::pi<Real>() / Real(2),
            tol
        },
        
        // ∫₂^∞ 1/x³ dx = 1/4
        {
            "1/x³ on [2,∞)",
            [](Real x) { return Real(1) / (x * x * x); },
            Real(2), inf,
            Real(0.25),
            tol
        }
    };
    
    std::cout << "\nSemi-infinite integral tests:\n";
    std::cout << std::string(80, '-') << "\n";
    
    for (const auto& tc : test_cases) {
        // Test with automatic transform selection
        auto result_auto = integrate_ray(tc.f, tc.lower, tol);
        Real error_auto = std::abs(result_auto.value - tc.exact);
        Real rel_error_auto = error_auto / std::abs(tc.exact);
        
        std::cout << std::setw(30) << tc.name << ": "
                  << "exact = " << std::scientific << std::setprecision(10) << tc.exact
                  << ", computed = " << result_auto.value
                  << ", rel_error = " << rel_error_auto
                  << ", evals = " << result_auto.evaluations
                  << "\n";
        
        CHECK_CLOSE_FRACTION(result_auto.value, tc.exact, tc.tolerance);
        CHECK(result_auto.converged);
        
        // Test with algebraic transform
        auto result_alg = integrate_ray(tc.f, tc.lower, tol, 
                                       transform_method::algebraic);
        CHECK_CLOSE_FRACTION(result_alg.value, tc.exact, tc.tolerance * 10);
        
        // Test with exp-sinh transform
        auto result_exp = integrate_ray(tc.f, tc.lower, tol,
                                       transform_method::exp_sinh);
        CHECK_CLOSE_FRACTION(result_exp.value, tc.exact, tc.tolerance * 10);
    }
}

/// \brief Test known integrals over the real line
void test_infinite_integrals()
{
    using Real = double;
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real tol = Real(1e-9);
    
    std::vector<test_case<Real>> test_cases = {
        // ∫_{-∞}^∞ e^(-x²) dx = √π
        {
            "exp(-x²) on (-∞,∞)",
            [](Real x) { return std::exp(-x * x); },
            -inf, inf,
            constants::root_pi<Real>(),
            tol
        },
        
        // ∫_{-∞}^∞ 1/(1+x²) dx = π
        {
            "1/(1+x²) on (-∞,∞)",
            [](Real x) { return Real(1) / (Real(1) + x * x); },
            -inf, inf,
            constants::pi<Real>(),
            tol
        },
        
        // ∫_{-∞}^∞ e^(-|x|) dx = 2
        {
            "exp(-|x|) on (-∞,∞)",
            [](Real x) { return std::exp(-std::abs(x)); },
            -inf, inf,
            Real(2),
            tol
        },
        
        // ∫_{-∞}^∞ sech(x) dx = π
        {
            "sech(x) on (-∞,∞)",
            [](Real x) { 
                Real c = std::cosh(x);
                return Real(1) / c; 
            },
            -inf, inf,
            constants::pi<Real>(),
            tol * 10  // sech is harder
        },
        
        // ∫_{-∞}^∞ e^(-x²/2) dx = √(2π)
        {
            "Gaussian on (-∞,∞)",
            [](Real x) { return std::exp(-x * x / Real(2)); },
            -inf, inf,
            constants::root_two_pi<Real>(),
            tol
        }
    };
    
    std::cout << "\nInfinite integral tests:\n";
    std::cout << std::string(80, '-') << "\n";
    
    for (const auto& tc : test_cases) {
        // Test with automatic transform selection
        auto result_auto = integrate_real_line(tc.f, tol);
        Real error_auto = std::abs(result_auto.value - tc.exact);
        Real rel_error_auto = error_auto / std::abs(tc.exact);
        
        std::cout << std::setw(30) << tc.name << ": "
                  << "exact = " << std::scientific << std::setprecision(10) << tc.exact
                  << ", computed = " << result_auto.value
                  << ", rel_error = " << rel_error_auto
                  << ", evals = " << result_auto.evaluations
                  << "\n";
        
        CHECK_CLOSE_FRACTION(result_auto.value, tc.exact, tc.tolerance);
        CHECK(result_auto.converged);
        
        // Test with tanh-sinh transform
        auto result_tanh = integrate_real_line(tc.f, tol,
                                              transform_method::tanh_sinh);
        CHECK_CLOSE_FRACTION(result_tanh.value, tc.exact, tc.tolerance * 10);
        
        // Test with double algebraic transform
        auto result_alg = integrate_real_line(tc.f, tol,
                                             transform_method::algebraic);
        CHECK_CLOSE_FRACTION(result_alg.value, tc.exact, tc.tolerance * 100);
    }
}

/// \brief Test oscillatory integrals
void test_oscillatory_integrals()
{
    using Real = double;
    const Real tol = Real(1e-7);  // Oscillatory integrals are harder
    
    std::cout << "\nOscillatory integral tests:\n";
    std::cout << std::string(80, '-') << "\n";
    
    // ∫₀^∞ sin(x)/x dx = π/2
    {
        auto f = [](Real x) { 
            if (x < Real(1e-10)) return Real(1);  // Handle x=0
            return Real(1) / x; 
        };
        
        auto result = integrate_oscillatory(f, Real(1), true, tol);
        Real exact = constants::pi<Real>() / Real(2);
        
        std::cout << "sin(x)/x on [0,∞): "
                  << "exact = " << exact
                  << ", computed = " << result.value
                  << ", rel_error = " << std::abs(result.value - exact) / exact
                  << ", evals = " << result.evaluations
                  << "\n";
        
        CHECK_CLOSE_FRACTION(result.value, exact, tol * 100);
    }
    
    // ∫₀^∞ cos(x)/(1+x²) dx = π/2 * e^(-1)
    {
        auto f = [](Real x) { return Real(1) / (Real(1) + x * x); };
        
        auto result = integrate_oscillatory(f, Real(1), false, tol);
        Real exact = constants::pi<Real>() / Real(2) * std::exp(-Real(1));
        
        std::cout << "cos(x)/(1+x²) on [0,∞): "
                  << "exact = " << exact
                  << ", computed = " << result.value
                  << ", rel_error = " << std::abs(result.value - exact) / exact
                  << ", evals = " << result.evaluations
                  << "\n";
        
        CHECK_CLOSE_FRACTION(result.value, exact, tol * 100);
    }
    
    // ∫₀^∞ e^(-x) sin(2x) dx = 2/5
    {
        auto f = [](Real x) { return std::exp(-x); };
        
        auto result = integrate_oscillatory(f, Real(2), true, tol);
        Real exact = Real(2) / Real(5);
        
        std::cout << "exp(-x)*sin(2x) on [0,∞): "
                  << "exact = " << exact
                  << ", computed = " << result.value
                  << ", rel_error = " << std::abs(result.value - exact) / exact
                  << ", evals = " << result.evaluations
                  << "\n";
        
        CHECK_CLOSE_FRACTION(result.value, exact, tol * 10);
    }
}

/// \brief Test transform comparison on different function types
void test_transform_comparison()
{
    using Real = double;
    const Real tol = Real(1e-8);
    
    std::cout << "\nTransform method comparison:\n";
    std::cout << std::string(80, '-') << "\n";
    
    // Function with exponential decay
    {
        auto f = [](Real x) { return std::exp(-x) * std::cos(x); };
        Real exact = Real(0.5);  // ∫₀^∞ e^(-x) cos(x) dx = 1/2
        
        std::cout << "exp(-x)*cos(x) on [0,∞):\n";
        
        auto result_auto = integrate_ray(f, Real(0), tol);
        std::cout << "  Automatic: value = " << result_auto.value
                  << ", error = " << result_auto.error
                  << ", evals = " << result_auto.evaluations << "\n";
        
        auto result_alg = integrate_ray(f, Real(0), tol, transform_method::algebraic);
        std::cout << "  Algebraic: value = " << result_alg.value
                  << ", error = " << result_alg.error
                  << ", evals = " << result_alg.evaluations << "\n";
        
        auto result_exp = integrate_ray(f, Real(0), tol, transform_method::exp_sinh);
        std::cout << "  Exp-sinh:  value = " << result_exp.value
                  << ", error = " << result_exp.error
                  << ", evals = " << result_exp.evaluations << "\n";
        
        CHECK_CLOSE_FRACTION(result_auto.value, exact, tol);
        CHECK_CLOSE_FRACTION(result_exp.value, exact, tol);
    }
    
    // Function with power-law decay
    {
        auto f = [](Real x) { return Real(1) / std::pow(Real(1) + x, Real(3)); };
        Real exact = Real(0.5);  // ∫₀^∞ 1/(1+x)³ dx = 1/2
        
        std::cout << "1/(1+x)³ on [0,∞):\n";
        
        auto result_auto = integrate_ray(f, Real(0), tol);
        std::cout << "  Automatic: value = " << result_auto.value
                  << ", error = " << result_auto.error
                  << ", evals = " << result_auto.evaluations << "\n";
        
        auto result_alg = integrate_ray(f, Real(0), tol, transform_method::algebraic);
        std::cout << "  Algebraic: value = " << result_alg.value
                  << ", error = " << result_alg.error
                  << ", evals = " << result_alg.evaluations << "\n";
        
        auto result_exp = integrate_ray(f, Real(0), tol, transform_method::exp_sinh);
        std::cout << "  Exp-sinh:  value = " << result_exp.value
                  << ", error = " << result_exp.error
                  << ", evals = " << result_exp.evaluations << "\n";
        
        CHECK_CLOSE_FRACTION(result_auto.value, exact, tol);
        CHECK_CLOSE_FRACTION(result_alg.value, exact, tol);
    }
}

/// \brief Test edge cases and error handling
void test_edge_cases()
{
    using Real = double;
    const Real inf = std::numeric_limits<Real>::infinity();
    const Real tol = Real(1e-8);
    
    // Test zero function
    {
        auto f = [](Real) { return Real(0); };
        auto result = integrate_ray(f, Real(0), tol);
        CHECK_SMALL(result.value, tol);
    }
    
    // Test constant function (should not converge on infinite domain)
    {
        auto f = [](Real) { return Real(1); };
        auto result = integrate_ray(f, Real(0), tol);
        CHECK(!result.converged);
    }
    
    // Test function with singularity at lower bound
    {
        auto f = [](Real x) { 
            if (x < Real(1e-10)) return Real(0);
            return Real(1) / std::sqrt(x) * std::exp(-x); 
        };
        // ∫₀^∞ x^(-1/2) e^(-x) dx = √π
        auto result = integrate_ray(f, Real(0), tol * 100);
        CHECK_CLOSE_FRACTION(result.value, constants::root_pi<Real>(), Real(1e-4));
    }
    
    // Test mixed infinite bounds
    {
        auto f = [](Real x) { return std::exp(-x * x); };
        
        // ∫₋₁^∞ e^(-x²) dx
        auto result = integrate_infinite(f, Real(-1), inf, tol);
        Real exact = constants::root_pi<Real>() / Real(2) * (Real(1) + boost::math::erf(Real(1)));
        CHECK_CLOSE_FRACTION(result.value, exact, tol * 10);
    }
}

/// \brief Test convergence rates for different transforms
void test_convergence_rates()
{
    using Real = double;
    
    std::cout << "\nConvergence rate analysis:\n";
    std::cout << std::string(80, '-') << "\n";
    
    // Test function: ∫₀^∞ e^(-x²) dx = √π/2
    auto f = [](Real x) { return std::exp(-x * x); };
    Real exact = constants::root_pi<Real>() / Real(2);
    
    std::vector<Real> tolerances = {1e-4, 1e-6, 1e-8, 1e-10};
    
    std::cout << "Function: exp(-x²) on [0,∞)\n";
    std::cout << "Tolerance    Algebraic Error    Exp-Sinh Error    Tanh-Sinh Error\n";
    
    for (Real tol : tolerances) {
        auto result_alg = integrate_ray(f, Real(0), tol, transform_method::algebraic);
        auto result_exp = integrate_ray(f, Real(0), tol, transform_method::exp_sinh);
        
        // For tanh-sinh, we need to use the real line version
        auto g = [f](Real x) { 
            if (x < Real(0)) return Real(0);
            return f(x);
        };
        auto result_tanh = integrate_real_line(g, tol, transform_method::tanh_sinh);
        // Adjust for half-line
        result_tanh.value /= Real(2);
        
        std::cout << std::scientific << std::setprecision(2)
                  << tol << "    "
                  << std::abs(result_alg.value - exact) << " (" 
                  << result_alg.evaluations << ")    "
                  << std::abs(result_exp.value - exact) << " ("
                  << result_exp.evaluations << ")    "
                  << std::abs(result_tanh.value - exact) << " ("
                  << result_tanh.evaluations << ")\n";
    }
}

/// \brief Performance benchmark
void test_performance_benchmark()
{
    using Real = double;
    const Real tol = Real(1e-10);
    
    std::cout << "\nPerformance benchmark:\n";
    std::cout << std::string(80, '-') << "\n";
    
    // Benchmark different integral types
    struct benchmark_case {
        std::string name;
        std::function<Real(Real)> f;
        Real exact;
        transform_method best_method;
    };
    
    std::vector<benchmark_case> cases = {
        {
            "Exponential decay",
            [](Real x) { return std::exp(-x); },
            Real(1),
            transform_method::exp_sinh
        },
        {
            "Gaussian decay",
            [](Real x) { return std::exp(-x * x); },
            constants::root_pi<Real>() / Real(2),
            transform_method::tanh_sinh
        },
        {
            "Power law decay",
            [](Real x) { return Real(1) / (Real(1) + x * x); },
            constants::pi<Real>() / Real(2),
            transform_method::algebraic
        }
    };
    
    for (const auto& bc : cases) {
        std::cout << bc.name << ":\n";
        
        auto start = std::chrono::high_resolution_clock::now();
        auto result = integrate_ray(bc.f, Real(0), tol, bc.best_method);
        auto end = std::chrono::high_resolution_clock::now();
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        std::cout << "  Result: " << result.value
                  << " (exact: " << bc.exact << ")\n"
                  << "  Error: " << std::scientific << result.error << "\n"
                  << "  Evaluations: " << result.evaluations << "\n"
                  << "  Time: " << duration.count() << " μs\n"
                  << "  Speed: " << result.evaluations * 1000000.0 / duration.count() 
                  << " evals/sec\n";
        
        CHECK_CLOSE_FRACTION(result.value, bc.exact, tol * 10);
    }
}

// Main function
int main()
{
    std::cout << "\n=== Infinite Domain Integration Tests ===\n\n";
    
    test_semi_infinite_integrals();
    test_infinite_integrals();
    test_oscillatory_integrals();
    test_transform_comparison();
    test_edge_cases();
    test_convergence_rates();
    test_performance_benchmark();
    
    std::cout << "\n=== Test Summary ===\n";
    std::cout << "Total tests: " << test_count << "\n";
    std::cout << "Passed: " << pass_count << "\n";
    std::cout << "Failed: " << fail_count << "\n";
    
    return fail_count > 0 ? 1 : 0;
}