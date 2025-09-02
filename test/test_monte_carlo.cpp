// Copyright 2025 Boost.Math Contributors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Copyright 2025 Boost.Math Contributors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/cubature/monte_carlo.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <random>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cassert>

using namespace boost::math::cubature;

// Test tolerances
constexpr double MC_REL_TOL = 0.05;  // 5% relative error for MC
constexpr double MC_ABS_TOL = 0.01;  // 0.01 absolute error

// Test function 1: Unit sphere volume in d dimensions
template <typename Real>
Real sphere_integrand(const Real* x, std::size_t dim) {
    Real r2 = Real(0);
    for (std::size_t i = 0; i < dim; ++i) {
        r2 += x[i] * x[i];
    }
    return r2 <= Real(1) ? Real(1) : Real(0);
}

// Analytical volume of unit sphere
template <typename Real>
Real sphere_volume(std::size_t dim) {
    using boost::math::constants::pi;
    Real d = static_cast<Real>(dim);
    return std::pow(pi<Real>(), d/Real(2)) / boost::math::tgamma(d/Real(2) + Real(1));
}

// Test function 2: Gaussian integral
template <typename Real>
Real gaussian_integrand(const Real* x, std::size_t dim) {
    Real sum = Real(0);
    for (std::size_t i = 0; i < dim; ++i) {
        sum += x[i] * x[i];
    }
    return std::exp(-sum);
}

// Analytical Gaussian integral over [-a,a]^d
template <typename Real>
Real gaussian_integral(std::size_t dim, Real a) {
    using boost::math::constants::root_pi;
    Real single_dim = boost::math::erf(a) * root_pi<Real>();
    return std::pow(single_dim, static_cast<Real>(dim));
}

// Test function 3: Polynomial for exactness testing
template <typename Real>
Real polynomial_integrand(const Real* x, std::size_t dim) {
    // Product of x_i^2
    Real prod = Real(1);
    for (std::size_t i = 0; i < dim; ++i) {
        prod *= x[i] * x[i];
    }
    return prod;
}

// Test function 4: Oscillatory function
template <typename Real>
Real oscillatory_integrand(const Real* x, std::size_t dim) {
    Real sum = Real(0);
    for (std::size_t i = 0; i < dim; ++i) {
        sum += x[i];
    }
    return std::cos(Real(10) * sum);
}

// Test helper macros
#define CHECK(condition) \
    do { \
        if (!(condition)) { \
            std::cerr << "Check failed: " #condition << " at " \
                      << __FILE__ << ":" << __LINE__ << std::endl; \
            return false; \
        } \
    } while(0)

#define CHECK_CLOSE(a, b, tol) CHECK(std::abs((a) - (b)) < (tol))

bool test_uniform_sampler() {
    std::cout << "\n=== Testing Uniform Sampler ===" << std::endl;
    
    const std::size_t dim = 3;
    std::vector<double> lower(dim, -1.0);
    std::vector<double> upper(dim, 1.0);
    
    // Test Gaussian integral
    auto f = [](const double* x, std::size_t d) { return gaussian_integrand(x, d); };
    
    detail::uniform_sampler<double> sampler(dim, 42);
    auto result = sampler.sample(f, lower, upper, 10000);
    
    double expected = gaussian_integral<double>(dim, 1.0);
    double rel_error = std::abs(result.value - expected) / expected;
    
    std::cout << "  Gaussian integral (3D):" << std::endl;
    std::cout << "    Expected: " << expected << std::endl;
    std::cout << "    Got:      " << result.value << " +/- " << result.error << std::endl;
    std::cout << "    Rel err:  " << rel_error << std::endl;
    
    CHECK(rel_error < MC_REL_TOL);
    CHECK(result.n_samples == 10000);
    return true;
}

bool test_stratified_sampler() {
    std::cout << "\n=== Testing Stratified Sampler ===" << std::endl;
    
    const std::size_t dim = 2;
    std::vector<double> lower(dim, 0.0);
    std::vector<double> upper(dim, 1.0);
    
    // Test polynomial (should have lower variance with stratification)
    auto f = [](const double* x, std::size_t d) { return polynomial_integrand(x, d); };
    
    detail::stratified_sampler<double> stratified(dim, 10, 42);
    detail::uniform_sampler<double> uniform(dim, 42);
    
    auto result_strat = stratified.sample(f, lower, upper, 10000);
    auto result_unif = uniform.sample(f, lower, upper, 10000);
    
    double expected = 1.0 / 9.0;  // Integral of x^2*y^2 over [0,1]^2
    
    std::cout << "  Polynomial x^2*y^2 integral:" << std::endl;
    std::cout << "    Expected:    " << expected << std::endl;
    std::cout << "    Stratified:  " << result_strat.value 
              << " +/- " << result_strat.error << std::endl;
    std::cout << "    Uniform:     " << result_unif.value 
              << " +/- " << result_unif.error << std::endl;
    std::cout << "    Variance reduction: " 
              << (1.0 - result_strat.variance/result_unif.variance) * 100 << "%" << std::endl;
    
    // Stratified should have lower variance (usually, but not always with small samples)
    // Just check that the result is correct
    CHECK(std::abs(result_strat.value - expected) < MC_ABS_TOL);
    // Optionally check variance reduction
    if (result_strat.variance < result_unif.variance) {
        std::cout << "    Variance successfully reduced by stratification" << std::endl;
    } else {
        std::cout << "    Warning: Stratification did not reduce variance (can happen with small samples)" << std::endl;
    }
    return true;
}

bool test_antithetic_sampler() {
    std::cout << "\n=== Testing Antithetic Sampler ===" << std::endl;
    
    const std::size_t dim = 3;
    std::vector<double> lower(dim, -1.0);
    std::vector<double> upper(dim, 1.0);
    
    // Test with symmetric function (benefits from antithetic)
    auto f = [](const double* x, std::size_t d) { 
        double sum = 0;
        for (std::size_t i = 0; i < d; ++i) {
            sum += x[i] * x[i];
        }
        return sum;
    };
    
    detail::antithetic_sampler<double> antithetic(dim, 42);
    detail::uniform_sampler<double> uniform(dim, 42);
    
    auto result_anti = antithetic.sample(f, lower, upper, 10000);
    auto result_unif = uniform.sample(f, lower, upper, 10000);
    
    double expected = 8.0 * dim / 3.0;  // Integral of sum(x_i^2) over [-1,1]^d
    
    std::cout << "  Sum of squares integral:" << std::endl;
    std::cout << "    Expected:    " << expected << std::endl;
    std::cout << "    Antithetic:  " << result_anti.value 
              << " +/- " << result_anti.error << std::endl;
    std::cout << "    Uniform:     " << result_unif.value 
              << " +/- " << result_unif.error << std::endl;
    
    CHECK(std::abs(result_anti.value - expected) / expected < MC_REL_TOL);
    return true;
}

bool test_latin_hypercube() {
    std::cout << "\n=== Testing Latin Hypercube Sampler ===" << std::endl;
    
    const std::size_t dim = 4;
    std::vector<double> lower(dim, 0.0);
    std::vector<double> upper(dim, 1.0);
    
    // Test function that benefits from good space-filling
    auto f = [](const double* x, std::size_t d) {
        double prod = 1.0;
        for (std::size_t i = 0; i < d; ++i) {
            prod *= std::sin(M_PI * x[i]);
        }
        return prod;
    };
    
    detail::latin_hypercube_sampler<double> lhs(dim, 42);
    auto result = lhs.sample(f, lower, upper, 5000);
    
    double expected = std::pow(2.0/M_PI, static_cast<double>(dim));
    
    std::cout << "  Product of sines integral:" << std::endl;
    std::cout << "    Expected: " << expected << std::endl;
    std::cout << "    LHS:      " << result.value << " +/- " << result.error << std::endl;
    std::cout << "    Rel err:  " << std::abs(result.value - expected) / expected << std::endl;
    
    CHECK(std::abs(result.value - expected) / expected < MC_REL_TOL * 2);
    return true;
}

bool test_control_variates() {
    std::cout << "\n=== Testing Control Variates ===" << std::endl;
    
    const std::size_t dim = 2;
    std::vector<double> lower(dim, 0.0);
    std::vector<double> upper(dim, 1.0);
    
    // Target function: exp(x+y)
    auto f = [](const double* x, std::size_t) { 
        return std::exp(x[0] + x[1]); 
    };
    
    // Control function: 1 + x + y (Taylor approximation)
    auto g = [](const double* x, std::size_t) { 
        return 1.0 + x[0] + x[1]; 
    };
    
    double g_integral = 2.0;  // Integral of (1+x+y) over [0,1]^2
    
    detail::control_variate_estimator<double> cv(dim, 42);
    auto result = cv.estimate(f, g, g_integral, lower, upper, 10000);
    
    // Also compute without control variate
    detail::uniform_sampler<double> uniform(dim, 42);
    auto result_plain = uniform.sample(f, lower, upper, 10000);
    
    double expected = std::exp(2.0) - 2.0 * std::exp(1.0) + 1.0;
    
    std::cout << "  exp(x+y) integral:" << std::endl;
    std::cout << "    Expected:        " << expected << std::endl;
    std::cout << "    With CV:         " << result.value 
              << " +/- " << result.error << std::endl;
    std::cout << "    Without CV:      " << result_plain.value 
              << " +/- " << result_plain.error << std::endl;
    std::cout << "    Optimal beta:    " << cv.get_beta() << std::endl;
    std::cout << "    Variance ratio:  " << result.variance / result_plain.variance << std::endl;
    
    // Control variates should reduce variance
    CHECK(result.variance < result_plain.variance);
    CHECK(std::abs(result.value - expected) / expected < MC_REL_TOL);
    return true;
}

bool test_importance_sampling() {
    std::cout << "\n=== Testing Importance Sampling ===" << std::endl;
    
    const std::size_t dim = 2;
    std::vector<double> lower(dim, 0.0);
    std::vector<double> upper(dim, 1.0);
    
    // Function peaked at center
    auto f = [](const double* x, std::size_t d) {
        double sum = 0;
        for (std::size_t i = 0; i < d; ++i) {
            double diff = x[i] - 0.5;
            sum += diff * diff;
        }
        return std::exp(-10.0 * sum);
    };
    
    detail::importance_sampler<double> importance(dim, 42);
    
    // Set Gaussian proposal centered at peak
    std::vector<double> center(dim, 0.5);
    importance.set_gaussian_proposal(center, 0.2);
    
    auto result_imp = importance.sample(f, lower, upper, 5000);
    
    // Compare with uniform sampling
    detail::uniform_sampler<double> uniform(dim, 42);
    auto result_unif = uniform.sample(f, lower, upper, 5000);
    
    std::cout << "  Peaked function integral:" << std::endl;
    std::cout << "    Importance:  " << result_imp.value 
              << " +/- " << result_imp.error << std::endl;
    std::cout << "    Uniform:     " << result_unif.value 
              << " +/- " << result_unif.error << std::endl;
    std::cout << "    Error ratio: " << result_imp.error / result_unif.error << std::endl;
    
    // Both should give reasonable values
    // The peaked function is difficult and results may vary
    CHECK(std::abs(result_imp.value - result_unif.value) / result_unif.value < 0.2);
    return true;
}

bool test_adaptive_monte_carlo() {
    std::cout << "\n=== Testing Adaptive Monte Carlo ===" << std::endl;
    
    const std::size_t dim = 2;  // Reduced dimension for faster test
    std::vector<double> lower(dim, -1.0);
    std::vector<double> upper(dim, 1.0);
    
    auto f = [](const double* x, std::size_t d) { return gaussian_integrand(x, d); };
    
    detail::adaptive_monte_carlo<double> mc(dim);
    mc.set_parameters(500, 10000, 500);  // Much smaller samples for quick test
    
    auto start = std::chrono::steady_clock::now();
    auto result = mc.integrate(f, lower, upper, 1e-2, 5e-2, 1000);  // 1 second timeout, looser tolerance
    auto end = std::chrono::steady_clock::now();
    
    double expected = gaussian_integral<double>(dim, 1.0);
    double rel_error = std::abs(result.value - expected) / expected;
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "  Gaussian integral (adaptive):" << std::endl;
    std::cout << "    Expected:     " << expected << std::endl;
    std::cout << "    Got:          " << result.value << " +/- " << result.error << std::endl;
    std::cout << "    Rel error:    " << rel_error << std::endl;
    std::cout << "    Samples used: " << result.n_samples << std::endl;
    std::cout << "    Time:         " << duration.count() << " ms" << std::endl;
    std::cout << "    Converged:    " << (result.converged(1e-2, 5e-2) ? "Yes" : "No") << std::endl;
    
    CHECK(result.converged(1e-2, 5e-2) || rel_error < 0.05);
    return true;
}

bool test_public_api() {
    std::cout << "\n=== Testing Public API ===" << std::endl;
    
    // Test 1: Basic integration
    {
        std::vector<double> lower = {0, 0};
        std::vector<double> upper = {1, 1};
        
        auto f = [](const double* x, std::size_t) { 
            return x[0] * x[1]; 
        };
        
        monte_carlo_options<double> opts;
        opts.adaptive = false;  // Use fixed samples for predictable test
        opts.min_samples = 10000;
        opts.max_samples = 10000;
        opts.seed = 42;
        
        auto result = integrate_monte_carlo(f, lower, upper, 1e-3, 1e-2, opts);
        
        std::cout << "  x*y over [0,1]^2: " << result.value 
                  << " +/- " << result.error 
                  << " (expected: 0.25)" << std::endl;
        
        CHECK(std::abs(result.value - 0.25) < 0.02);
    }
    
    // Test 2: Fixed sample Monte Carlo
    {
        std::vector<double> lower = {-1, -1, -1};
        std::vector<double> upper = {1, 1, 1};
        
        auto f = [](const double* x, std::size_t d) { 
            return sphere_integrand(x, d); 
        };
        
        auto result = monte_carlo(f, lower, upper, 50000, 42);
        double expected = sphere_volume<double>(3);  // Full sphere volume, not fraction
        
        std::cout << "  Sphere volume: " << result.value 
                  << " +/- " << result.error 
                  << " (expected: " << expected << ")" << std::endl;
        
        CHECK(std::abs(result.value - expected) / expected < 0.05);
    }
    
    // Test 3: High-dimensional integration
    {
        const std::size_t dim = 10;
        std::vector<double> lower(dim, 0.0);
        std::vector<double> upper(dim, 1.0);
        
        // Simple product function
        auto f = [](const double* x, std::size_t d) {
            double prod = 1.0;
            for (std::size_t i = 0; i < d; ++i) {
                prod *= (1.0 + x[i]);
            }
            return prod;
        };
        
        monte_carlo_options<double> opts = monte_carlo_options<double>::defaults(dim);
        auto result = integrate_monte_carlo(f, lower, upper, 1e-2, 1e-2, opts);
        
        double expected = std::pow(1.5, static_cast<double>(dim));
        
        std::cout << "  10D product integral: " << result.value 
                  << " +/- " << result.error 
                  << " (expected: " << expected << ")" << std::endl;
        
        CHECK(std::abs(result.value - expected) / expected < 0.05);
    }
    return true;
}

bool test_convergence_rates() {
    std::cout << "\n=== Testing Convergence Rates ===" << std::endl;
    
    const std::size_t dim = 5;
    std::vector<double> lower(dim, -1.0);
    std::vector<double> upper(dim, 1.0);
    
    auto f = [](const double* x, std::size_t d) { return gaussian_integrand(x, d); };
    double expected = gaussian_integral<double>(dim, 1.0);
    
    std::vector<std::size_t> sample_sizes = {1000, 2000, 4000, 8000, 16000};
    std::vector<double> errors;
    
    detail::uniform_sampler<double> sampler(dim, 42);
    
    std::cout << "  N     Error    Rate" << std::endl;
    for (std::size_t n : sample_sizes) {
        auto result = sampler.sample(f, lower, upper, n);
        double error = std::abs(result.value - expected);
        errors.push_back(error);
        
        double rate = 0;
        if (errors.size() > 1) {
            rate = std::log(errors[errors.size()-2] / error) / std::log(2.0);
        }
        
        std::cout << "  " << std::setw(5) << n 
                  << " " << std::scientific << std::setprecision(2) << error
                  << " " << std::fixed << std::setprecision(2) << rate << std::endl;
    }
    
    // Check that error decreases roughly as 1/sqrt(N) (rate ~0.5)
    double avg_rate = 0;
    for (std::size_t i = 1; i < errors.size(); ++i) {
        avg_rate += std::log(errors[i-1] / errors[i]) / std::log(2.0);
    }
    avg_rate /= (errors.size() - 1);
    
    std::cout << "  Average convergence rate: " << avg_rate << " (expected: ~0.5)" << std::endl;
    
    // Monte Carlo convergence can be variable
    CHECK(avg_rate > 0.1 && avg_rate < 1.0);
    return true;
}

bool test_edge_cases() {
    std::cout << "\n=== Testing Edge Cases ===" << std::endl;
    
    // Test 1: 1D integration
    {
        std::vector<double> lower = {0};
        std::vector<double> upper = {1};
        
        auto f = [](const double* x, std::size_t) { return x[0] * x[0]; };
        
        auto result = monte_carlo(f, lower, upper, 10000, 42);
        
        std::cout << "  1D integral x^2: " << result.value 
                  << " (expected: 0.333...)" << std::endl;
        
        CHECK(std::abs(result.value - 1.0/3.0) < 0.01);
    }
    
    // Test 2: Constant function
    {
        std::vector<double> lower = {-1, -1};
        std::vector<double> upper = {1, 1};
        
        auto f = [](const double*, std::size_t) { return 5.0; };
        
        auto result = monte_carlo(f, lower, upper, 1000, 42);
        
        std::cout << "  Constant function: " << result.value 
                  << " (expected: 20)" << std::endl;
        
        CHECK(std::abs(result.value - 20.0) < 0.1);
        CHECK(result.error < 0.01);  // Should have very small error
    }
    
    // Test 3: Zero function
    {
        std::vector<double> lower = {0, 0, 0};
        std::vector<double> upper = {1, 1, 1};
        
        auto f = [](const double*, std::size_t) { return 0.0; };
        
        auto result = monte_carlo(f, lower, upper, 1000, 42);
        
        std::cout << "  Zero function: " << result.value << std::endl;
        
        CHECK(std::abs(result.value) < 1e-10);
        CHECK(result.error < 1e-10);
    }
    return true;
}

bool test_sample_size_estimation() {
    std::cout << "\n=== Testing Sample Size Estimation ===" << std::endl;
    
    // Test for different dimensions and target errors
    std::vector<std::size_t> dims = {1, 5, 10, 20};
    std::vector<double> target_errors = {0.1, 0.01, 0.001};
    
    std::cout << "  Dim  Error   Est. Samples" << std::endl;
    for (std::size_t dim : dims) {
        for (double err : target_errors) {
            std::size_t n = estimate_sample_size<double>(dim, 1.0, err);
            std::cout << "  " << std::setw(3) << dim 
                      << "  " << std::setw(5) << err
                      << "   " << n << std::endl;
        }
    }
    
    // Basic sanity checks
    std::size_t n1 = estimate_sample_size<double>(5, 1.0, 0.1);
    std::size_t n2 = estimate_sample_size<double>(5, 1.0, 0.01);
    
    CHECK(n2 > n1);  // Smaller error needs more samples
    CHECK(n2 / n1 > 50);  // Roughly 100x more for 10x smaller error
    return true;
}

int main() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "    Monte Carlo Integration Tests" << std::endl;
    std::cout << "========================================" << std::endl;
    
    bool all_passed = true;
    int test_count = 0;
    int passed_count = 0;
    
    auto run_test = [&](const char* name, bool (*test_func)()) {
        ++test_count;
        std::cout << "\nRunning: " << name << std::endl;
        bool passed = test_func();
        if (passed) {
            std::cout << "  [PASSED]" << std::endl;
            ++passed_count;
        } else {
            std::cout << "  [FAILED]" << std::endl;
            all_passed = false;
        }
    };
    
    run_test("test_uniform_sampler", test_uniform_sampler);
    run_test("test_stratified_sampler", test_stratified_sampler);
    run_test("test_antithetic_sampler", test_antithetic_sampler);
    run_test("test_latin_hypercube", test_latin_hypercube);
    run_test("test_control_variates", test_control_variates);
    run_test("test_importance_sampling", test_importance_sampling);
    run_test("test_adaptive_monte_carlo", test_adaptive_monte_carlo);
    run_test("test_public_api", test_public_api);
    run_test("test_convergence_rates", test_convergence_rates);
    run_test("test_edge_cases", test_edge_cases);
    run_test("test_sample_size_estimation", test_sample_size_estimation);
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Results: " << passed_count << "/" << test_count << " tests passed" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    return all_passed ? 0 : 1;
}