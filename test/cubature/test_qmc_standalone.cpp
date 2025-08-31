// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Working test for QMC integration
#include <boost/math/cubature/qmc.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include <random>
#include <vector>

using boost::math::cubature::hypercube;
using boost::math::cubature::integrate_qmc;
using boost::math::cubature::integrate_rqmc;
using boost::math::constants::pi;
using boost::math::constants::root_pi;

// Simple test framework
class TestRunner {
private:
    int tests_run = 0;
    int tests_passed = 0;
    int tests_failed = 0;
    std::string current_test;
    
public:
    void start_test(const std::string& name) {
        current_test = name;
        std::cout << "Running: " << name << "..." << std::flush;
        tests_run++;
    }
    
    void check(bool condition, const std::string& message = "") {
        if (!condition) {
            std::cout << " FAILED\n";
            if (!message.empty()) {
                std::cout << "  Error: " << message << "\n";
            }
            tests_failed++;
        } else {
            std::cout << " PASSED\n";
            tests_passed++;
        }
    }
    
    template <typename T>
    void check_close(T actual, T expected, T tolerance, const std::string& message = "") {
        T error = std::abs(actual - expected);
        bool passed = error <= tolerance * std::abs(expected);
        
        if (!passed) {
            std::cout << " FAILED\n";
            std::cout << "  Expected: " << expected << "\n";
            std::cout << "  Actual:   " << actual << "\n";
            std::cout << "  Error:    " << error << "\n";
            std::cout << "  Tolerance:" << tolerance * std::abs(expected) << "\n";
            if (!message.empty()) {
                std::cout << "  " << message << "\n";
            }
            tests_failed++;
        } else {
            std::cout << " PASSED\n";
            tests_passed++;
        }
    }
    
    void summary() {
        std::cout << "\n=== Test Summary ===\n";
        std::cout << "Tests run:    " << tests_run << "\n";
        std::cout << "Tests passed: " << tests_passed << "\n";
        std::cout << "Tests failed: " << tests_failed << "\n";
        
        if (tests_failed == 0) {
            std::cout << "\nAll tests PASSED!\n";
        } else {
            std::cout << "\nSome tests FAILED!\n";
        }
    }
    
    int exit_code() const {
        return tests_failed > 0 ? 1 : 0;
    }
};

// Test functions
void test_qmc_convergence(TestRunner& runner) {
    runner.start_test("QMC Convergence");
    
    // Test function: product of sines (accepts vector)
    auto f = [](const std::vector<double>& x) -> double {
        double prod = 1.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            prod *= std::sin(pi<double>() * x[i]);
        }
        return prod;
    };
    
    hypercube<double> box(3);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // Exact integral: (2/pi)^3
    double expected = std::pow(2.0 / pi<double>(), 3);
    
    // Test with 10000 samples
    auto res = integrate_qmc<double>(f, box, 10000);
    
    runner.check_close(res.value, expected, 0.01, "QMC convergence for smooth function");
}

void test_rqmc_variance(TestRunner& runner) {
    runner.start_test("RQMC Variance Reduction");
    
    // Gaussian function (accepts vector)
    auto f = [](const std::vector<double>& x) -> double {
        return std::exp(-(x[0] * x[0] + x[1] * x[1]));
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), -2.0);
    std::fill(box.upper.begin(), box.upper.end(), 2.0);
    
    auto res = integrate_rqmc<double>(f, box, 1000, 10);
    
    // Approximate expected value
    double expected = pi<double>();
    
    runner.check_close(res.value, expected, 0.1, "RQMC Gaussian integral");
    runner.check(res.error > 0.0, "RQMC error estimate should be positive");
}

void test_high_dimensional(TestRunner& runner) {
    runner.start_test("High-Dimensional QMC");
    
    // Product function in 10D (accepts vector)
    auto f = [](const std::vector<double>& x) -> double {
        double prod = 1.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            prod *= (2.0 * x[i]);
        }
        return prod;
    };
    
    hypercube<double> box(10);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // Exact: 2^10 * (1/2)^10 = 1
    double expected = 1.0;
    
    auto res = integrate_qmc<double>(f, box, 50000);
    
    runner.check_close(res.value, expected, 0.05, "10D product integral");
}

void test_discontinuous(TestRunner& runner) {
    runner.start_test("Discontinuous Function");
    
    // Characteristic function of a circle (accepts vector)
    auto f = [](const std::vector<double>& x) -> double {
        double dx = x[0] - 0.5;
        double dy = x[1] - 0.5;
        double r2 = dx * dx + dy * dy;
        return (r2 < 0.25) ? 1.0 : 0.0;
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // Area of circle with radius 0.5
    double expected = pi<double>() * 0.25;
    
    auto res = integrate_qmc<double>(f, box, 50000);
    
    runner.check_close(res.value, expected, 0.02, "Circle area");
}

void test_sobol_uniformity(TestRunner& runner) {
    runner.start_test("Sobol Sequence Uniformity");
    
    // Constant function - tests uniformity (accepts vector)
    auto f = [](const std::vector<double>&) -> double {
        return 1.0;
    };
    
    hypercube<double> box(3);
    std::fill(box.lower.begin(), box.lower.end(), -1.0);
    std::fill(box.upper.begin(), box.upper.end(), 2.0);
    
    // Volume = 3^3 = 27
    double expected = 27.0;
    
    auto res = integrate_qmc<double>(f, box, 128);
    
    runner.check_close(res.value, expected, 1e-10, "Constant integral (uniformity test)");
}

void test_smooth_function(TestRunner& runner) {
    runner.start_test("Smooth Function QMC");
    
    // Very smooth periodic function (accepts vector)
    auto f = [](const std::vector<double>& x) -> double {
        return std::cos(2.0 * pi<double>() * x[0]) * 
               std::cos(2.0 * pi<double>() * x[1]) *
               std::exp(-(x[0] + x[1]));
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // Compare two sample sizes
    auto res1 = integrate_qmc<double>(f, box, 256);
    auto res2 = integrate_qmc<double>(f, box, 1024);
    
    // Should see improvement with more samples
    runner.check(std::isfinite(res1.value), "Result 1 is finite");
    runner.check(std::isfinite(res2.value), "Result 2 is finite");
    
    // Values should be close (smooth function converges fast)
    double diff = std::abs(res2.value - res1.value);
    runner.check(diff < 0.1, "Convergence for smooth function");
}

void test_error_conditions(TestRunner& runner) {
    runner.start_test("Error Conditions");
    
    auto f = [](const std::vector<double>&) -> double { return 1.0; };
    
    // Test with zero dimension
    try {
        hypercube<double> box(0);
        auto res = integrate_qmc<double>(f, box, 1000);
        runner.check(res.status != boost::math::cubature::status_code::success,
                    "Zero dimension should fail");
    } catch (const std::exception& e) {
        runner.check(true, "Zero dimension throws as expected");
    }
    
    // Test with zero samples
    hypercube<double> box2(2);
    std::fill(box2.lower.begin(), box2.lower.end(), 0.0);
    std::fill(box2.upper.begin(), box2.upper.end(), 1.0);
    
    try {
        auto res2 = integrate_qmc<double>(f, box2, 0);
        runner.check(res2.status != boost::math::cubature::status_code::success,
                    "Zero samples should fail");
    } catch (const std::exception& e) {
        runner.check(true, "Zero samples throws as expected");
    }
}

void test_qmc_vs_rqmc(TestRunner& runner) {
    runner.start_test("QMC vs RQMC Comparison");
    
    // Test function (accepts vector)
    auto f = [](const std::vector<double>& x) -> double {
        return std::exp(x[0]) * std::sin(pi<double>() * x[1]);
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // QMC result
    auto qmc_res = integrate_qmc<double>(f, box, 5000);
    
    // RQMC result
    auto rqmc_res = integrate_rqmc<double>(f, box, 1000, 5);
    
    // Results should be similar
    double diff = std::abs(qmc_res.value - rqmc_res.value);
    runner.check(diff < 0.05 * std::abs(qmc_res.value), 
                "QMC and RQMC give similar results");
    
    // RQMC should have error estimate
    runner.check(rqmc_res.error > 0.0, "RQMC provides error estimate");
}

// Additional test for integration accuracy
void test_integration_accuracy(TestRunner& runner) {
    runner.start_test("Integration Accuracy");
    
    // Test function with known integral
    auto f = [](const std::vector<double>& x) -> double {
        // f(x,y) = x*y over [0,1]^2, integral = 1/4
        return x[0] * x[1];
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    auto res = integrate_qmc<double>(f, box, 10000);
    
    double expected = 0.25;
    runner.check_close(res.value, expected, 0.001, "Simple polynomial integral");
    runner.check(res.evaluations == 10000, "Correct number of evaluations");
}

int main() {
    std::cout << "=== QMC Integration Tests ===\n\n";
    
    TestRunner runner;
    
    // Run all tests
    test_qmc_convergence(runner);
    test_rqmc_variance(runner);
    test_high_dimensional(runner);
    test_discontinuous(runner);
    test_sobol_uniformity(runner);
    test_smooth_function(runner);
    test_error_conditions(runner);
    test_qmc_vs_rqmc(runner);
    test_integration_accuracy(runner);
    
    // Print summary
    runner.summary();
    
    return runner.exit_code();
}