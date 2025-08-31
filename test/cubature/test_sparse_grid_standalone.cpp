// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Standalone test for sparse grid integration
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>

using boost::math::cubature::hypercube;
using boost::math::cubature::integrate_sparse_grid;
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
void test_polynomial_exactness(TestRunner& runner) {
    runner.start_test("Polynomial Exactness");
    
    // Test with 2D monomial x^2 * y^2 (total degree 4)
    auto f = [](const double* x, std::size_t) -> double {
        return x[0] * x[0] * x[1] * x[1];
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // Level 3 should handle total degree 5 exactly
    auto res = integrate_sparse_grid<double>(f, box, 3);
    
    // Exact integral: (1/3) * (1/3) = 1/9
    double expected = 1.0 / 9.0;
    runner.check_close(res.value, expected, 1e-10, "Polynomial integral");
}

void test_smooth_function(TestRunner& runner) {
    runner.start_test("Smooth Function Convergence");
    
    // Test with smooth exponential function
    auto f = [](const double* x, std::size_t) -> double {
        return std::exp(x[0] + x[1]);
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    double expected = (std::exp(1.0) - 1.0) * (std::exp(1.0) - 1.0);
    
    // Test with level 4
    auto res = integrate_sparse_grid<double>(f, box, 4);
    runner.check_close(res.value, expected, 1e-4, "Smooth function convergence");
}

void test_node_count(TestRunner& runner) {
    runner.start_test("Node Count Growth");
    
    // Simple constant function to count evaluations
    std::size_t eval_count = 0;
    auto f = [&eval_count](const double*, std::size_t) -> double {
        eval_count++;
        return 1.0;
    };
    
    hypercube<double> box(3);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // Test level 2
    eval_count = 0;
    auto res = integrate_sparse_grid<double>(f, box, 2);
    
    // Should integrate to 1 exactly
    runner.check_close(res.value, 1.0, 1e-10, "Constant integral");
    
    // Verify evaluation count matches
    runner.check(res.evaluations == eval_count, "Evaluation count matches");
    runner.check(eval_count > 0, "Has evaluations");
}

void test_high_dimension(TestRunner& runner) {
    runner.start_test("High Dimensional Integration");
    
    // Product of cosines - smooth periodic function
    auto f = [](const double* x, std::size_t dim) -> double {
        double prod = 1.0;
        for (std::size_t i = 0; i < dim; ++i) {
            prod *= std::cos(pi<double>() * x[i]);
        }
        return prod;
    };
    
    // Test in 4D
    hypercube<double> box(4);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    auto res = integrate_sparse_grid<double>(f, box, 3);
    
    // Integral should be 0 (odd function integrated symmetrically)
    runner.check(std::abs(res.value) < 1e-6, "Oscillatory integral near zero");
}

void test_weight_diagnostics(TestRunner& runner) {
    runner.start_test("Weight Diagnostics");
    
    // Sharp Gaussian that may trigger weight cancellation
    auto f = [](const double* x, std::size_t) -> double {
        double r2 = 100.0 * ((x[0] - 0.5) * (x[0] - 0.5) + 
                             (x[1] - 0.5) * (x[1] - 0.5));
        return std::exp(-r2);
    };
    
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    // High level may cause issues with localized function
    auto res = integrate_sparse_grid<double>(f, box, 5);
    
    // Just check that it doesn't crash and gives some result
    runner.check(std::isfinite(res.value), "Result is finite");
}

void test_edge_cases(TestRunner& runner) {
    runner.start_test("Edge Cases");
    
    auto f = [](const double* x, std::size_t) -> double {
        return x[0] + x[1];
    };
    
    // Test with level 0 (minimal grid)
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    auto res = integrate_sparse_grid<double>(f, box, 0);
    // Level 0 uses minimal grid with 1 evaluation at center
    runner.check(res.evaluations == 1,
                "Level 0 should use minimal grid (1 evaluation)");
    
    // Test with 0 dimensions
    try {
        hypercube<double> box0(0);
        auto res0 = integrate_sparse_grid<double>(f, box0, 3);
        runner.check(res0.status != boost::math::cubature::status_code::success,
                    "Zero dimension should fail");
    } catch (const std::exception& e) {
        runner.check(true, "Zero dimension throws as expected");
    }
}

void test_cc_nodes(TestRunner& runner) {
    runner.start_test("Clenshaw-Curtis Nodes");
    
    // Test that CC nodes are symmetric and nested
    auto f = [](const double* x, std::size_t) -> double {
        return 1.0;
    };
    
    hypercube<double> box(1);
    box.lower = {-1.0};
    box.upper = {1.0};
    
    // Level 2 should include x = -1, 0, 1
    auto res = integrate_sparse_grid<double>(f, box, 2);
    
    // Check integral of constant
    runner.check_close(res.value, 2.0, 1e-10, "Constant integral over [-1,1]");
}

int main() {
    std::cout << "=== Sparse Grid Integration Standalone Tests ===\n\n";
    
    TestRunner runner;
    
    // Run all tests
    test_polynomial_exactness(runner);
    test_smooth_function(runner);
    test_node_count(runner);
    test_high_dimension(runner);
    test_weight_diagnostics(runner);
    test_edge_cases(runner);
    test_cc_nodes(runner);
    
    // Print summary
    runner.summary();
    
    return runner.exit_code();
}