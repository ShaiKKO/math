// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Standalone test for transforms
#include <boost/math/cubature/transforms.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include <limits>

using boost::math::cubature::rational_transform;
using boost::math::cubature::tangent_transform;
using boost::math::cubature::exponential_transform;
using boost::math::cubature::algebraic_transform;
using boost::math::cubature::duffy_transform;
using boost::math::cubature::tent;
using boost::math::constants::pi;

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
void test_rational_transform(TestRunner& runner) {
    runner.start_test("Rational Transform");
    
    rational_transform<double> transform;
    
    // Test forward transform
    auto result = transform.forward(0.0);
    runner.check(std::abs(result.first) < 1e-10, "u=0 maps to x=0");
    
    result = transform.forward(0.5);
    runner.check_close(result.first, 1.0, 1e-10, "u=0.5 maps to x=1");
    runner.check_close(result.second, 4.0, 1e-10, "Jacobian at u=0.5");
    
    // Test inverse transform
    double u = transform.inverse(0.0);
    runner.check(std::abs(u) < 1e-10, "x=0 maps to u=0");
    
    u = transform.inverse(1.0);
    runner.check_close(u, 0.5, 1e-10, "x=1 maps to u=0.5");
    
    // Test round-trip consistency
    for (double u_test : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        auto [x, jac] = transform.forward(u_test);
        double u_back = transform.inverse(x);
        runner.check_close(u_back, u_test, 1e-10, "Round-trip consistency");
    }
}

void test_tangent_transform(TestRunner& runner) {
    runner.start_test("Tangent Transform");
    
    tangent_transform<double> transform;
    
    // Test forward transform
    auto result = transform.forward(0.5);
    runner.check(std::abs(result.first) < 1e-10, "u=0.5 maps to x=0");
    runner.check_close(result.second, pi<double>(), 1e-10, "Jacobian at u=0.5");
    
    result = transform.forward(0.25);
    runner.check_close(result.first, -1.0, 1e-8, "u=0.25 maps to x=-1");
    
    // Test inverse transform
    double u = transform.inverse(0.0);
    runner.check_close(u, 0.5, 1e-10, "x=0 maps to u=0.5");
    
    // Test round-trip
    for (double u_test : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        auto [x, jac] = transform.forward(u_test);
        double u_back = transform.inverse(x);
        runner.check_close(u_back, u_test, 1e-8, "Round-trip consistency");
    }
}

void test_exponential_transform(TestRunner& runner) {
    runner.start_test("Exponential Transform");
    
    exponential_transform<double> transform;
    
    // Test forward transform
    auto result = transform.forward(0.0);
    runner.check(std::abs(result.first) < 1e-10, "u=0 maps to x=0");
    runner.check_close(result.second, 1.0, 1e-10, "Jacobian at u=0");
    
    // Test inverse transform
    double u = transform.inverse(0.0);
    runner.check(std::abs(u) < 1e-10, "x=0 maps to u=0");
    
    u = transform.inverse(1.0);
    double expected_u = 1.0 - std::exp(-1.0);
    runner.check_close(u, expected_u, 1e-10, "x=1 maps correctly");
    
    // Test round-trip
    for (double u_test : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        auto [x, jac] = transform.forward(u_test);
        double u_back = transform.inverse(x);
        runner.check_close(u_back, u_test, 1e-10, "Round-trip consistency");
    }
}

void test_algebraic_transform(TestRunner& runner) {
    runner.start_test("Algebraic Transform");
    
    // Test with α = β = 1 (reduces to rational transform)
    algebraic_transform<double> transform1(1.0, 1.0);
    
    auto result = transform1.forward(0.5);
    runner.check_close(result.first, 1.0, 1e-10, "α=β=1, u=0.5 maps to x=1");
    
    double u = transform1.inverse(1.0);
    runner.check_close(u, 0.5, 1e-10, "α=β=1, x=1 maps to u=0.5");
    
    // Test with α = 2, β = 1
    algebraic_transform<double> transform2(2.0, 1.0);
    
    result = transform2.forward(0.5);
    runner.check_close(result.first, 0.5, 1e-10, "α=2,β=1, u=0.5 maps to x=0.5");
    
    // Test boundary behavior
    result = transform2.forward(0.0);
    runner.check(std::abs(result.first) < 1e-10, "u=0 maps to x=0");
    runner.check(std::abs(result.second) < 1e-10, "Jacobian at u=0");
}

void test_duffy_transform(TestRunner& runner) {
    runner.start_test("Duffy Transform");
    
    // Test 1D (identity)
    {
        double u[1] = {0.5};
        double x[1];
        double jac = duffy_transform<double>::apply(u, x, 1);
        runner.check_close(x[0], 0.5, 1e-10, "1D Duffy is identity");
        runner.check_close(jac, 1.0, 1e-10, "1D Jacobian is 1");
    }
    
    // Test 2D
    {
        double u[2] = {0.5, 0.5};
        double x[2];
        double jac = duffy_transform<double>::apply(u, x, 2);
        
        runner.check_close(x[0], 0.25, 1e-10, "2D x[0] = u[0]*(1-u[1])");
        runner.check_close(x[1], 0.25, 1e-10, "2D x[1] = u[0]*u[1]");
        runner.check_close(jac, 0.5, 1e-10, "2D Jacobian = u[0]");
        
        // Test that points map to simplex
        runner.check(x[0] >= 0.0 && x[1] >= 0.0 && x[0] + x[1] <= 1.0 + 1e-10,
                    "Points map to simplex");
    }
    
    // Test inverse
    {
        double x[2] = {0.3, 0.2};
        double u[2];
        duffy_transform<double>::inverse(x, u, 2);
        
        // Correct inverse: u[0] = x[0] + x[1], u[1] = x[1]/(x[0] + x[1])
        runner.check_close(u[0], 0.5, 1e-10, "Inverse u[0] = sum(x)");
        runner.check_close(u[1], 0.4, 1e-10, "Inverse u[1] = x[1]/sum(x)");
        
        // Verify round-trip
        double x_back[2];
        duffy_transform<double>::apply(u, x_back, 2);
        runner.check_close(x_back[0], x[0], 1e-10, "Round-trip x[0]");
        runner.check_close(x_back[1], x[1], 1e-10, "Round-trip x[1]");
    }
}

void test_tent_transform(TestRunner& runner) {
    runner.start_test("Tent Transform");
    
    // Test tent transform: tent(u) = 1 - 2*|u - 0.5|
    runner.check_close(tent(0.0), 0.0, 1e-10, "tent(0) = 0");
    runner.check_close(tent(0.25), 0.5, 1e-10, "tent(0.25) = 0.5");
    runner.check_close(tent(0.5), 1.0, 1e-10, "tent(0.5) = 1");
    runner.check_close(tent(0.75), 0.5, 1e-10, "tent(0.75) = 0.5");
    runner.check_close(tent(1.0), 0.0, 1e-10, "tent(1) = 0");
    
    // Test symmetry
    for (double u : {0.1, 0.2, 0.3, 0.4}) {
        double t1 = tent(u);
        double t2 = tent(1.0 - u);
        runner.check_close(t1, t2, 1e-10, "Tent symmetry");
    }
    
    // Test that it's always in [0,1]
    for (double u : {0.0, 0.1, 0.5, 0.9, 1.0}) {
        double t = tent(u);
        runner.check(t >= 0.0 && t <= 1.0, "Tent output in [0,1]");
    }
}

int main() {
    std::cout << "=== Transforms Standalone Tests ===\n\n";
    
    TestRunner runner;
    
    // Run all tests
    test_rational_transform(runner);
    test_tangent_transform(runner);
    test_exponential_transform(runner);
    test_algebraic_transform(runner);
    test_duffy_transform(runner);
    test_tent_transform(runner);
    
    // Print summary
    runner.summary();
    
    return runner.exit_code();
}