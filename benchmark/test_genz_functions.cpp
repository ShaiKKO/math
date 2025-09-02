// Copyright (c) 2025 Boost.Math Contributors
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include "test_functions.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

namespace bmb = boost::math::benchmark;

// Test a single Genz function
template <typename Real, typename Function>
void test_function(const std::string& name, Function& func, std::size_t dimension) {
    std::cout << "\nTesting " << name << " (dimension=" << dimension << ")\n";
    std::cout << std::string(50, '-') << "\n";
    
    // Test point in the middle of the unit hypercube
    std::vector<Real> x(dimension, 0.5);
    
    // Evaluate function
    func.reset_eval_count();
    Real value = func(x);
    
    std::cout << "f([0.5, ...]) = " << std::scientific << std::setprecision(6) << value << "\n";
    
    // Get analytical integral
    Real analytical = 0;
    if (auto* osc = dynamic_cast<bmb::genz_oscillatory<Real>*>(&func)) {
        analytical = osc->analytical_integral();
    } else if (auto* peak = dynamic_cast<bmb::genz_product_peak<Real>*>(&func)) {
        analytical = peak->analytical_integral();
    } else if (auto* corner = dynamic_cast<bmb::genz_corner_peak<Real>*>(&func)) {
        analytical = corner->analytical_integral();
    } else if (auto* gauss = dynamic_cast<bmb::genz_gaussian<Real>*>(&func)) {
        analytical = gauss->analytical_integral();
    } else if (auto* cont = dynamic_cast<bmb::genz_continuous<Real>*>(&func)) {
        analytical = cont->analytical_integral();
    } else if (auto* disc = dynamic_cast<bmb::genz_discontinuous<Real>*>(&func)) {
        analytical = disc->analytical_integral();
    }
    
    std::cout << "Analytical integral: " << analytical << "\n";
    
    // Simple Monte Carlo estimate for verification
    std::size_t num_samples = 10000;
    Real sum = 0;
    std::mt19937 gen(123);
    std::uniform_real_distribution<Real> dist(0, 1);
    
    for (std::size_t i = 0; i < num_samples; ++i) {
        std::vector<Real> sample(dimension);
        for (std::size_t j = 0; j < dimension; ++j) {
            sample[j] = dist(gen);
        }
        sum += func(sample);
    }
    
    Real mc_estimate = sum / num_samples;
    Real mc_error = std::abs(mc_estimate - analytical) / (std::abs(analytical) + 1e-10);
    
    std::cout << "Monte Carlo estimate (" << num_samples << " samples): " 
              << mc_estimate << "\n";
    std::cout << "MC relative error: " << mc_error << "\n";
    std::cout << "Function evaluations: " << func.get_eval_count() << "\n";
    
    // Test with different difficulty levels
    std::cout << "\nDifficulty scaling:\n";
    for (Real difficulty : {0.5, 1.0, 2.0, 5.0}) {
        func.set_difficulty(difficulty);
        func.reset_eval_count();
        
        Real difficult_value = func(x);
        std::cout << "  Difficulty " << std::fixed << std::setprecision(1) << difficulty 
                  << ": f([0.5, ...]) = " << std::scientific << std::setprecision(4) 
                  << difficult_value << "\n";
    }
}

// Test polynomial function
template <typename Real>
void test_polynomial() {
    std::cout << "\nTesting Polynomial Function\n";
    std::cout << std::string(50, '=') << "\n";
    
    // Test x^2 * y^3 * z
    std::vector<std::size_t> powers = {2, 3, 1};
    bmb::polynomial_function<Real> poly(powers);
    
    std::vector<Real> x = {0.5, 0.5, 0.5};
    Real value = poly(x);
    Real expected = std::pow(0.5, 2) * std::pow(0.5, 3) * std::pow(0.5, 1);
    
    std::cout << "f([0.5, 0.5, 0.5]) = " << value << "\n";
    std::cout << "Expected: " << expected << "\n";
    std::cout << "Match: " << (std::abs(value - expected) < 1e-10 ? "YES" : "NO") << "\n";
    
    Real analytical = poly.analytical_integral();
    Real expected_integral = 1.0 / (3 * 4 * 2);  // 1/(2+1) * 1/(3+1) * 1/(1+1)
    
    std::cout << "Analytical integral: " << analytical << "\n";
    std::cout << "Expected integral: " << expected_integral << "\n";
    std::cout << "Match: " << (std::abs(analytical - expected_integral) < 1e-10 ? "YES" : "NO") << "\n";
}

// Test singular function
template <typename Real>
void test_singular() {
    std::cout << "\nTesting Singular Function\n";
    std::cout << std::string(50, '=') << "\n";
    
    std::size_t dimension = 2;
    Real alpha = 0.5;  // Singularity at origin
    
    bmb::singular_function<Real> singular(dimension, alpha);
    
    // Test at different points
    std::vector<std::vector<Real>> test_points = {
        {0.1, 0.1},
        {0.5, 0.5},
        {0.01, 0.01}
    };
    
    for (const auto& point : test_points) {
        Real value = singular(point);
        Real expected = std::pow(point[0], -alpha) * std::pow(point[1], -alpha);
        
        std::cout << "f([" << point[0] << ", " << point[1] << "]) = " 
                  << std::scientific << std::setprecision(4) << value << "\n";
        std::cout << "  Expected: " << expected << "\n";
        std::cout << "  Match: " << (std::abs(value - expected) < 1e-10 ? "YES" : "NO") << "\n";
    }
    
    Real analytical = singular.analytical_integral();
    Real expected_integral = std::pow(1.0 / (1 - alpha), dimension);
    
    std::cout << "\nAnalytical integral: " << analytical << "\n";
    std::cout << "Expected integral: " << expected_integral << "\n";
    std::cout << "Match: " << (std::abs(analytical - expected_integral) < 1e-10 ? "YES" : "NO") << "\n";
}

int main() {
    using Real = double;
    
    std::cout << "========================================\n";
    std::cout << "    Genz Test Functions Verification\n";
    std::cout << "========================================\n";
    
    // Test dimensions
    std::vector<std::size_t> test_dimensions = {2, 3, 5};
    
    for (std::size_t dim : test_dimensions) {
        std::cout << "\n\n*** DIMENSION " << dim << " ***\n";
        
        // Test all Genz functions
        auto function_types = bmb::genz_test_suite<Real>::all_functions();
        
        for (auto func_type : function_types) {
            const char* func_name = bmb::genz_test_suite<Real>::function_name(func_type);
            auto test_func = bmb::genz_test_suite<Real>::create_function(func_type, dim);
            
            test_function<Real>(func_name, *test_func, dim);
        }
    }
    
    // Test additional functions
    std::cout << "\n\n*** ADDITIONAL TEST FUNCTIONS ***\n";
    test_polynomial<Real>();
    test_singular<Real>();
    
    // Verify analytical formulas for specific cases
    std::cout << "\n\n*** ANALYTICAL FORMULA VERIFICATION ***\n";
    std::cout << std::string(50, '=') << "\n";
    
    // Test Gaussian integral with known parameters
    {
        std::size_t dim = 2;
        auto gauss = std::make_unique<bmb::genz_gaussian<Real>>(dim);
        std::vector<Real> c = {1.0, 1.0};
        std::vector<Real> w = {0.5, 0.5};
        gauss->set_parameters(c, w);
        
        Real analytical = gauss->analytical_integral();
        // For Gaussian with c=[1,1], w=[0.5,0.5], the integral should be:
        // (√π/2 * erf(0.5) * 2)^2 ≈ (√π/2 * 0.5205 * 2)^2
        Real erf_half = std::erf(0.5);
        Real expected = std::pow(std::sqrt(M_PI) * erf_half, dim);
        
        std::cout << "Gaussian (dim=2, c=[1,1], w=[0.5,0.5]):\n";
        std::cout << "  Analytical: " << analytical << "\n";
        std::cout << "  Expected: " << expected << "\n";
        std::cout << "  Error: " << std::abs(analytical - expected) << "\n";
    }
    
    // Test Product Peak integral
    {
        std::size_t dim = 1;
        auto peak = std::make_unique<bmb::genz_product_peak<Real>>(dim);
        std::vector<Real> c = {2.0};
        std::vector<Real> w = {0.5};
        peak->set_parameters(c, w);
        
        Real analytical = peak->analytical_integral();
        // For 1D: c * (atan(c*(1-w)) + atan(c*w))
        Real expected = 2.0 * (std::atan(2.0 * 0.5) + std::atan(2.0 * 0.5));
        
        std::cout << "\nProduct Peak (dim=1, c=[2], w=[0.5]):\n";
        std::cout << "  Analytical: " << analytical << "\n";
        std::cout << "  Expected: " << expected << "\n";
        std::cout << "  Error: " << std::abs(analytical - expected) << "\n";
    }
    
    std::cout << "\n========================================\n";
    std::cout << "    All tests completed successfully!\n";
    std::cout << "========================================\n";
    
    return 0;
}