// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/// \file genz_integration_test.cpp
/// \brief Integration tests using Genz test functions
/// \details Tests adaptive, sparse grid, and QMC integrators against
///          the standard Genz test function suite with known exact integrals.

#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/qmc.hpp>
#include "genz_test_functions.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace boost::math::cubature;

/// \brief Test harness for cubature algorithms
template <typename Real>
class cubature_test_harness {
public:
    /// \brief Integrator wrapper for adaptive algorithm
    struct adaptive_integrator {
        std::size_t max_evals;
        
        adaptive_integrator(std::size_t max_evals = 100000) 
            : max_evals(max_evals) {}
        
        template <typename F>
        result<Real> operator()(const F& f, Real abs_tol, Real rel_tol) {
            hypercube<Real> unit_cube(f.dimension());
            std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
            std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
            
            return integrate_adaptive<Real>(f, unit_cube, abs_tol, rel_tol, max_evals);
        }
    };
    
    /// \brief Integrator wrapper for sparse grid algorithm
    struct sparse_grid_integrator {
        unsigned level;
        
        sparse_grid_integrator(unsigned level = 5) : level(level) {}
        
        template <typename F>
        result<Real> operator()(const F& f, Real /*abs_tol*/, Real /*rel_tol*/) {
            hypercube<Real> unit_cube(f.dimension());
            std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
            std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
            
            return integrate_sparse_grid<Real>(f, unit_cube, level);
        }
    };
    
    /// \brief Integrator wrapper for QMC algorithm
    struct qmc_integrator {
        std::size_t n_points;
        
        qmc_integrator(std::size_t n_points = 10000) : n_points(n_points) {}
        
        template <typename F>
        result<Real> operator()(const F& f, Real /*abs_tol*/, Real /*rel_tol*/) {
            hypercube<Real> unit_cube(f.dimension());
            std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
            std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
            
            return integrate_qmc<Real>(f, unit_cube, n_points);
        }
    };
    
    /// \brief Print test results in tabular format
    static void print_results(
        const std::string& algorithm_name,
        const std::vector<test::test_result<Real>>& results)
    {
        std::cout << "\n=== " << algorithm_name << " Results ===\n\n";
        std::cout << std::setw(15) << "Function"
                  << std::setw(6) << "Dim"
                  << std::setw(15) << "Exact"
                  << std::setw(15) << "Computed"
                  << std::setw(15) << "Abs Error"
                  << std::setw(15) << "Rel Error"
                  << std::setw(10) << "Evals"
                  << std::setw(10) << "Status"
                  << "\n";
        std::cout << std::string(110, '-') << "\n";
        
        for (const auto& r : results) {
            std::cout << std::setw(15) << r.function_name
                      << std::setw(6) << r.dimension
                      << std::setw(15) << std::scientific << std::setprecision(6) 
                      << r.exact_value
                      << std::setw(15) << r.computed_value
                      << std::setw(15) << r.absolute_error
                      << std::setw(15) << r.relative_error
                      << std::setw(10) << std::fixed << r.evaluations
                      << std::setw(10) << (r.passed ? "PASS" : "FAIL")
                      << "\n";
        }
        
        // Summary statistics
        int passed = 0;
        for (const auto& r : results) {
            if (r.passed) passed++;
        }
        
        std::cout << "\nSummary: " << passed << "/" << results.size() 
                  << " tests passed\n";
    }
    
    /// \brief Run comprehensive test suite
    static void run_comprehensive_tests() {
        const Real abs_tol = 1e-6;
        const Real rel_tol = 1e-6;
        
        // Test dimensions
        std::vector<std::size_t> dimensions = {2, 3, 4, 5};
        
        // Test difficulties
        std::vector<int> difficulties = {0, 1, 2};  // easy, moderate, hard
        
        for (std::size_t dim : dimensions) {
            std::cout << "\n" << std::string(50, '=') << "\n";
            std::cout << "Testing Dimension " << dim << "\n";
            std::cout << std::string(50, '=') << "\n";
            
            for (int diff : difficulties) {
                std::cout << "\nDifficulty: " 
                          << (diff == 0 ? "Easy" : diff == 1 ? "Moderate" : "Hard")
                          << "\n";
                
                // Test adaptive integration
                {
                    adaptive_integrator integrator(100000);
                    auto results = test::genz_test_suite<Real>::run_suite(
                        dim, integrator, abs_tol, rel_tol, diff);
                    print_results("Adaptive", results);
                }
                
                // Test sparse grid (only for moderate difficulty)
                if (diff == 1 && dim <= 4) {
                    sparse_grid_integrator integrator(5);
                    auto results = test::genz_test_suite<Real>::run_suite(
                        dim, integrator, abs_tol, rel_tol, diff);
                    print_results("Sparse Grid", results);
                }
                
                // Test QMC (only for moderate difficulty)
                if (diff == 1 && dim <= 10) {
                    qmc_integrator integrator(10000);
                    auto results = test::genz_test_suite<Real>::run_suite(
                        dim, integrator, abs_tol, rel_tol, diff);
                    print_results("QMC", results);
                }
            }
        }
    }
    
    /// \brief Run convergence study for a specific function
    static void run_convergence_study() {
        std::cout << "\n" << std::string(50, '=') << "\n";
        std::cout << "Convergence Study: Gaussian Function\n";
        std::cout << std::string(50, '=') << "\n\n";
        
        const std::size_t dim = 3;
        std::vector<Real> a = {5.0, 5.0, 5.0};
        std::vector<Real> u = {0.5, 0.5, 0.5};
        
        test::genz_gaussian<Real> f(dim, a, u);
        Real exact = f.exact_integral();
        
        std::cout << "Exact integral: " << exact << "\n\n";
        
        // Adaptive convergence
        std::cout << "Adaptive Integration Convergence:\n";
        std::cout << std::setw(10) << "Max Evals"
                  << std::setw(15) << "Result"
                  << std::setw(15) << "Error"
                  << std::setw(10) << "Used"
                  << "\n";
        
        for (std::size_t max_evals : {100, 500, 1000, 5000, 10000, 50000}) {
            hypercube<Real> unit_cube(dim);
            std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
            std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
            
            auto result = integrate_adaptive<Real>(f, unit_cube, 1e-10, 1e-10, max_evals);
            Real error = std::abs(result.value - exact);
            
            std::cout << std::setw(10) << max_evals
                      << std::setw(15) << std::scientific << std::setprecision(8) 
                      << result.value
                      << std::setw(15) << error
                      << std::setw(10) << std::fixed << result.evaluations
                      << "\n";
        }
        
        // Sparse grid convergence
        std::cout << "\nSparse Grid Convergence:\n";
        std::cout << std::setw(10) << "Level"
                  << std::setw(15) << "Result"
                  << std::setw(15) << "Error"
                  << std::setw(10) << "Nodes"
                  << "\n";
        
        for (unsigned level = 2; level <= 7; ++level) {
            hypercube<Real> unit_cube(dim);
            std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
            std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
            
            auto result = integrate_sparse_grid<Real>(f, unit_cube, level);
            Real error = std::abs(result.value - exact);
            
            std::cout << std::setw(10) << level
                      << std::setw(15) << std::scientific << std::setprecision(8) 
                      << result.value
                      << std::setw(15) << error
                      << std::setw(10) << std::fixed << result.evaluations
                      << "\n";
        }
        
        // QMC convergence
        std::cout << "\nQMC Convergence:\n";
        std::cout << std::setw(10) << "Points"
                  << std::setw(15) << "Result"
                  << std::setw(15) << "Error"
                  << "\n";
        
        for (std::size_t n_points : {100, 500, 1000, 5000, 10000, 50000}) {
            hypercube<Real> unit_cube(dim);
            std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
            std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
            
            auto result = integrate_qmc<Real>(f, unit_cube, n_points);
            Real error = std::abs(result.value - exact);
            
            std::cout << std::setw(10) << n_points
                      << std::setw(15) << std::scientific << std::setprecision(8) 
                      << result.value
                      << std::setw(15) << error
                      << "\n";
        }
    }
};

int main() {
    std::cout << "Boost.Math Cubature Library - Genz Function Test Suite\n";
    std::cout << "========================================================\n";
    
    // Run basic tests
    cubature_test_harness<double>::run_comprehensive_tests();
    
    // Run convergence study
    cubature_test_harness<double>::run_convergence_study();
    
    return 0;
}