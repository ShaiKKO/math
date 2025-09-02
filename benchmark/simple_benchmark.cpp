// Simple benchmark to test our cubature library methods
#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <iomanip>

// Include our cubature headers
#include "../include/boost/math/cubature/adaptive.hpp"
#include "../include/boost/math/cubature/sparse_grid.hpp"
#include "../include/boost/math/cubature/qmc.hpp"
#include "../include/boost/math/cubature/monte_carlo.hpp"

using namespace boost::math::cubature;

// Simple test functions
template <typename Real>
Real gaussian_2d(const Real* x, std::size_t) {
    return std::exp(-(x[0]*x[0] + x[1]*x[1]));
}

template <typename Real>
Real polynomial_3d(const Real* x, std::size_t) {
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[0]*x[1] + x[1]*x[2];
}

template <typename Real>
Real oscillatory_2d(const Real* x, std::size_t) {
    return std::cos(10*x[0]) * std::sin(10*x[1]);
}

template <typename Real>
Real unit_sphere_5d(const Real* x, std::size_t) {
    Real sum = 0;
    for (int i = 0; i < 5; ++i) {
        sum += x[i] * x[i];
    }
    return sum <= 1.0 ? 1.0 : 0.0;
}

// Benchmark a single method
template <typename Func>
void benchmark_method(const std::string& method_name,
                      const std::string& function_name,
                      Func integrator,
                      double expected) {
    auto start = std::chrono::high_resolution_clock::now();
    
    auto [result, error] = integrator();
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    double relative_error = std::abs(result - expected) / std::abs(expected);
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << method_name << " | " 
              << function_name << " | "
              << "Result: " << result << " | "
              << "Expected: " << expected << " | "
              << "Rel Error: " << std::scientific << relative_error << " | "
              << "Time: " << std::fixed << elapsed.count() << "s\n";
}

int main() {
    std::cout << "\n=== Boost.Math Cubature Library Benchmarks ===\n\n";
    
    // Test 1: 2D Gaussian integral
    std::cout << "Test 1: 2D Gaussian on [-2,2]×[-2,2]\n";
    std::cout << "----------------------------------------\n";
    {
        std::vector<double> lower = {-2, -2};
        std::vector<double> upper = {2, 2};
        double expected = M_PI * (1 - std::exp(-4));  // Analytical value
        double tol = 1e-4;
        
        // Adaptive integration
        benchmark_method("Adaptive", "Gaussian 2D",
            [&]() { return integrate_adaptive(gaussian_2d<double>, lower, upper, tol, 0); },
            expected);
        
        // Sparse grid
        benchmark_method("Sparse Grid", "Gaussian 2D",
            [&]() { 
                hypercube<double> box(lower, upper);
                return integrate_sparse_grid(gaussian_2d<double>, box, 5);
            },
            expected);
        
        // QMC (Sobol)
        benchmark_method("QMC (Sobol)", "Gaussian 2D",
            [&]() { 
                hypercube<double> box(lower, upper);
                return integrate_qmc(gaussian_2d<double>, box, 10000);
            },
            expected);
        
        // Monte Carlo
        benchmark_method("Monte Carlo", "Gaussian 2D",
            [&]() { 
                return monte_carlo(gaussian_2d<double>, lower, upper, 10000);
            },
            expected);
    }
    
    std::cout << "\n";
    
    // Test 2: 3D Polynomial
    std::cout << "Test 2: 3D Polynomial on [0,1]³\n";
    std::cout << "----------------------------------------\n";
    {
        std::vector<double> lower = {0, 0, 0};
        std::vector<double> upper = {1, 1, 1};
        double expected = 2.0;  // Analytical value for this polynomial
        double tol = 1e-6;
        
        // Adaptive integration
        benchmark_method("Adaptive", "Polynomial 3D",
            [&]() { return integrate_adaptive(polynomial_3d<double>, lower, upper, tol, 0); },
            expected);
        
        // Sparse grid
        benchmark_method("Sparse Grid", "Polynomial 3D",
            [&]() { 
                hypercube<double> box(lower, upper);
                return integrate_sparse_grid(polynomial_3d<double>, box, 4);
            },
            expected);
        
        // QMC
        benchmark_method("QMC (Sobol)", "Polynomial 3D",
            [&]() { 
                hypercube<double> box(lower, upper);
                return integrate_qmc(polynomial_3d<double>, box, 5000);
            },
            expected);
        
        // Monte Carlo
        benchmark_method("Monte Carlo", "Polynomial 3D",
            [&]() { 
                return monte_carlo(polynomial_3d<double>, lower, upper, 10000);
            },
            expected);
    }
    
    std::cout << "\n";
    
    // Test 3: 2D Oscillatory function
    std::cout << "Test 3: 2D Oscillatory on [0,1]²\n";
    std::cout << "----------------------------------------\n";
    {
        std::vector<double> lower = {0, 0};
        std::vector<double> upper = {1, 1};
        double expected = 0.0;  // Due to symmetry
        double tol = 1e-4;
        
        // Adaptive integration
        benchmark_method("Adaptive", "Oscillatory 2D",
            [&]() { return integrate_adaptive(oscillatory_2d<double>, lower, upper, tol, 0); },
            expected);
        
        // Sparse grid
        benchmark_method("Sparse Grid", "Oscillatory 2D",
            [&]() { 
                hypercube<double> box(lower, upper);
                return integrate_sparse_grid(oscillatory_2d<double>, box, 7);
            },
            expected);
        
        // QMC
        benchmark_method("QMC (Sobol)", "Oscillatory 2D",
            [&]() { 
                hypercube<double> box(lower, upper);
                return integrate_qmc(oscillatory_2d<double>, box, 20000);
            },
            expected);
        
        // Monte Carlo
        benchmark_method("Monte Carlo", "Oscillatory 2D",
            [&]() { 
                return monte_carlo(oscillatory_2d<double>, lower, upper, 50000);
            },
            expected);
    }
    
    std::cout << "\n";
    
    // Test 4: 5D Unit sphere volume
    std::cout << "Test 4: 5D Unit Sphere Volume\n";
    std::cout << "----------------------------------------\n";
    {
        std::vector<double> lower = {-1, -1, -1, -1, -1};
        std::vector<double> upper = {1, 1, 1, 1, 1};
        double expected = 8 * M_PI * M_PI / 15.0;  // Volume of 5D unit sphere
        double tol = 1e-2;
        
        // Adaptive integration (might be slow)
        benchmark_method("Adaptive", "Sphere 5D",
            [&]() { return integrate_adaptive(unit_sphere_5d<double>, lower, upper, tol, 0); },
            expected);
        
        // Sparse grid
        benchmark_method("Sparse Grid", "Sphere 5D",
            [&]() { 
                hypercube<double> box(lower, upper);
                return integrate_sparse_grid(unit_sphere_5d<double>, box, 3);
            },
            expected);
        
        // QMC
        benchmark_method("QMC (Sobol)", "Sphere 5D",
            [&]() { 
                hypercube<double> box(lower, upper);
                return integrate_qmc(unit_sphere_5d<double>, box, 50000);
            },
            expected);
        
        // Monte Carlo
        benchmark_method("Monte Carlo", "Sphere 5D",
            [&]() { 
                return monte_carlo(unit_sphere_5d<double>, lower, upper, 100000);
            },
            expected);
    }
    
    std::cout << "\n=== Benchmark Complete ===\n\n";
    
    // Summary of methods
    std::cout << "Method Characteristics:\n";
    std::cout << "----------------------\n";
    std::cout << "Adaptive: Best for low dimensions (d ≤ 5), smooth functions\n";
    std::cout << "Sparse Grid: Good for moderate dimensions (5 < d ≤ 20)\n";
    std::cout << "QMC (Sobol): Excellent for high dimensions, better than MC\n";
    std::cout << "Monte Carlo: Dimension-independent, good with variance reduction\n";
    
    return 0;
}