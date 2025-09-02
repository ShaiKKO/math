// Working benchmark for Boost.Math Cubature Library
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

// Test functions
template <typename Real>
Real gaussian_2d(const Real* x, std::size_t) {
    return std::exp(-(x[0]*x[0] + x[1]*x[1]));
}

template <typename Real>
Real polynomial_3d(const Real* x, std::size_t) {
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[0]*x[1] + x[1]*x[2];
}

template <typename Real>
Real unit_sphere_5d(const Real* x, std::size_t) {
    Real sum = 0;
    for (int i = 0; i < 5; ++i) {
        sum += x[i] * x[i];
    }
    return sum <= 1.0 ? 1.0 : 0.0;
}

template <typename Real>
Real exponential_10d(const Real* x, std::size_t) {
    Real sum = 0;
    for (int i = 0; i < 10; ++i) {
        sum += x[i];
    }
    return std::exp(-sum);
}

// Timer utility
class Timer {
    std::chrono::high_resolution_clock::time_point start;
public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    double elapsed() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start).count();
    }
};

void print_result(const std::string& method, double result, double expected, 
                  double time, double error = 0) {
    double rel_error = std::abs(result - expected) / std::abs(expected + 1e-10);
    
    std::cout << std::left << std::setw(15) << method 
              << " | Result: " << std::fixed << std::setprecision(6) << result
              << " | Expected: " << expected
              << " | RelErr: " << std::scientific << std::setprecision(2) << rel_error
              << " | Time: " << std::fixed << std::setprecision(4) << time << "s";
    
    if (error > 0) {
        std::cout << " | EstErr: " << std::scientific << std::setprecision(2) << error;
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "\n============================================\n";
    std::cout << "   Boost.Math Cubature Library Benchmarks\n";
    std::cout << "============================================\n\n";
    
    // Test 1: 2D Gaussian
    {
        std::cout << "Test 1: 2D Gaussian integral on [-2,2]×[-2,2]\n";
        std::cout << "----------------------------------------------\n";
        
        hypercube<double> box(2);
        box.lower[0] = -2; box.lower[1] = -2;
        box.upper[0] = 2; box.upper[1] = 2;
        
        double expected = M_PI * (1 - std::exp(-4));
        
        // Adaptive
        {
            Timer t;
            auto [result, error] = integrate_adaptive(gaussian_2d<double>, box, 1e-6, 0);
            print_result("Adaptive", result, expected, t.elapsed(), error);
        }
        
        // Sparse Grid
        {
            Timer t;
            auto [result, error] = integrate_sparse_grid(gaussian_2d<double>, box, 5);
            print_result("Sparse Grid", result, expected, t.elapsed(), error);
        }
        
        // QMC
        {
            Timer t;
            auto [result, error] = integrate_qmc(gaussian_2d<double>, box, 10000);
            print_result("QMC (10k pts)", result, expected, t.elapsed(), error);
        }
        
        // Monte Carlo
        {
            Timer t;
            std::vector<double> lower = {-2, -2};
            std::vector<double> upper = {2, 2};
            auto mc_result = monte_carlo(gaussian_2d<double>, lower, upper, 10000);
            print_result("Monte Carlo", mc_result.value, expected, t.elapsed(), mc_result.error);
        }
        
        std::cout << std::endl;
    }
    
    // Test 2: 3D Polynomial
    {
        std::cout << "Test 2: 3D Polynomial on [0,1]³\n";
        std::cout << "--------------------------------\n";
        
        hypercube<double> box(3);
        for (int i = 0; i < 3; ++i) {
            box.lower[i] = 0;
            box.upper[i] = 1;
        }
        
        double expected = 2.0;
        
        // Adaptive
        {
            Timer t;
            auto [result, error] = integrate_adaptive(polynomial_3d<double>, box, 1e-8, 0);
            print_result("Adaptive", result, expected, t.elapsed(), error);
        }
        
        // Sparse Grid
        {
            Timer t;
            auto [result, error] = integrate_sparse_grid(polynomial_3d<double>, box, 4);
            print_result("Sparse Grid", result, expected, t.elapsed(), error);
        }
        
        // QMC
        {
            Timer t;
            auto [result, error] = integrate_qmc(polynomial_3d<double>, box, 5000);
            print_result("QMC (5k pts)", result, expected, t.elapsed(), error);
        }
        
        // Monte Carlo
        {
            Timer t;
            std::vector<double> lower = {0, 0, 0};
            std::vector<double> upper = {1, 1, 1};
            auto mc_result = monte_carlo(polynomial_3d<double>, lower, upper, 10000);
            print_result("Monte Carlo", mc_result.value, expected, t.elapsed(), mc_result.error);
        }
        
        std::cout << std::endl;
    }
    
    // Test 3: 5D Unit Sphere
    {
        std::cout << "Test 3: 5D Unit Sphere Volume\n";
        std::cout << "------------------------------\n";
        
        hypercube<double> box(5);
        for (int i = 0; i < 5; ++i) {
            box.lower[i] = -1;
            box.upper[i] = 1;
        }
        
        double expected = 8 * M_PI * M_PI / 15.0;
        
        // Adaptive (may be slow)
        {
            Timer t;
            auto [result, error] = integrate_adaptive(unit_sphere_5d<double>, box, 1e-2, 0);
            print_result("Adaptive", result, expected, t.elapsed(), error);
        }
        
        // Sparse Grid
        {
            Timer t;
            auto [result, error] = integrate_sparse_grid(unit_sphere_5d<double>, box, 3);
            print_result("Sparse Grid", result, expected, t.elapsed(), error);
        }
        
        // QMC
        {
            Timer t;
            auto [result, error] = integrate_qmc(unit_sphere_5d<double>, box, 50000);
            print_result("QMC (50k pts)", result, expected, t.elapsed(), error);
        }
        
        // Monte Carlo
        {
            Timer t;
            std::vector<double> lower(5, -1);
            std::vector<double> upper(5, 1);
            auto mc_result = monte_carlo(unit_sphere_5d<double>, lower, upper, 100000);
            print_result("Monte Carlo", mc_result.value, expected, t.elapsed(), mc_result.error);
        }
        
        std::cout << std::endl;
    }
    
    // Test 4: 10D Exponential (high dimension test)
    {
        std::cout << "Test 4: 10D Exponential on [0,1]¹⁰\n";
        std::cout << "-----------------------------------\n";
        
        hypercube<double> box(10);
        for (int i = 0; i < 10; ++i) {
            box.lower[i] = 0;
            box.upper[i] = 1;
        }
        
        double expected = std::pow(1 - std::exp(-1), 10);
        
        // Sparse Grid (adaptive would be too slow)
        {
            Timer t;
            auto [result, error] = integrate_sparse_grid(exponential_10d<double>, box, 2);
            print_result("Sparse Grid", result, expected, t.elapsed(), error);
        }
        
        // QMC - should perform well
        {
            Timer t;
            auto [result, error] = integrate_qmc(exponential_10d<double>, box, 100000);
            print_result("QMC (100k pts)", result, expected, t.elapsed(), error);
        }
        
        // Monte Carlo
        {
            Timer t;
            std::vector<double> lower(10, 0);
            std::vector<double> upper(10, 1);
            auto mc_result = monte_carlo(exponential_10d<double>, lower, upper, 100000);
            print_result("Monte Carlo", mc_result.value, expected, t.elapsed(), mc_result.error);
        }
        
        // Monte Carlo with variance reduction
        {
            Timer t;
            std::vector<double> lower(10, 0);
            std::vector<double> upper(10, 1);
            auto mc_result = integrate_monte_carlo(exponential_10d<double>, lower, upper, 1e-3, 1e-2);
            print_result("MC Adaptive", mc_result.value, expected, t.elapsed(), mc_result.error);
        }
        
        std::cout << std::endl;
    }
    
    // Summary
    std::cout << "============================================\n";
    std::cout << "                 Summary\n";
    std::cout << "============================================\n";
    std::cout << "• Adaptive: Best for low dimensions (d ≤ 5)\n";
    std::cout << "• Sparse Grid: Good for moderate dimensions\n";
    std::cout << "• QMC: Excellent for high dimensions\n";
    std::cout << "• Monte Carlo: Dimension-independent convergence\n";
    std::cout << "============================================\n\n";
    
    return 0;
}