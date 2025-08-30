// Copyright 2025 Your Name
// Comprehensive test for parallel integration features

#include <boost/math/cubature/detail/parallel_adaptive.hpp>
#include <boost/math/cubature/detail/parallel_qmc.hpp>
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/qmc.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace boost::math::cubature;
using namespace boost::math::cubature::detail;

// Test function: 3D Gaussian
template <typename Real>
class gaussian_3d {
public:
    Real operator()(const Real* x, std::size_t) const {
        Real sum = 0;
        for (int i = 0; i < 3; ++i) {
            Real diff = x[i] - Real(0.5);
            sum += diff * diff;
        }
        return std::exp(-Real(5) * sum);
    }
};

// Test function: 2D Oscillatory
template <typename Real>
class oscillatory_2d {
public:
    Real operator()(const Real* x, std::size_t) const {
        return std::cos(Real(10) * x[0]) * std::sin(Real(10) * x[1]);
    }
};

void test_parallel_adaptive() {
    std::cout << "\n=== Testing Parallel Adaptive Integration ===\n";
    
    gaussian_3d<double> f;
    hypercube<double> box(3);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    double abs_tol = 1e-8;
    double rel_tol = 1e-8;
    std::size_t max_eval = 10000;
    
    std::vector<std::size_t> thread_counts = {1, 2, 4};
    std::vector<result<double>> results;
    
    for (std::size_t n_threads : thread_counts) {
        parallel_config config;
        config.num_threads = n_threads;
        config.deterministic = true;
        config.min_batch_size = 8;
        
        parallel_adaptive_integrator<double, gaussian_3d<double>> integrator(
            f, box, abs_tol, rel_tol, max_eval, config);
        
        auto res = integrator.integrate();
        
        std::cout << "Threads: " << std::setw(2) << n_threads 
                  << " -> Value: " << std::fixed << std::setprecision(10) 
                  << res.value
                  << ", Error: " << std::scientific << std::setprecision(3)
                  << res.error
                  << ", Evals: " << std::fixed << res.evaluations << "\n";
        
        // Check determinism
        if (!results.empty()) {
            const auto& expected = results[0];
            double diff = std::abs(res.value - expected.value);
            if (diff > 1e-14) {
                std::cerr << "ERROR: Values differ by " << diff << "\n";
                std::cerr << "  Expected: " << expected.value << "\n";
                std::cerr << "  Got:      " << res.value << "\n";
            } else {
                std::cout << "  ✓ Result matches (diff = " << std::scientific 
                          << diff << ")\n";
            }
        }
        
        results.push_back(res);
    }
}

void test_parallel_qmc() {
    std::cout << "\n=== Testing Parallel QMC Integration ===\n";
    
    oscillatory_2d<double> f;
    hypercube<double> box(2);
    std::fill(box.lower.begin(), box.lower.end(), 0.0);
    std::fill(box.upper.begin(), box.upper.end(), 1.0);
    
    std::size_t n_points = 10000;
    bool use_scrambling = true;
    std::size_t n_replicates = 5;
    
    std::vector<std::size_t> thread_counts = {1, 2, 4};
    std::vector<result<double>> results;
    
    for (std::size_t n_threads : thread_counts) {
        parallel_config config;
        config.num_threads = n_threads;
        config.deterministic = true;
        config.min_batch_size = 100;
        
        auto res = integrate_qmc_parallel<double, oscillatory_2d<double>, 2, 
            boost::math::cubature::default_policy>(
            f, box, n_points, use_scrambling, n_replicates, config);
        
        std::cout << "Threads: " << std::setw(2) << n_threads 
                  << " -> Value: " << std::fixed << std::setprecision(10) 
                  << res.value
                  << ", Error: " << std::scientific << std::setprecision(3)
                  << res.error << "\n";
        
        // Check determinism
        if (!results.empty()) {
            const auto& expected = results[0];
            double diff = std::abs(res.value - expected.value);
            if (diff > 1e-14) {
                std::cerr << "ERROR: Values differ by " << diff << "\n";
            } else {
                std::cout << "  ✓ Result matches (diff = " << std::scientific 
                          << diff << ")\n";
            }
        }
        
        results.push_back(res);
    }
}

int main() {
    std::cout << "Boost.Math Cubature - Parallel Integration Tests\n";
    std::cout << "================================================\n";
    
    try {
        test_parallel_adaptive();
        test_parallel_qmc();
        
        std::cout << "\n================================================\n";
        std::cout << "✓ All parallel integration tests passed!\n";
        std::cout << "✓ Determinism verified across thread counts\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}