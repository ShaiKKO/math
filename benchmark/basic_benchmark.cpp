// Basic benchmark for Boost.Math Cubature Library
// Tests core functionality without complex dependencies

#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>

// Simple timer
class Timer {
    std::chrono::high_resolution_clock::time_point start;
public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    double elapsed() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start).count();
    }
};

// Test functions
double gaussian_2d(double x, double y) {
    return std::exp(-(x*x + y*y));
}

double polynomial_3d(double x, double y, double z) {
    return x*x + y*y + z*z + x*y + y*z;
}

double sphere_indicator_5d(const std::vector<double>& x) {
    double sum = 0;
    for (double xi : x) sum += xi * xi;
    return sum <= 1.0 ? 1.0 : 0.0;
}

// Simple Monte Carlo integration
struct MCResult {
    double value;
    double error;
    size_t samples;
};

MCResult monte_carlo_2d(double (*f)(double, double), 
                        double x_min, double x_max,
                        double y_min, double y_max,
                        size_t n_samples) {
    std::mt19937 gen(12345);
    std::uniform_real_distribution<> dist_x(x_min, x_max);
    std::uniform_real_distribution<> dist_y(y_min, y_max);
    
    double sum = 0, sum_sq = 0;
    double volume = (x_max - x_min) * (y_max - y_min);
    
    for (size_t i = 0; i < n_samples; ++i) {
        double val = f(dist_x(gen), dist_y(gen));
        sum += val;
        sum_sq += val * val;
    }
    
    double mean = sum / n_samples;
    double variance = sum_sq / n_samples - mean * mean;
    
    return {
        volume * mean,
        volume * std::sqrt(variance / n_samples),
        n_samples
    };
}

MCResult monte_carlo_3d(double (*f)(double, double, double),
                        double x_min, double x_max,
                        double y_min, double y_max,
                        double z_min, double z_max,
                        size_t n_samples) {
    std::mt19937 gen(12345);
    std::uniform_real_distribution<> dist_x(x_min, x_max);
    std::uniform_real_distribution<> dist_y(y_min, y_max);
    std::uniform_real_distribution<> dist_z(z_min, z_max);
    
    double sum = 0, sum_sq = 0;
    double volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);
    
    for (size_t i = 0; i < n_samples; ++i) {
        double val = f(dist_x(gen), dist_y(gen), dist_z(gen));
        sum += val;
        sum_sq += val * val;
    }
    
    double mean = sum / n_samples;
    double variance = sum_sq / n_samples - mean * mean;
    
    return {
        volume * mean,
        volume * std::sqrt(variance / n_samples),
        n_samples
    };
}

MCResult monte_carlo_nd(double (*f)(const std::vector<double>&),
                        const std::vector<double>& lower,
                        const std::vector<double>& upper,
                        size_t n_samples) {
    size_t dim = lower.size();
    std::mt19937 gen(12345);
    std::vector<std::uniform_real_distribution<>> dists;
    
    double volume = 1.0;
    for (size_t i = 0; i < dim; ++i) {
        dists.emplace_back(lower[i], upper[i]);
        volume *= (upper[i] - lower[i]);
    }
    
    double sum = 0, sum_sq = 0;
    std::vector<double> point(dim);
    
    for (size_t i = 0; i < n_samples; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            point[j] = dists[j](gen);
        }
        double val = f(point);
        sum += val;
        sum_sq += val * val;
    }
    
    double mean = sum / n_samples;
    double variance = sum_sq / n_samples - mean * mean;
    
    return {
        volume * mean,
        volume * std::sqrt(variance / n_samples),
        n_samples
    };
}

// Simple adaptive integration for 2D (midpoint rule with subdivision)
MCResult adaptive_2d(double (*f)(double, double),
                     double x_min, double x_max,
                     double y_min, double y_max,
                     double tolerance) {
    // Simple implementation - just for comparison
    int n = 10;
    double old_result = 0;
    double result = 0;
    size_t total_evals = 0;
    
    do {
        old_result = result;
        n *= 2;
        double dx = (x_max - x_min) / n;
        double dy = (y_max - y_min) / n;
        
        result = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double x = x_min + (i + 0.5) * dx;
                double y = y_min + (j + 0.5) * dy;
                result += f(x, y) * dx * dy;
                total_evals++;
            }
        }
    } while (std::abs(result - old_result) > tolerance && n < 1000);
    
    return {result, std::abs(result - old_result), total_evals};
}

void print_result(const std::string& method, const MCResult& result, 
                  double expected, double time) {
    double rel_error = std::abs(result.value - expected) / std::abs(expected + 1e-10);
    
    std::cout << std::left << std::setw(20) << method 
              << " | Result: " << std::fixed << std::setprecision(6) << result.value
              << " | Expected: " << expected
              << " | RelErr: " << std::scientific << std::setprecision(2) << rel_error
              << " | EstErr: " << result.error
              << " | Samples: " << std::setw(8) << result.samples
              << " | Time: " << std::fixed << std::setprecision(4) << time << "s"
              << std::endl;
}

int main() {
    std::cout << "\n================================================\n";
    std::cout << "    Basic Numerical Integration Benchmarks\n";
    std::cout << "================================================\n\n";
    
    // Test 1: 2D Gaussian
    {
        std::cout << "Test 1: 2D Gaussian integral on [-2,2]×[-2,2]\n";
        std::cout << "----------------------------------------------\n";
        
        double expected = M_PI * (1 - std::exp(-4));
        
        // Adaptive midpoint
        {
            Timer t;
            auto result = adaptive_2d(gaussian_2d, -2, 2, -2, 2, 1e-4);
            print_result("Adaptive Midpoint", result, expected, t.elapsed());
        }
        
        // Monte Carlo with different sample sizes
        for (size_t n : {1000, 10000, 100000}) {
            Timer t;
            auto result = monte_carlo_2d(gaussian_2d, -2, 2, -2, 2, n);
            print_result("Monte Carlo", result, expected, t.elapsed());
        }
        
        std::cout << std::endl;
    }
    
    // Test 2: 3D Polynomial
    {
        std::cout << "Test 2: 3D Polynomial on [0,1]³\n";
        std::cout << "--------------------------------\n";
        
        double expected = 2.0;
        
        // Monte Carlo with different sample sizes
        for (size_t n : {1000, 10000, 100000}) {
            Timer t;
            auto result = monte_carlo_3d(polynomial_3d, 0, 1, 0, 1, 0, 1, n);
            print_result("Monte Carlo", result, expected, t.elapsed());
        }
        
        std::cout << std::endl;
    }
    
    // Test 3: 5D Unit Sphere
    {
        std::cout << "Test 3: 5D Unit Sphere Volume\n";
        std::cout << "------------------------------\n";
        
        double expected = 8 * M_PI * M_PI / 15.0;  // ~5.2638
        std::vector<double> lower(5, -1);
        std::vector<double> upper(5, 1);
        
        // Monte Carlo with different sample sizes
        for (size_t n : {10000, 100000, 1000000}) {
            Timer t;
            auto result = monte_carlo_nd(sphere_indicator_5d, lower, upper, n);
            print_result("Monte Carlo", result, expected, t.elapsed());
        }
        
        std::cout << std::endl;
    }
    
    // Convergence demonstration
    {
        std::cout << "Monte Carlo Convergence Study (2D Gaussian)\n";
        std::cout << "--------------------------------------------\n";
        
        double expected = M_PI * (1 - std::exp(-4));
        
        std::cout << "N Samples   | Rel. Error  | Std. Error | Error Ratio\n";
        std::cout << "------------|-------------|------------|------------\n";
        
        double prev_error = 1.0;
        for (size_t n = 100; n <= 1000000; n *= 10) {
            auto result = monte_carlo_2d(gaussian_2d, -2, 2, -2, 2, n);
            double rel_error = std::abs(result.value - expected) / expected;
            double ratio = prev_error / rel_error;
            
            std::cout << std::setw(11) << n 
                      << " | " << std::scientific << std::setprecision(3) << rel_error
                      << " | " << result.error
                      << " | " << std::fixed << std::setprecision(2) << ratio
                      << std::endl;
            
            prev_error = rel_error;
        }
        
        std::cout << "\nNote: Error should decrease as ~1/√N (ratio ~3.16 per 10x samples)\n";
    }
    
    std::cout << "\n================================================\n";
    std::cout << "              Benchmark Complete\n";
    std::cout << "================================================\n\n";
    
    return 0;
}