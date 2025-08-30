// Copyright 2025 Your Name
// Test for parallel executor with Boost.Asio

#include <boost/math/cubature/detail/parallel_executor.hpp>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>

using namespace boost::math::cubature::detail;

// Simple test function
double compute_heavy(int x) {
    double result = 0;
    for (int i = 0; i < 1000; ++i) {
        result += std::sin(x * 0.001 * i) * std::cos(x * 0.002 * i);
    }
    return result;
}

int main() {
    std::cout << "Testing Parallel Executor with Boost.Asio\n";
    std::cout << "==========================================\n\n";
    
    // Test data
    std::vector<int> data(100);
    std::iota(data.begin(), data.end(), 0);
    
    // Test with different thread counts
    std::vector<std::size_t> thread_counts = {1, 2, 4};
    std::vector<double> results;
    
    for (std::size_t n_threads : thread_counts) {
        parallel_config config;
        config.num_threads = n_threads;
        config.deterministic = true;
        config.min_batch_size = 10;
        
        parallel_executor<double> executor(config);
        
        // Map-reduce operation
        auto map_func = [](int x) -> double { return compute_heavy(x); };
        auto reduce_op = [](double a, double b) -> double { return a + b; };
        
        double result = executor.map_reduce(data, map_func, reduce_op, 0.0);
        
        std::cout << "Threads: " << n_threads << " -> Result: " << result << "\n";
        
        // Check determinism
        if (!results.empty()) {
            double expected = results[0];
            double diff = std::abs(result - expected);
            if (diff > 1e-10) {
                std::cerr << "ERROR: Results differ by " << diff << "\n";
                return 1;
            }
        }
        results.push_back(result);
    }
    
    std::cout << "\n✓ All thread counts produce identical results\n";
    std::cout << "✓ Parallel executor is deterministic\n";
    
    return 0;
}