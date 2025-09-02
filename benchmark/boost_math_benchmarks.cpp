// Copyright (c) 2025 Boost.Math Contributors
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include "benchmark_base.hpp"
#include "test_functions.hpp"

// Include Boost.Math cubature headers
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/qmc.hpp>
#include <boost/math/cubature/monte_carlo.hpp>

#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <memory>

namespace bmb = boost::math::benchmark;
namespace bmc = boost::math::cubature;

template <typename Real>
struct integrator_result {
    Real value;
    Real error;
    std::size_t evaluations;
};

// Wrapper for adaptive integration
template <typename Real, typename Function>
integrator_result<Real> benchmark_adaptive(Function& f, 
                                          std::size_t dimension,
                                          Real tolerance) {
    integrator_result<Real> result;
    
    // Create bounds for unit hypercube
    std::vector<Real> lower(dimension, 0);
    std::vector<Real> upper(dimension, 1);
    
    // Reset function evaluation counter
    f.reset_eval_count();
    
    try {
        // Use adaptive integration
        auto [value, error] = bmc::adaptive(
            [&f](const std::vector<Real>& x) { return f(x); },
            lower, upper,
            bmc::tolerance<Real>(tolerance, 0)  // relative tolerance only
        );
        
        result.value = value;
        result.error = error;
        result.evaluations = f.get_eval_count();
    } catch (const std::exception& e) {
        std::cerr << "Adaptive integration failed: " << e.what() << "\n";
        throw;
    }
    
    return result;
}

// Wrapper for sparse grid integration
template <typename Real, typename Function>
integrator_result<Real> benchmark_sparse_grid(Function& f,
                                             std::size_t dimension,
                                             Real tolerance) {
    integrator_result<Real> result;
    
    // Reset function evaluation counter
    f.reset_eval_count();
    
    try {
        // Determine appropriate level based on tolerance
        std::size_t level = 3;  // Start with level 3
        if (tolerance < 1e-3) level = 4;
        if (tolerance < 1e-4) level = 5;
        if (tolerance < 1e-6) level = 6;
        
        // Use Clenshaw-Curtis sparse grid
        auto sg = bmc::make_sparse_grid<Real>(
            dimension, level,
            bmc::sparse_grid_rule::clenshaw_curtis
        );
        
        // Integrate over unit hypercube
        auto [value, error] = sg.integrate(
            [&f](const std::vector<Real>& x) { return f(x); }
        );
        
        result.value = value;
        result.error = error;
        result.evaluations = f.get_eval_count();
    } catch (const std::exception& e) {
        std::cerr << "Sparse grid integration failed: " << e.what() << "\n";
        throw;
    }
    
    return result;
}

// Wrapper for QMC integration
template <typename Real, typename Function>
integrator_result<Real> benchmark_qmc(Function& f,
                                     std::size_t dimension,
                                     Real tolerance) {
    integrator_result<Real> result;
    
    // Reset function evaluation counter
    f.reset_eval_count();
    
    try {
        // Determine number of points based on tolerance
        std::size_t num_points = 1000;
        if (tolerance < 1e-3) num_points = 10000;
        if (tolerance < 1e-4) num_points = 100000;
        if (tolerance < 1e-6) num_points = 1000000;
        
        // Use Sobol sequence
        auto qmc = bmc::qmc<Real>(dimension, bmc::qmc_sequence::sobol);
        
        // Integrate over unit hypercube
        auto [value, error] = qmc.integrate(
            [&f](const std::vector<Real>& x) { return f(x); },
            num_points
        );
        
        result.value = value;
        result.error = error;
        result.evaluations = f.get_eval_count();
    } catch (const std::exception& e) {
        std::cerr << "QMC integration failed: " << e.what() << "\n";
        throw;
    }
    
    return result;
}

// Wrapper for Monte Carlo integration
template <typename Real, typename Function>
integrator_result<Real> benchmark_monte_carlo(Function& f,
                                             std::size_t dimension,
                                             Real tolerance) {
    integrator_result<Real> result;
    
    // Reset function evaluation counter
    f.reset_eval_count();
    
    try {
        // Determine number of samples based on tolerance
        std::size_t num_samples = 10000;
        if (tolerance < 1e-2) num_samples = 100000;
        if (tolerance < 1e-3) num_samples = 1000000;
        if (tolerance < 1e-4) num_samples = 10000000;
        
        // Use standard Monte Carlo
        auto mc = bmc::monte_carlo<Real>(dimension);
        
        // Integrate over unit hypercube
        auto [value, error] = mc.integrate(
            [&f](const std::vector<Real>& x) { return f(x); },
            num_samples
        );
        
        result.value = value;
        result.error = error;
        result.evaluations = f.get_eval_count();
    } catch (const std::exception& e) {
        std::cerr << "Monte Carlo integration failed: " << e.what() << "\n";
        throw;
    }
    
    return result;
}

// Main benchmark function for a single configuration
template <typename Real>
bmb::benchmark_result run_single_benchmark(
    const std::string& method_name,
    const std::string& function_name,
    std::size_t dimension,
    Real tolerance,
    std::unique_ptr<bmb::genz_function_base<Real>>& test_func,
    std::size_t num_runs = 3) {
    
    bmb::benchmark_result result;
    result.method_name = method_name;
    result.function_name = function_name;
    result.dimension = dimension;
    result.tolerance = tolerance;
    
    // Get analytical solution
    Real analytical = 0;
    if (auto* osc = dynamic_cast<bmb::genz_oscillatory<Real>*>(test_func.get())) {
        analytical = osc->analytical_integral();
    } else if (auto* peak = dynamic_cast<bmb::genz_product_peak<Real>*>(test_func.get())) {
        analytical = peak->analytical_integral();
    } else if (auto* corner = dynamic_cast<bmb::genz_corner_peak<Real>*>(test_func.get())) {
        analytical = corner->analytical_integral();
    } else if (auto* gauss = dynamic_cast<bmb::genz_gaussian<Real>*>(test_func.get())) {
        analytical = gauss->analytical_integral();
    } else if (auto* cont = dynamic_cast<bmb::genz_continuous<Real>*>(test_func.get())) {
        analytical = cont->analytical_integral();
    } else if (auto* disc = dynamic_cast<bmb::genz_discontinuous<Real>*>(test_func.get())) {
        analytical = disc->analytical_integral();
    }
    
    std::vector<double> times;
    std::vector<double> errors;
    std::vector<std::size_t> calls;
    
    for (std::size_t run = 0; run < num_runs; ++run) {
        bmb::high_resolution_timer timer;
        integrator_result<Real> int_result;
        
        try {
            timer.start();
            
            if (method_name == "adaptive") {
                int_result = benchmark_adaptive<Real>(*test_func, dimension, tolerance);
            } else if (method_name == "sparse_grid") {
                int_result = benchmark_sparse_grid<Real>(*test_func, dimension, tolerance);
            } else if (method_name == "qmc") {
                int_result = benchmark_qmc<Real>(*test_func, dimension, tolerance);
            } else if (method_name == "monte_carlo") {
                int_result = benchmark_monte_carlo<Real>(*test_func, dimension, tolerance);
            } else {
                throw std::runtime_error("Unknown method: " + method_name);
            }
            
            timer.stop();
            
            Real rel_error = std::abs(int_result.value - analytical) / 
                           (std::abs(analytical) + 1e-10);
            
            times.push_back(timer.elapsed_seconds());
            errors.push_back(rel_error);
            calls.push_back(int_result.evaluations);
            
        } catch (const std::exception& e) {
            result.converged = false;
            result.error_message = e.what();
            break;
        }
    }
    
    if (result.converged && !times.empty()) {
        // Use median values
        std::sort(times.begin(), times.end());
        std::sort(errors.begin(), errors.end());
        std::sort(calls.begin(), calls.end());
        
        size_t mid = times.size() / 2;
        result.execution_time = times[mid];
        result.relative_error = errors[mid];
        result.function_calls = calls[mid];
        
        // Calculate efficiency score
        if (result.relative_error > 0 && result.execution_time > 0) {
            result.efficiency_score = 1.0 / (result.execution_time * result.relative_error);
        }
    }
    
    return result;
}

// Run comprehensive benchmarks
template <typename Real>
void run_boost_math_benchmarks(bmb::benchmark_suite& suite) {
    // Test configurations
    std::vector<std::size_t> dimensions = {2, 3, 5, 8, 10, 15, 20, 30};
    std::vector<Real> tolerances = {1e-2, 1e-3, 1e-4, 1e-6};
    std::vector<std::string> methods = {"adaptive", "sparse_grid", "qmc", "monte_carlo"};
    
    // All Genz test functions
    auto function_types = bmb::genz_test_suite<Real>::all_functions();
    
    std::size_t total_benchmarks = dimensions.size() * tolerances.size() * 
                                  methods.size() * function_types.size();
    std::size_t current = 0;
    
    std::cout << "Running Boost.Math cubature benchmarks...\n";
    std::cout << "Total configurations: " << total_benchmarks << "\n\n";
    
    for (auto func_type : function_types) {
        const char* func_name = bmb::genz_test_suite<Real>::function_name(func_type);
        
        for (std::size_t dim : dimensions) {
            // Skip very high dimensions for some methods
            if (dim > 20 && (func_type == bmb::genz_test_suite<Real>::function_type::oscillatory ||
                           func_type == bmb::genz_test_suite<Real>::function_type::discontinuous)) {
                continue;  // These are particularly challenging in high dimensions
            }
            
            for (Real tol : tolerances) {
                // Skip very tight tolerances in high dimensions for Monte Carlo
                if (dim > 15 && tol < 1e-4) {
                    continue;  // Would take too long
                }
                
                for (const auto& method : methods) {
                    ++current;
                    
                    // Skip certain combinations that would take too long
                    if (method == "monte_carlo" && tol < 1e-4 && dim > 10) {
                        continue;
                    }
                    if (method == "adaptive" && dim > 15) {
                        continue;  // Adaptive struggles in very high dimensions
                    }
                    if (method == "sparse_grid" && dim > 20) {
                        continue;  // Sparse grid becomes impractical
                    }
                    
                    std::cout << "Progress: " << current << "/" << total_benchmarks 
                              << " - " << method << " on " << func_name 
                              << " (dim=" << dim << ", tol=" << tol << ")...";
                    std::cout.flush();
                    
                    // Create test function
                    auto test_func = bmb::genz_test_suite<Real>::create_function(func_type, dim);
                    
                    // Run benchmark
                    try {
                        auto result = run_single_benchmark<Real>(
                            method, func_name, dim, tol, test_func, 3
                        );
                        
                        suite.add_result(result);
                        
                        std::cout << " done (error=" << std::scientific 
                                  << std::setprecision(2) << result.relative_error 
                                  << ", time=" << std::fixed << std::setprecision(3) 
                                  << result.execution_time << "s)\n";
                    } catch (const std::exception& e) {
                        std::cout << " FAILED: " << e.what() << "\n";
                    }
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    try {
        // Create benchmark suite
        bmb::benchmark_suite suite("Boost.Math Cubature Benchmarks");
        
        // Run benchmarks with double precision
        run_boost_math_benchmarks<double>(suite);
        
        // Export results
        std::string output_prefix = "boost_math_benchmark";
        if (argc > 1) {
            output_prefix = argv[1];
        }
        
        suite.export_csv(output_prefix + ".csv");
        suite.export_json(output_prefix + ".json");
        
        // Print summary
        suite.print_summary();
        
        // Analyze dimension scaling for each method
        std::cout << "\n\n=== DIMENSION SCALING ANALYSIS ===\n";
        for (const auto& method : {"adaptive", "sparse_grid", "qmc", "monte_carlo"}) {
            suite.analyze_dimension_scaling(method, "gaussian", 1e-3);
        }
        
        std::cout << "\nBenchmark results saved to:\n";
        std::cout << "  - " << output_prefix << ".csv\n";
        std::cout << "  - " << output_prefix << ".json\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}