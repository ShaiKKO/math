// Copyright 2025 Boost Contributors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/cubature/qmc.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <numeric>

using namespace boost::math::cubature;

BOOST_AUTO_TEST_SUITE(test_qmc_sequences)

// Test function 1: Constant function
template <typename Real>
Real constant_function(const Real*, std::size_t) {
    return Real(1);
}

// Test function 2: Linear function
template <typename Real>
Real linear_function(const Real* x, std::size_t dim) {
    Real sum = 0;
    for (std::size_t i = 0; i < dim; ++i) {
        sum += x[i];
    }
    return sum;
}

// Test function 3: Product function
template <typename Real>
Real product_function(const Real* x, std::size_t dim) {
    Real prod = 1;
    for (std::size_t i = 0; i < dim; ++i) {
        prod *= x[i];
    }
    return prod;
}

// Test function 4: Gaussian
template <typename Real>
Real gaussian_function(const Real* x, std::size_t dim) {
    Real sum_sq = 0;
    for (std::size_t i = 0; i < dim; ++i) {
        Real centered = x[i] - Real(0.5);
        sum_sq += centered * centered;
    }
    return std::exp(-Real(10) * sum_sq);
}

// Test function 5: Oscillatory
template <typename Real>
Real oscillatory_function(const Real* x, std::size_t dim) {
    Real sum = 0;
    for (std::size_t i = 0; i < dim; ++i) {
        sum += x[i];
    }
    return std::cos(Real(2) * M_PI * sum);
}

// Helper function to compute analytical integral values
template <typename Real>
Real analytical_integral(const std::string& func_name, std::size_t dim) {
    if (func_name == "constant") {
        return Real(1);  // Integral of 1 over [0,1]^d
    } else if (func_name == "linear") {
        return Real(dim) / Real(2);  // Sum of d integrals of x_i over [0,1]
    } else if (func_name == "product") {
        return std::pow(Real(0.5), Real(dim));  // Product of integrals
    } else if (func_name == "gaussian") {
        // Approximate for centered Gaussian
        Real sigma = Real(1) / std::sqrt(Real(20));
        return std::pow(sigma * std::sqrt(Real(2) * M_PI), Real(dim)) * 
               std::erf(Real(0.5) / (sigma * std::sqrt(Real(2)))) * Real(2);
    }
    return Real(0);
}

// Test all sequence types on a given function
template <typename Real>
void test_sequence_type(detail::qmc_sequence_type seq_type, 
                        const std::string& seq_name,
                        std::size_t dim,
                        std::size_t n_points) {
    
    std::cout << "\nTesting " << seq_name << " sequence in " << dim << "D with " 
              << n_points << " points:" << std::endl;
    
    // Create integration domain [0,1]^d
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Configure QMC options
    qmc_options<Real> opts;
    opts.sequence_type = seq_type;
    opts.use_scrambling = false;
    opts.n_replicates = 1;
    
    // For lattice rules, set the number of points
    if (seq_type == detail::qmc_sequence_type::lattice_rank1) {
        opts.lattice_points = n_points;
    }
    
    // Test different functions
    struct TestCase {
        std::function<Real(const Real*, std::size_t)> func;
        std::string name;
        Real tolerance;
    };
    
    std::vector<TestCase> test_cases = {
        {constant_function<Real>, "constant", Real(1e-10)},
        {linear_function<Real>, "linear", Real(1e-3)},
        {product_function<Real>, "product", Real(1e-2)},
        {gaussian_function<Real>, "gaussian", Real(1e-2)}
    };
    
    for (const auto& test : test_cases) {
        auto result = integrate_qmc_ex<Real>(test.func, box, n_points, opts);
        
        BOOST_CHECK_EQUAL(result.status, status_code::success);
        
        // For simple test functions, get analytical value
        Real exact = analytical_integral<Real>(test.name, dim);
        Real rel_error = std::abs(result.value - exact) / (std::abs(exact) + Real(1e-10));
        
        std::cout << std::fixed << std::setprecision(8)
                  << "  " << test.name << ": value = " << result.value 
                  << ", exact = " << exact
                  << ", rel_error = " << std::scientific << rel_error 
                  << std::endl;
        
        // Check accuracy (relaxed for QMC without many points)
        if (test.name == "constant") {
            // Constant function should be exact for any QMC sequence
            BOOST_CHECK_CLOSE(result.value, exact, Real(1e-8));
        } else if (n_points >= 1000) {
            // For other functions, expect reasonable accuracy with enough points
            BOOST_CHECK_CLOSE(result.value, exact, Real(10)); // 10% tolerance
        }
    }
}

BOOST_AUTO_TEST_CASE(test_all_sequences_2d) {
    using Real = double;
    const std::size_t dim = 2;
    const std::size_t n_points = 1024;
    
    // Test Sobol sequence
    test_sequence_type<Real>(detail::qmc_sequence_type::sobol, "Sobol", dim, n_points);
    
    // Test Halton sequence  
    test_sequence_type<Real>(detail::qmc_sequence_type::halton, "Halton", dim, n_points);
    
    // Test Faure sequence (if available)
    try {
        test_sequence_type<Real>(detail::qmc_sequence_type::faure, "Faure", dim, n_points);
    } catch (const std::runtime_error& e) {
        std::cout << "Faure sequence not available: " << e.what() << std::endl;
    }
    
    // Test Niederreiter sequence (if available)
    try {
        test_sequence_type<Real>(detail::qmc_sequence_type::niederreiter_base2, 
                                "Niederreiter", dim, n_points);
    } catch (const std::runtime_error& e) {
        std::cout << "Niederreiter sequence not available: " << e.what() << std::endl;
    }
    
    // Test Lattice rules
    test_sequence_type<Real>(detail::qmc_sequence_type::lattice_rank1, 
                            "Lattice", dim, 1021); // Use prime number
}

BOOST_AUTO_TEST_CASE(test_scrambling) {
    using Real = double;
    const std::size_t dim = 3;
    const std::size_t n_points = 512;
    const std::size_t n_replicates = 10;
    
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Test with scrambling for randomized QMC
    qmc_options<Real> opts;
    opts.sequence_type = detail::qmc_sequence_type::sobol;
    opts.use_scrambling = true;
    opts.n_replicates = n_replicates;
    opts.scramble_seed = 12345;
    
    auto result = integrate_qmc_ex<Real>(product_function<Real>, box, n_points, opts);
    
    BOOST_CHECK_EQUAL(result.status, status_code::success);
    BOOST_CHECK_GT(result.error, Real(0)); // Should have error estimate with RQMC
    
    Real exact = std::pow(Real(0.5), Real(dim));
    std::cout << "\nRQMC with " << n_replicates << " replicates:" << std::endl;
    std::cout << "  Value: " << result.value << " ± " << result.error << std::endl;
    std::cout << "  Exact: " << exact << std::endl;
    
    // Check that result is within error bounds
    BOOST_CHECK_LE(std::abs(result.value - exact), Real(3) * result.error);
}

BOOST_AUTO_TEST_CASE(test_high_dimension) {
    using Real = double;
    const std::size_t dim = 10;
    const std::size_t n_points = 8192;
    
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Test different sequences in high dimension
    std::vector<detail::qmc_sequence_type> sequences = {
        detail::qmc_sequence_type::sobol,
        detail::qmc_sequence_type::halton
    };
    
    std::cout << "\nHigh-dimensional integration (" << dim << "D):" << std::endl;
    
    for (auto seq_type : sequences) {
        qmc_options<Real> opts;
        opts.sequence_type = seq_type;
        opts.use_scrambling = false;
        
        if (seq_type == detail::qmc_sequence_type::lattice_rank1) {
            opts.lattice_points = 8191; // Prime number
        }
        
        auto result = integrate_qmc_ex<Real>(constant_function<Real>, box, n_points, opts);
        
        BOOST_CHECK_EQUAL(result.status, status_code::success);
        BOOST_CHECK_CLOSE(result.value, Real(1), Real(1e-6));
        
        std::string seq_name = (seq_type == detail::qmc_sequence_type::sobol) ? "Sobol" :
                               (seq_type == detail::qmc_sequence_type::halton) ? "Halton" :
                               "Unknown";
        std::cout << "  " << seq_name << ": " << result.value 
                  << " (exact = 1.0)" << std::endl;
    }
}

BOOST_AUTO_TEST_CASE(test_convergence_rates) {
    using Real = double;
    const std::size_t dim = 4;
    
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Test convergence with increasing number of points
    std::vector<std::size_t> point_counts = {64, 256, 1024, 4096, 16384};
    
    std::cout << "\nConvergence study for product function in " << dim << "D:" << std::endl;
    std::cout << std::setw(10) << "N" << std::setw(15) << "Sobol Error" 
              << std::setw(15) << "Halton Error" << std::setw(15) << "MC Error" << std::endl;
    
    Real exact = std::pow(Real(0.5), Real(dim));
    
    for (std::size_t n : point_counts) {
        // Sobol sequence
        qmc_options<Real> sobol_opts;
        sobol_opts.sequence_type = detail::qmc_sequence_type::sobol;
        auto sobol_result = integrate_qmc_ex<Real>(product_function<Real>, box, n, sobol_opts);
        Real sobol_error = std::abs(sobol_result.value - exact);
        
        // Halton sequence
        qmc_options<Real> halton_opts;
        halton_opts.sequence_type = detail::qmc_sequence_type::halton;
        auto halton_result = integrate_qmc_ex<Real>(product_function<Real>, box, n, halton_opts);
        Real halton_error = std::abs(halton_result.value - exact);
        
        // Monte Carlo for comparison (using scrambled Sobol with 1 replicate as proxy)
        qmc_options<Real> mc_opts;
        mc_opts.sequence_type = detail::qmc_sequence_type::sobol;
        mc_opts.use_scrambling = true;
        mc_opts.n_replicates = 10;
        auto mc_result = integrate_qmc_ex<Real>(product_function<Real>, box, n, mc_opts);
        Real mc_error = std::abs(mc_result.value - exact);
        
        std::cout << std::setw(10) << n 
                  << std::scientific << std::setprecision(3)
                  << std::setw(15) << sobol_error
                  << std::setw(15) << halton_error
                  << std::setw(15) << mc_error << std::endl;
        
        // Check that QMC sequences converge
        if (n >= 1024) {
            BOOST_CHECK_LT(sobol_error, Real(0.01));
            BOOST_CHECK_LT(halton_error, Real(0.1));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_tent_transform) {
    using Real = double;
    const std::size_t dim = 2;
    const std::size_t n_points = 1000;
    
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Function with boundary singularity
    auto boundary_func = [](const Real* x, std::size_t) -> Real {
        return Real(1) / std::sqrt(x[0] * x[1] + Real(1e-10));
    };
    
    // Test with and without tent transform
    qmc_options<Real> opts_with_tent;
    opts_with_tent.sequence_type = detail::qmc_sequence_type::sobol;
    opts_with_tent.use_tent_transform = true;
    
    qmc_options<Real> opts_no_tent;
    opts_no_tent.sequence_type = detail::qmc_sequence_type::sobol;
    opts_no_tent.use_tent_transform = false;
    
    auto result_with = integrate_qmc_ex<Real>(boundary_func, box, n_points, opts_with_tent);
    auto result_without = integrate_qmc_ex<Real>(boundary_func, box, n_points, opts_no_tent);
    
    std::cout << "\nTent transform effect on boundary singularity:" << std::endl;
    std::cout << "  With tent: " << result_with.value << std::endl;
    std::cout << "  Without tent: " << result_without.value << std::endl;
    
    // Both should succeed
    BOOST_CHECK_EQUAL(result_with.status, status_code::success);
    BOOST_CHECK_EQUAL(result_without.status, status_code::success);
}

BOOST_AUTO_TEST_CASE(test_vector_integration) {
    using Real = double;
    const std::size_t dim = 3;
    const std::size_t n_points = 1024;
    const std::size_t n_components = 3;
    
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Vector-valued function
    auto vector_func = [](const Real* x, Real* f, std::size_t ncomp) {
        f[0] = Real(1);  // Constant
        f[1] = x[0] + x[1] + x[2];  // Linear
        if (ncomp > 2) {
            f[2] = x[0] * x[1] * x[2];  // Product
        }
    };
    
    // Test plain QMC
    auto results = integrate_qmc_vector<Real>(vector_func, box, n_components, n_points);
    
    BOOST_CHECK_EQUAL(results.size(), n_components);
    BOOST_CHECK_CLOSE(results[0].value, Real(1), Real(1e-6));
    BOOST_CHECK_CLOSE(results[1].value, Real(1.5), Real(0.1));
    BOOST_CHECK_CLOSE(results[2].value, Real(0.125), Real(1));
    
    // Test RQMC
    auto rqmc_results = integrate_rqmc_vector<Real>(vector_func, box, n_components, 
                                                    n_points, 10);
    
    BOOST_CHECK_EQUAL(rqmc_results.size(), n_components);
    for (std::size_t i = 0; i < n_components; ++i) {
        BOOST_CHECK_GT(rqmc_results[i].error, Real(0));
    }
}

BOOST_AUTO_TEST_CASE(test_edge_cases) {
    using Real = double;
    
    // Test 1D case
    {
        hypercube<Real> box(1);
        box.lower[0] = Real(0);
        box.upper[0] = Real(1);
        
        auto f = [](const Real* x, std::size_t) { return x[0] * x[0]; };
        
        qmc_options<Real> opts;
        opts.sequence_type = detail::qmc_sequence_type::sobol;
        
        auto result = integrate_qmc_ex<Real>(f, box, 100, opts);
        BOOST_CHECK_EQUAL(result.status, status_code::success);
        BOOST_CHECK_CLOSE(result.value, Real(1.0/3.0), Real(1));
    }
    
    // Test with different bounds
    {
        hypercube<Real> box(2);
        box.lower[0] = Real(-1);
        box.upper[0] = Real(2);
        box.lower[1] = Real(0);
        box.upper[1] = Real(1);
        
        auto f = [](const Real* x, std::size_t) { return Real(1); };
        
        qmc_options<Real> opts;
        opts.sequence_type = detail::qmc_sequence_type::halton;
        
        auto result = integrate_qmc_ex<Real>(f, box, 500, opts);
        BOOST_CHECK_EQUAL(result.status, status_code::success);
        BOOST_CHECK_CLOSE(result.value, Real(3), Real(1e-3)); // Volume = 3*1 = 3
    }
}

BOOST_AUTO_TEST_CASE(test_halton_leap) {
    using Real = double;
    const std::size_t dim = 5;
    const std::size_t n_points = 1000;
    
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    // Test leaped Halton sequence
    qmc_options<Real> opts;
    opts.sequence_type = detail::qmc_sequence_type::halton;
    opts.halton_leap = 409;  // Common leap value for better uniformity
    
    auto result = integrate_qmc_ex<Real>(linear_function<Real>, box, n_points, opts);
    
    BOOST_CHECK_EQUAL(result.status, status_code::success);
    Real exact = Real(dim) / Real(2);
    BOOST_CHECK_CLOSE(result.value, exact, Real(5)); // 5% tolerance
    
    std::cout << "\nLeaped Halton (leap=" << opts.halton_leap << ") in " << dim 
              << "D: " << result.value << " (exact = " << exact << ")" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()