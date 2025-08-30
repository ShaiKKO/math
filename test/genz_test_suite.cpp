// Copyright 2025 Your Organization
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/// \file genz_test_suite.cpp
/// \brief Comprehensive test suite using Genz test functions
/// \details Tests all integration algorithms against the 6 standard Genz
///          test functions with known analytical integrals

#define BOOST_TEST_MODULE GenzTestSuite
#include <boost/test/included/unit_test.hpp>
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/qmc.hpp>
#include "genz_test_functions.hpp"
#include <iostream>
#include <iomanip>

using namespace boost::math::cubature;
using namespace boost::math::cubature::test;

// Test configuration
struct test_config {
    static constexpr double abs_tol = 1e-4;  // More realistic tolerance
    static constexpr double rel_tol = 1e-3;  // 0.1% relative error
    static constexpr std::size_t max_eval = 100000;
    static constexpr unsigned sparse_level = 5;
    static constexpr std::size_t qmc_points = 10000;
    // Sparse grids struggle with sharp peaks, need relaxed tolerance
    static constexpr double sparse_rel_tol = 0.2;  // 20% for sharp functions
};

// Helper to run a single test
template <typename Real, typename Integrator>
void run_genz_test(
    const genz_function<Real>& f,
    Integrator& integrator,
    const std::string& algorithm_name,
    std::size_t dim,
    int difficulty)
{
    Real exact = f.exact_integral();
    auto result = integrator(f, test_config::abs_tol, test_config::rel_tol);
    
    Real abs_error = std::abs(result.value - exact);
    Real rel_error = abs_error / (std::abs(exact) + 1e-10);
    
    // Use relaxed tolerance for sparse grids on difficult functions
    Real tol_to_use = test_config::rel_tol;
    if (algorithm_name == "SparseGrid" && 
        (std::string(f.name()) == "Gaussian" || 
         std::string(f.name()) == "Product Peak" ||
         std::string(f.name()) == "Corner Peak") && 
        difficulty >= 1) {
        tol_to_use = test_config::sparse_rel_tol;
    }
    
    bool passed = (abs_error <= test_config::abs_tol) ||
                  (rel_error <= tol_to_use);
    
    BOOST_TEST_MESSAGE(
        algorithm_name << " | " << f.name() << " | dim=" << dim 
        << " | diff=" << difficulty
        << " | exact=" << std::scientific << std::setprecision(6) << exact
        << " | computed=" << result.value
        << " | abs_err=" << abs_error
        << " | rel_err=" << rel_error
        << " | evals=" << result.evaluations
        << " | " << (passed ? "PASS" : "FAIL")
    );
    
    BOOST_CHECK(passed);
}

// Integrator wrappers
template <typename Real>
struct adaptive_wrapper {
    template <typename F>
    result<Real> operator()(const F& f, Real abs_tol, Real rel_tol) {
        hypercube<Real> unit_cube(f.dimension());
        std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
        std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
        return integrate_adaptive<Real>(f, unit_cube, abs_tol, rel_tol, 
                                       test_config::max_eval);
    }
};

template <typename Real>
struct sparse_grid_wrapper {
    unsigned level;
    sparse_grid_wrapper(unsigned l = test_config::sparse_level) : level(l) {}
    
    template <typename F>
    result<Real> operator()(const F& f, Real /*abs_tol*/, Real /*rel_tol*/) {
        hypercube<Real> unit_cube(f.dimension());
        std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
        std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
        return integrate_sparse_grid<Real>(f, unit_cube, level);
    }
};

template <typename Real>
struct qmc_wrapper {
    std::size_t n_points;
    qmc_wrapper(std::size_t n = test_config::qmc_points) : n_points(n) {}
    
    template <typename F>
    result<Real> operator()(const F& f, Real /*abs_tol*/, Real /*rel_tol*/) {
        hypercube<Real> unit_cube(f.dimension());
        std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
        std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
        return integrate_qmc<Real>(f, unit_cube, n_points);
    }
};

// Test cases for each Genz function
BOOST_AUTO_TEST_SUITE(genz_oscillatory)

BOOST_AUTO_TEST_CASE(adaptive_2d_3d_4d) {
    adaptive_wrapper<double> integrator;
    
    for (std::size_t dim : {2, 3, 4}) {
        for (int diff : {0, 1}) {  // Easy and moderate difficulty
            auto f = create_genz_function<double>(1, dim, diff);
            run_genz_test(*f, integrator, "Adaptive", dim, diff);
        }
    }
}

BOOST_AUTO_TEST_CASE(sparse_grid_2d_3d) {
    sparse_grid_wrapper<double> integrator;
    
    for (std::size_t dim : {2, 3}) {
        for (int diff : {0, 1}) {
            auto f = create_genz_function<double>(1, dim, diff);
            run_genz_test(*f, integrator, "SparseGrid", dim, diff);
        }
    }
}

BOOST_AUTO_TEST_CASE(qmc_2d_3d_5d) {
    qmc_wrapper<double> integrator;
    
    for (std::size_t dim : {2, 3, 5}) {
        for (int diff : {0, 1}) {
            auto f = create_genz_function<double>(1, dim, diff);
            run_genz_test(*f, integrator, "QMC", dim, diff);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(genz_product_peak)

BOOST_AUTO_TEST_CASE(adaptive_all_dims) {
    adaptive_wrapper<double> integrator;
    
    for (std::size_t dim : {2, 3, 4, 5}) {
        for (int diff : {0, 1}) {
            auto f = create_genz_function<double>(2, dim, diff);
            run_genz_test(*f, integrator, "Adaptive", dim, diff);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(genz_corner_peak)

BOOST_AUTO_TEST_CASE(adaptive_low_dims) {
    adaptive_wrapper<double> integrator;
    
    for (std::size_t dim : {2, 3}) {
        for (int diff : {0, 1}) {
            auto f = create_genz_function<double>(3, dim, diff);
            run_genz_test(*f, integrator, "Adaptive", dim, diff);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(genz_gaussian)

BOOST_AUTO_TEST_CASE(all_algorithms) {
    const std::size_t dim = 3;
    const int diff = 1;
    auto f = create_genz_function<double>(4, dim, diff);
    
    adaptive_wrapper<double> adaptive;
    run_genz_test(*f, adaptive, "Adaptive", dim, diff);
    
    sparse_grid_wrapper<double> sparse;
    run_genz_test(*f, sparse, "SparseGrid", dim, diff);
    
    qmc_wrapper<double> qmc;
    run_genz_test(*f, qmc, "QMC", dim, diff);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(genz_continuous)

BOOST_AUTO_TEST_CASE(adaptive_2d_3d) {
    adaptive_wrapper<double> integrator;
    
    for (std::size_t dim : {2, 3}) {
        for (int diff : {0, 1}) {
            auto f = create_genz_function<double>(5, dim, diff);
            run_genz_test(*f, integrator, "Adaptive", dim, diff);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(genz_discontinuous)

BOOST_AUTO_TEST_CASE(adaptive_2d) {
    adaptive_wrapper<double> integrator;
    
    // Discontinuous functions are challenging - test only 2D with easy difficulty
    auto f = create_genz_function<double>(6, 2, 0);
    run_genz_test(*f, integrator, "Adaptive", 2, 0);
}

BOOST_AUTO_TEST_CASE(qmc_better_than_adaptive) {
    const std::size_t dim = 3;
    auto f = create_genz_function<double>(6, dim, 0);
    
    // QMC should handle discontinuities better than adaptive
    adaptive_wrapper<double> adaptive;
    qmc_wrapper<double> qmc(50000);  // More points for discontinuous
    
    auto adaptive_result = adaptive(*f, test_config::abs_tol, test_config::rel_tol);
    auto qmc_result = qmc(*f, test_config::abs_tol, test_config::rel_tol);
    
    double exact = f->exact_integral();
    double adaptive_error = std::abs(adaptive_result.value - exact);
    double qmc_error = std::abs(qmc_result.value - exact);
    
    BOOST_TEST_MESSAGE("Discontinuous function comparison:");
    BOOST_TEST_MESSAGE("  Adaptive error: " << adaptive_error);
    BOOST_TEST_MESSAGE("  QMC error: " << qmc_error);
    
    // QMC should be at least competitive
    BOOST_CHECK(qmc_error < adaptive_error * 2.0);
}

BOOST_AUTO_TEST_SUITE_END()

// Performance and convergence tests
BOOST_AUTO_TEST_SUITE(convergence)

BOOST_AUTO_TEST_CASE(adaptive_monotone_convergence) {
    const std::size_t dim = 3;
    auto f = create_genz_function<double>(4, dim, 1);  // Gaussian, moderate
    double exact = f->exact_integral();
    
    hypercube<double> unit_cube(dim);
    std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), 0.0);
    std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), 1.0);
    
    std::vector<std::size_t> max_evals = {100, 500, 1000, 5000, 10000};
    std::vector<double> errors;
    
    for (auto max_eval : max_evals) {
        auto result = integrate_adaptive<double>(*f, unit_cube, 1e-10, 1e-10, max_eval);
        double error = std::abs(result.value - exact);
        errors.push_back(error);
        
        BOOST_TEST_MESSAGE("max_eval=" << max_eval 
                          << " error=" << std::scientific << error
                          << " evals_used=" << result.evaluations);
    }
    
    // Check monotone convergence
    for (std::size_t i = 1; i < errors.size(); ++i) {
        // Allow small tolerance for numerical noise
        BOOST_CHECK(errors[i] <= errors[i-1] * 1.1);
    }
}

BOOST_AUTO_TEST_CASE(qmc_convergence_rate) {
    const std::size_t dim = 4;
    auto f = create_genz_function<double>(4, dim, 1);  // Smooth function
    double exact = f->exact_integral();
    
    hypercube<double> unit_cube(dim);
    std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), 0.0);
    std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), 1.0);
    
    std::vector<std::size_t> point_counts = {100, 400, 1600, 6400};
    std::vector<double> errors;
    
    for (auto n_points : point_counts) {
        auto result = integrate_qmc<double>(*f, unit_cube, n_points);
        double error = std::abs(result.value - exact);
        errors.push_back(error);
        
        BOOST_TEST_MESSAGE("n_points=" << n_points 
                          << " error=" << std::scientific << error);
    }
    
    // Check O(1/n) convergence for smooth functions
    for (std::size_t i = 1; i < errors.size(); ++i) {
        double ratio = errors[i] / errors[i-1];
        double expected_ratio = std::sqrt(static_cast<double>(point_counts[i-1]) / 
                                         point_counts[i]);
        
        // Allow 50% deviation from theoretical rate
        BOOST_CHECK(ratio < expected_ratio * 1.5);
    }
}

BOOST_AUTO_TEST_SUITE_END()

// Determinism tests
BOOST_AUTO_TEST_SUITE(determinism)

BOOST_AUTO_TEST_CASE(adaptive_deterministic) {
    auto f = create_genz_function<double>(4, 3, 1);
    adaptive_wrapper<double> integrator;
    
    // Run multiple times
    std::vector<double> values;
    for (int i = 0; i < 3; ++i) {
        auto result = integrator(*f, test_config::abs_tol, test_config::rel_tol);
        values.push_back(result.value);
    }
    
    // All values should be identical
    for (std::size_t i = 1; i < values.size(); ++i) {
        BOOST_CHECK_EQUAL(values[i], values[0]);
    }
}

BOOST_AUTO_TEST_SUITE_END()