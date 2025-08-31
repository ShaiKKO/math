// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_CONSTANTS_HPP
#define BOOST_MATH_CUBATURE_CONSTANTS_HPP

#include <cstddef>

namespace boost { namespace math { namespace cubature {

/// \file constants.hpp
/// \brief Central location for cubature library constants
/// \details Provides compile-time constants used throughout the cubature library
///          to avoid magic numbers and improve maintainability.

namespace constants {

// Dimension limits
constexpr std::size_t minimum_dimension = 1;
constexpr std::size_t maximum_dimension = 15;
constexpr std::size_t maximum_qmc_dimension = 40;
constexpr std::size_t maximum_sparse_grid_dimension = 20;

// Adaptive integration defaults
constexpr std::size_t default_max_regions = 10000;
constexpr std::size_t default_max_depth = 50;
constexpr std::size_t default_initial_regions = 100;
constexpr std::size_t default_max_evaluations = 1000000;

// Numerical safety factors
constexpr double machine_epsilon_safety_factor = 100.0;
constexpr double error_safeguard_factor = 1.5;
constexpr double roundoff_detection_threshold = 0.1;
constexpr double dcuhre_empirical_safety_factor = 50.0;  // From DCUHRE paper

// Sparse grid parameters
constexpr unsigned default_sparse_grid_level = 5;
constexpr unsigned maximum_sparse_grid_level = 12;
constexpr std::size_t sparse_grid_cache_size = 1000;

// QMC parameters
constexpr std::size_t default_qmc_points = 10000;
constexpr std::size_t minimum_qmc_points = 16;
constexpr std::size_t default_rqmc_replicates = 10;
constexpr std::size_t maximum_rqmc_replicates = 100;

// Memory pool configuration
constexpr std::size_t region_pool_initial_size = 1000;
constexpr std::size_t region_pool_growth_factor = 2;

// Convergence criteria
constexpr std::size_t minimum_iterations_for_convergence_check = 3;
constexpr std::size_t convergence_history_window = 5;
constexpr double convergence_rate_threshold = -0.5;

// Performance thresholds
constexpr std::size_t vectorization_threshold = 4;
constexpr std::size_t parallel_execution_threshold = 1000;
constexpr std::size_t cache_line_size = 64;

} // namespace constants

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_CONSTANTS_HPP