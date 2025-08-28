// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_POLICIES_HPP
#define BOOST_MATH_CUBATURE_POLICIES_HPP

// Experimental toggle for this preview module
#ifndef BOOST_MATH_CUBATURE_EXPERIMENTAL
#define BOOST_MATH_CUBATURE_EXPERIMENTAL 1
#endif

#include <boost/math/policies/policy.hpp>
#include <cstdint>

namespace boost { namespace math { namespace cubature {

/// \file policies.hpp
/// \brief Policy aliasing and result/status types for cubature front-ends.

/// Status codes for integration routines
enum class status_code { 
  success,           // Integration converged to requested tolerance
  maxeval_reached,   // Maximum evaluations reached
  cancelled,         // User cancellation requested  
  maxregions_reached,// Maximum regions limit reached (adaptive)
  dimension_error,   // Rules not available for this dimension
  not_implemented,   // Feature not available (e.g., missing dependencies)
  invalid_input      // Invalid input parameters
};

/// \brief Error reliability metrics for integration results
/// \details Provides confidence indicators and convergence diagnostics
template <typename Real>
struct reliability_metrics {
  Real chi2_probability{};       // Chi-squared probability of error reliability (0-1)
  Real convergence_rate{};        // Estimated convergence rate (negative = faster)
  Real error_ratio{};            // Ratio of final to initial error estimate
  Real reliability_factor{};     // Overall reliability score (0-1, higher is better)
  std::size_t refinement_depth{}; // Maximum refinement depth reached
  std::size_t regions_processed{};// Total regions processed (adaptive)
  bool monotone_convergence{};    // Whether error decreased monotonically
  Real condition_estimate{};      // Estimate of problem conditioning (lower is better)
  
  /// \brief Check if error estimate is reliable
  bool is_reliable() const noexcept {
    return reliability_factor > Real(0.5) && chi2_probability > Real(0.05);
  }
  
  /// \brief Get convergence quality assessment
  const char* convergence_quality() const noexcept {
    if (convergence_rate < Real(-2)) return "excellent";
    if (convergence_rate < Real(-1)) return "good";
    if (convergence_rate < Real(-0.5)) return "fair";
    if (convergence_rate < Real(0)) return "slow";
    return "poor";
  }
};

/// Result type returned by all algorithms
/// - value: integral estimate
/// - error: absolute error bound or standard error (QMC)
/// - evaluations: number of integrand evaluations performed
/// - status: status_code indicating outcome
/// - reliability: error reliability metrics and convergence diagnostics
template <class Real>
struct result {
  Real value{};
  Real error{};
  std::uint64_t evaluations{0};
  status_code status{status_code::success};
  reliability_metrics<Real> reliability{};
  
  /// \brief Check if result meets requested tolerance reliably
  bool is_converged(Real abs_tol, Real rel_tol) const noexcept {
    return status == status_code::success && 
           error <= abs_tol + rel_tol * std::abs(value) &&
           reliability.is_reliable();
  }
};

/// Policy alias forwarded from Boost.Math policies
using default_policy = boost::math::policies::policy<>;

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_POLICIES_HPP

