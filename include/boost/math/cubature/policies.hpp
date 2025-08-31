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
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/cubature/constants.hpp>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace boost { namespace math { namespace cubature {

/// \file policies.hpp
/// \brief Policy support and result types for cubature integration
/// \details Provides comprehensive policy integration with Boost.Math,
///          including error handling, precision control, and evaluation limits.
/// 
/// \section policies_usage Policy Usage
/// \code
/// // Use default policies
/// auto result = integrate_adaptive(f, box, 1e-6, 1e-6);
/// 
/// // Custom policy with specific error handling
/// using my_policy = boost::math::policies::policy<
///     boost::math::policies::domain_error<boost::math::policies::throw_on_error>,
///     boost::math::policies::evaluation_error<boost::math::policies::errno_on_error>
/// >;
/// auto result = integrate_adaptive(f, box, 1e-6, 1e-6, 0, my_policy{});
/// 
/// // Policy for high precision
/// using high_prec_policy = boost::math::policies::policy<
///     boost::math::policies::precision<50>
/// >;
/// \endcode

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
  
  /// \brief Add iteration data for tracking convergence
  void add_iteration(Real error, std::size_t /*evaluations*/) {
    // Simple tracking - can be enhanced later
    error_ratio = error;
    regions_processed++;
  }
};

/// \brief Result type returned by all integration algorithms
/// 
/// \details Encapsulates the integration result along with error estimates,
///          performance metrics, and reliability indicators. This structure
///          provides comprehensive diagnostics to assess integration quality.
///
/// \tparam Real Floating-point type used for computation
///
/// Members:
/// - value: Estimated integral value
/// - error: Conservative absolute error bound (adaptive) or standard error (QMC)
/// - evaluations: Total number of integrand evaluations performed  
/// - status: Integration outcome (success, maxeval_reached, etc.)
/// - reliability: Advanced metrics for assessing error estimate reliability
///
/// \invariant error >= 0 (non-negative error estimate)
/// \invariant evaluations > 0 if status == success
/// 
/// \see reliability_metrics for detailed convergence diagnostics
/// \see status_code for possible integration outcomes
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

/// \brief Policy alias forwarded from Boost.Math policies
using default_policy = boost::math::policies::policy<>;

/// \brief Helper to extract precision from policy
template <typename Real, typename Policy = default_policy>
struct precision_traits {
  using precision_type = typename boost::math::policies::precision<Real, Policy>::type;
  static constexpr int digits = precision_type::value;
  static constexpr bool is_high_precision = digits > 53;  // More than double
};

/// \brief Policy-aware numerical accumulator with Kahan summation
/// \tparam Real Floating point type
/// \tparam Policy Boost.Math policy type
template <typename Real, typename Policy = default_policy>
class policy_accumulator {
private:
  Real sum_;
  Real correction_;  // Kahan correction term
  const Policy& policy_;
  
public:
  explicit policy_accumulator(const Policy& pol = Policy{})
    : sum_(Real(0)), correction_(Real(0)), policy_(pol) {}
  
  /// \brief Add value with Kahan summation for better precision
  void add(Real value) {
    // Kahan summation algorithm
    Real y = value - correction_;
    Real t = sum_ + y;
    correction_ = (t - sum_) - y;
    sum_ = t;
    
    // Check for overflow/underflow according to policy
    if (!std::isfinite(sum_)) {
      boost::math::policies::raise_overflow_error<Real>(
        "cubature::policy_accumulator", 
        "Accumulator overflow", policy_);
    }
  }
  
  /// \brief Get accumulated sum
  Real sum() const { return sum_; }
  
  /// \brief Reset accumulator
  void reset() { 
    sum_ = Real(0); 
    correction_ = Real(0); 
  }
};

/// \brief Policy-controlled evaluation limits
template <typename Policy>
struct evaluation_limits {
  static std::size_t max_evaluations(const Policy& /*pol*/) {
    // Could be customized by policy in future
    return std::size_t(1000000);
  }
  
  static std::size_t max_regions(const Policy& /*pol*/) {
    return std::size_t(100000);
  }
  
  static std::size_t max_depth(const Policy& /*pol*/) {
    return std::size_t(50);
  }
};

/// \brief Helper to handle errors according to policy
template <typename Real, typename Policy>
inline void raise_integration_error(
    const char* function,
    const char* message, 
    Real error_estimate,
    const Policy& pol)
{
  // Use Boost.Math error handling based on policy
  boost::math::policies::raise_evaluation_error<Real>(
    function, message, error_estimate, pol);
}

/// \brief Check and handle dimension errors
template <typename Policy>
inline bool check_dimension(
    std::size_t dim,
    std::size_t max_dim,
    const char* function,
    const Policy& pol)
{
  if (dim > max_dim) {
    boost::math::policies::raise_domain_error<std::size_t>(
      function,
      "Dimension %1% exceeds maximum supported dimension",
      dim, pol);
    return false;
  }
  return true;
}

/// \brief Policy-aware tolerance adjustment for multiprecision
template <typename Real, typename Policy>
struct tolerance_adjuster {
  template <typename P = Policy>
  static typename std::enable_if<precision_traits<P>::is_high_precision, Real>::type
  adjust_absolute(Real tol, const Policy& /*pol*/) {
    // For high precision types, ensure tolerance is meaningful
    return std::max(tol, std::numeric_limits<Real>::epsilon() * Real(constants::machine_epsilon_safety_factor));
  }
  
  template <typename P = Policy>
  static typename std::enable_if<!precision_traits<P>::is_high_precision, Real>::type
  adjust_absolute(Real tol, const Policy& /*pol*/) {
    return std::max(tol, std::numeric_limits<Real>::epsilon());
  }
  
  static Real adjust_relative(Real tol, const Policy& /*pol*/) {
    // Relative tolerance should be at least sqrt(epsilon)
    Real min_rel = std::sqrt(std::numeric_limits<Real>::epsilon());
    return std::max(tol, min_rel);
  }
};

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_POLICIES_HPP

