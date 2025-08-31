// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_RELIABILITY_METRICS_HPP
#define BOOST_MATH_CUBATURE_DETAIL_RELIABILITY_METRICS_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Convergence history tracker for reliability analysis
/// \details Implementation for tracking integration convergence behavior
template <typename Real>
class convergence_history {
public:
  std::vector<Real> error_sequence;      // Error estimates at each iteration
  std::vector<Real> integral_sequence;   // Integral estimates at each iteration
  std::vector<std::size_t> eval_sequence;// Cumulative evaluations at each step
  std::vector<std::size_t> region_sequence; // Number of regions at each step
  
  /// \brief Record a new iteration
  void record(Real integral, Real error, std::size_t evals, std::size_t regions) {
    integral_sequence.push_back(integral);
    error_sequence.push_back(error);
    eval_sequence.push_back(evals);
    region_sequence.push_back(regions);
  }
  
  /// \brief Compute convergence rate using linear regression on log(error) vs log(evaluations)
  Real compute_convergence_rate() const {
    if (error_sequence.size() < 3) {
      return Real(0);  // Not enough data
    }
    
    // Use last 5 points or all available
    std::size_t n = std::min(std::size_t(5), error_sequence.size());
    std::size_t start = error_sequence.size() - n;
    
    // Compute linear regression on log-log scale
    Real sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;
    std::size_t valid_points = 0;
    
    for (std::size_t i = start; i < error_sequence.size(); ++i) {
      if (error_sequence[i] > 0 && eval_sequence[i] > 0) {
        Real log_evals = std::log(static_cast<Real>(eval_sequence[i]));
        Real log_error = std::log(error_sequence[i]);
        
        sum_x += log_evals;
        sum_y += log_error;
        sum_xx += log_evals * log_evals;
        sum_xy += log_evals * log_error;
        valid_points++;
      }
    }
    
    if (valid_points < 2) {
      return Real(0);
    }
    
    // Slope of log(error) vs log(evaluations) gives convergence rate
    Real denom = valid_points * sum_xx - sum_x * sum_x;
    if (std::abs(denom) < std::numeric_limits<Real>::epsilon()) {
      return Real(0);
    }
    
    return (valid_points * sum_xy - sum_x * sum_y) / denom;
  }
  
  /// \brief Check if error decreases monotonically
  bool is_monotone() const {
    if (error_sequence.size() < 2) {
      return true;
    }
    
    for (std::size_t i = 1; i < error_sequence.size(); ++i) {
      if (error_sequence[i] > error_sequence[i-1] * Real(1.1)) {
        // Allow 10% tolerance for numerical noise
        return false;
      }
    }
    return true;
  }
  
  /// \brief Compute error ratio (final/initial)
  Real compute_error_ratio() const {
    if (error_sequence.empty()) {
      return Real(1);
    }
    if (error_sequence.size() == 1) {
      return Real(1);
    }
    
    Real initial = error_sequence[0];
    Real final = error_sequence.back();
    
    if (initial <= 0) {
      return Real(1);
    }
    
    return final / initial;
  }
  
  /// \brief Estimate condition number from integral variation
  Real estimate_condition() const {
    if (integral_sequence.size() < 2) {
      return Real(1);
    }
    
    // Find min and max integral estimates
    auto minmax = std::minmax_element(integral_sequence.begin(), integral_sequence.end());
    Real min_val = *minmax.first;
    Real max_val = *minmax.second;
    
    // Compute relative variation
    Real avg_val = std::accumulate(integral_sequence.begin(), integral_sequence.end(), Real(0)) 
                   / integral_sequence.size();
    
    if (std::abs(avg_val) < std::numeric_limits<Real>::epsilon()) {
      // Near zero integral - use absolute variation
      return (max_val - min_val) > Real(0.1) ? Real(10) : Real(1);
    }
    
    Real relative_var = (max_val - min_val) / std::abs(avg_val);
    
    // Map to condition estimate (1 = well-conditioned, >10 = ill-conditioned)
    return Real(1) + Real(9) * std::min(relative_var, Real(1));
  }
};

/// \brief Chi-squared reliability calculator
/// \details Implements CUBA-style chi-squared probability calculation
template <typename Real>
class chi2_reliability {
public:
  /// \brief Compute chi-squared probability from error distribution
  /// \param region_errors Vector of error estimates from different regions
  /// \param total_error Overall error estimate
  /// \return Chi-squared probability (0-1), higher indicates more reliable
  static Real compute(const std::vector<Real>& region_errors, Real total_error) {
    if (region_errors.size() < 3 || total_error <= 0) {
      return Real(0.5);  // Default neutral reliability
    }
    
    // Compute chi-squared statistic
    Real chi2 = 0;
    Real mean_error = total_error / region_errors.size();
    
    for (const auto& err : region_errors) {
      if (err > 0) {
        Real normalized = (err - mean_error) / mean_error;
        chi2 += normalized * normalized;
      }
    }
    
    // Degrees of freedom
    std::size_t dof = region_errors.size() - 1;
    
    // Convert to probability using incomplete gamma function
    // P(X > chi2) where X ~ Chi2(dof)
    try {
      Real prob = boost::math::gamma_q(Real(dof) / 2, chi2 / 2);
      return std::max(Real(0), std::min(Real(1), prob));
    } catch (...) {
      // Fallback for extreme values
      return chi2 < Real(dof) ? Real(0.8) : Real(0.2);
    }
  }
};

/// \brief Overall reliability factor calculator
/// \details Combines multiple reliability indicators into single score
template <typename Real>
class reliability_calculator {
public:
  /// \brief Compute overall reliability factor
  /// \param chi2_prob Chi-squared probability
  /// \param convergence_rate Convergence rate (negative is better)
  /// \param error_ratio Final/initial error ratio
  /// \param is_monotone Whether convergence was monotone
  /// \param condition_est Condition number estimate
  /// \return Overall reliability score (0-1)
  static Real compute_reliability_factor(
      Real chi2_prob,
      Real convergence_rate,
      Real error_ratio,
      bool is_monotone,
      Real condition_est)
  {
    // Individual component scores (0-1)
    
    // Chi-squared score: direct use of probability
    Real chi2_score = chi2_prob;
    
    // Convergence score: better for more negative rates
    Real conv_score;
    if (convergence_rate >= 0) {
      conv_score = Real(0.1);  // Poor convergence
    } else if (convergence_rate <= -2) {
      conv_score = Real(1.0);  // Excellent convergence
    } else {
      // Linear interpolation between -2 and 0
      conv_score = Real(0.1) + Real(0.9) * (-convergence_rate / 2);
    }
    
    // Error ratio score: better for smaller ratios
    Real ratio_score;
    if (error_ratio >= 1) {
      ratio_score = Real(0.1);  // No improvement
    } else if (error_ratio <= 0.001) {
      ratio_score = Real(1.0);  // Excellent improvement
    } else {
      // Logarithmic scale
      ratio_score = Real(1) - (std::log10(error_ratio) + 3) / 3;
      ratio_score = std::max(Real(0.1), std::min(Real(1), ratio_score));
    }
    
    // Monotone bonus
    Real monotone_score = is_monotone ? Real(1.0) : Real(0.5);
    
    // Condition score: inverse relationship
    Real cond_score;
    if (condition_est <= 1) {
      cond_score = Real(1.0);  // Well-conditioned
    } else if (condition_est >= 10) {
      cond_score = Real(0.2);  // Ill-conditioned
    } else {
      cond_score = Real(1) - Real(0.8) * (condition_est - 1) / 9;
    }
    
    // Weighted average
    const Real w_chi2 = Real(0.3);
    const Real w_conv = Real(0.25);
    const Real w_ratio = Real(0.2);
    const Real w_monotone = Real(0.15);
    const Real w_cond = Real(0.1);
    
    Real weighted_sum = w_chi2 * chi2_score +
                       w_conv * conv_score +
                       w_ratio * ratio_score +
                       w_monotone * monotone_score +
                       w_cond * cond_score;
    
    return std::max(Real(0), std::min(Real(1), weighted_sum));
  }
  
  /// \brief Compute reliability metrics from convergence history
  static reliability_metrics<Real> compute_metrics(
      const convergence_history<Real>& history,
      const std::vector<Real>& region_errors,
      Real final_error,
      std::size_t max_depth)
  {
    reliability_metrics<Real> metrics;
    
    // Basic statistics
    metrics.refinement_depth = max_depth;
    metrics.regions_processed = history.region_sequence.empty() ? 0 : history.region_sequence.back();
    
    // Convergence analysis
    metrics.convergence_rate = history.compute_convergence_rate();
    metrics.error_ratio = history.compute_error_ratio();
    metrics.monotone_convergence = history.is_monotone();
    metrics.condition_estimate = history.estimate_condition();
    
    // Chi-squared reliability
    metrics.chi2_probability = chi2_reliability<Real>::compute(region_errors, final_error);
    
    // Overall reliability
    metrics.reliability_factor = compute_reliability_factor(
        metrics.chi2_probability,
        metrics.convergence_rate,
        metrics.error_ratio,
        metrics.monotone_convergence,
        metrics.condition_estimate);
    
    return metrics;
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_RELIABILITY_METRICS_HPP