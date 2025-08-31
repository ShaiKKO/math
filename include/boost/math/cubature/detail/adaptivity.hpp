// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_ADAPTIVITY_HPP
#define BOOST_MATH_CUBATURE_DETAIL_ADAPTIVITY_HPP

// Include STL headers before opening namespace to avoid lookup issues
#include <vector>
#include <queue>
#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <type_traits>
#include <memory>

// Reopen namespace after STL includes to ensure proper name lookup
namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Region structure for adaptive cubature (DCUHRE-style)
/// \details Represents a hypercube subregion with bounds, estimates, and error tracking.
///  Following ALGORITHMS.md §1.2 and TECHNICAL_BLUEPRINT.md §2.1
template <typename Real>
struct region {
  std::vector<Real> a;           // Lower bounds (size d)
  std::vector<Real> b;           // Upper bounds (size d)
  Real estimate_coarse;          // Coarse rule estimate (degree 5 or 7)
  Real estimate_fine;            // Fine rule estimate (degree 7 or 9)
  Real error;                    // Error estimate |fine - coarse| (safeguarded)
  std::size_t split_dim;         // Dimension to split on next refinement
  std::size_t evaluations;       // Function evaluations used for this region
  std::shared_ptr<void> cached_values;  // Cached function evaluations for reuse
  
  // Constructor
  region(std::size_t dim = 0) 
    : a(dim), b(dim), estimate_coarse(0), estimate_fine(0), 
      error(0), split_dim(dim), evaluations(0), cached_values(nullptr) {}  // Initialize split_dim to dim (invalid)
  
  // Get dimension
  std::size_t dimension() const { return a.size(); }
  
  // Get volume
  Real volume() const {
    Real vol = Real(1);
    for (std::size_t i = 0; i < a.size(); ++i) {
      vol *= (b[i] - a[i]);
    }
    return vol;
  }
  
  // Get center point
  std::vector<Real> center() const {
    std::vector<Real> c(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      c[i] = (a[i] + b[i]) / Real(2);
    }
    return c;
  }
  
  // Get half-widths
  std::vector<Real> half_widths() const {
    std::vector<Real> hw(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      hw[i] = (b[i] - a[i]) / Real(2);
    }
    return hw;
  }
};

/// \brief Comparator for priority queue (max-heap by error)
template <typename Real>
struct region_error_comparator {
  bool operator()(const region<Real>& r1, const region<Real>& r2) const {
    return r1.error < r2.error; // Max-heap: largest error first
  }
};

/// \brief Priority queue type for regions
template <typename Real>
using region_priority_queue = std::priority_queue<
  region<Real>,
  std::vector<region<Real>>,
  region_error_comparator<Real>
>;

/// \brief Compute directional variation for axis selection
/// \details Uses symmetric node pairs to estimate variation along each axis
///  Following ALGORITHMS.md §1.2 - axis selection criterion
template <typename Real, typename NodeContainer, typename ValueContainer>
std::size_t compute_split_dimension(
    const std::vector<Real>& /*center*/,
    const std::vector<Real>& half_widths,
    const NodeContainer& nodes,
    const ValueContainer& values,
    std::size_t dim)
{
  std::vector<Real> variations(dim, Real(0));
  
  // Compute directional variations using symmetric node pairs
  // For Genz-Malik rules, we have specific orbit structures
  
  // Analyze function variation along each axis using node pairs
  if (!nodes.empty() && !values.empty()) {
    // For each dimension, look at symmetric pairs along that axis
    for (std::size_t d = 0; d < dim; ++d) {
      Real sum_diff_sq = Real(0);
      std::size_t pair_count = 0;
      
      // Find symmetric pairs along dimension d
      for (std::size_t i = 0; i < nodes.size(); ++i) {
        for (std::size_t j = i + 1; j < nodes.size(); ++j) {
          // Check if nodes differ only in dimension d
          bool is_symmetric_pair = true;
          for (std::size_t k = 0; k < dim; ++k) {
            if (k != d) {
              // Other dimensions should match (within tolerance)
              if (std::abs(nodes[i][k] - nodes[j][k]) > Real(1e-10)) {
                is_symmetric_pair = false;
                break;
              }
            }
          }
          
          if (is_symmetric_pair) {
            // This is a symmetric pair along dimension d
            Real diff = values[i] - values[j];
            Real dist = std::abs(nodes[i][d] - nodes[j][d]);
            if (dist > Real(1e-10)) {
              // Normalized difference quotient
              Real grad_est = std::abs(diff) / dist;
              sum_diff_sq += grad_est * grad_est;
              ++pair_count;
            }
          }
        }
      }
      
      // Average variation estimate for dimension d
      if (pair_count > 0) {
        variations[d] = std::sqrt(sum_diff_sq / Real(pair_count));
      }
    }
  }
  
  // Find axis with maximum weighted variation
  std::size_t max_dim = 0;
  Real max_var = variations[0] * half_widths[0];
  
  for (std::size_t i = 1; i < dim; ++i) {
    // Weight by axis length^5 (as per Genz-Malik paper)
    Real width = half_widths[i];
    Real weighted_var = variations[i] * width * width * width * width * width;
    if (weighted_var > max_var) {
      max_var = weighted_var;
      max_dim = i;
    }
  }
  
  return max_dim;
}

/// \brief Bisect a region along specified dimension
/// \details Creates two child regions by splitting at the midpoint
template <typename Real>
std::pair<region<Real>, region<Real>> bisect_region(
    const region<Real>& parent,
    std::size_t split_dim)
{
  region<Real> left(parent.dimension());
  region<Real> right(parent.dimension());
  
  // Copy bounds
  left.a = parent.a;
  left.b = parent.b;
  right.a = parent.a;
  right.b = parent.b;
  
  // Split at midpoint
  Real midpoint = (parent.a[split_dim] + parent.b[split_dim]) / Real(2);
  left.b[split_dim] = midpoint;
  right.a[split_dim] = midpoint;
  
  return {left, right};
}

/// \brief Error safeguard using spread estimate
/// \details Guards against spuriously small embedded differences
///  Following ALGORITHMS.md §1.3 - safeguards
template <typename Real>
Real safeguard_error(Real embedded_diff, Real spread_estimate, Real /*tolerance*/) {
  // If embedded difference is too small relative to tolerance,
  // use spread estimate instead (max - min of partial sums)
  const Real machine_eps = std::numeric_limits<Real>::epsilon();
  
  if (std::abs(embedded_diff) < machine_eps * std::abs(spread_estimate)) {
    return spread_estimate;
  }
  
  return std::abs(embedded_diff);
}

/// \brief Kahan summation algorithm for accurate accumulation
/// \details Compensated summation following Boost.Math patterns for numerical stability
/// \note This implementation follows the algorithm from Boost.Accumulators
template <typename Real>
class kahan_accumulator {
private:
  Real sum_;
  Real compensation_;
  
public:
  using value_type = Real;
  
  kahan_accumulator() noexcept : sum_(Real(0)), compensation_(Real(0)) {}
  
  explicit kahan_accumulator(Real initial) noexcept 
    : sum_(initial), compensation_(Real(0)) {}
  
  void add(Real value) noexcept {
    // Kahan's compensated summation algorithm
    Real y = value - compensation_;
    Real t = sum_ + y;
    compensation_ = (t - sum_) - y;
    sum_ = t;
  }
  
  kahan_accumulator& operator+=(Real value) noexcept {
    add(value);
    return *this;
  }
  
  Real sum() const noexcept { return sum_; }
  Real value() const noexcept { return sum_; }
  Real compensation() const noexcept { return compensation_; }
  
  void reset() noexcept {
    sum_ = Real(0);
    compensation_ = Real(0);
  }
  
  void reset(Real initial) noexcept {
    sum_ = initial;
    compensation_ = Real(0);
  }
};

/// \brief Neumaier summation (improved Kahan) for even better accuracy
/// \details Handles cases where compensation becomes larger than sum
template <typename Real>
class neumaier_accumulator {
private:
  Real sum_;
  Real compensation_;
  
public:
  using value_type = Real;
  
  neumaier_accumulator() noexcept : sum_(Real(0)), compensation_(Real(0)) {}
  
  explicit neumaier_accumulator(Real initial) noexcept
    : sum_(initial), compensation_(Real(0)) {}
  
  void add(Real value) noexcept {
    Real t = sum_ + value;
    if (std::abs(sum_) >= std::abs(value)) {
      compensation_ += (sum_ - t) + value;
    } else {
      compensation_ += (value - t) + sum_;
    }
    sum_ = t;
  }
  
  neumaier_accumulator& operator+=(Real value) noexcept {
    add(value);
    return *this;
  }
  
  Real sum() const noexcept { return sum_ + compensation_; }
  Real value() const noexcept { return sum_ + compensation_; }
  
  void reset() noexcept {
    sum_ = Real(0);
    compensation_ = Real(0);
  }
};

// Alias for backward compatibility
template <typename Real>
using two_sum_accumulator = kahan_accumulator<Real>;

/// \brief Check if error is within requested tolerance
template <typename Real>
bool is_converged_within_tolerance(Real error, Real integral, 
                      Real abs_tol, Real rel_tol) {
  return error <= std::max(abs_tol, rel_tol * std::abs(integral));
}

/// \brief Transform point from unit hypercube to region
template <typename Real>
void transform_point(const std::vector<Real>& unit_point,
                    const region<Real>& reg,
                    std::vector<Real>& transformed) {
  const std::size_t dim = unit_point.size();
  transformed.resize(dim);
  
  for (std::size_t i = 0; i < dim; ++i) {
    transformed[i] = reg.a[i] + unit_point[i] * (reg.b[i] - reg.a[i]);
  }
}

/// \brief Compute Jacobian for region transformation
template <typename Real>
Real compute_jacobian(const region<Real>& reg) {
  return reg.volume();
}

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_ADAPTIVITY_HPP