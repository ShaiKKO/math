// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_ORBIT_EVALUATOR_HPP
#define BOOST_MATH_CUBATURE_DETAIL_ORBIT_EVALUATOR_HPP

#include <boost/math/cubature/detail/genz_malik_9_7_tables.hpp>
#include <boost/math/cubature/detail/adaptivity.hpp>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Orbit-based function evaluations for Genz-Malik rules
/// \details Evaluates integrands at symmetric orbit groups and computes
///          directional fourth differences for optimal axis selection
template <typename Real>
class orbit_evaluator {
public:
  
  /// \brief Cached evaluations organized by orbit structure
  template <std::size_t Dim>
  struct orbit_values {
    Real f_center;                                    // Center point
    std::array<Real, Dim> f_axis_l1_plus;            // +λ₁ along each axis
    std::array<Real, Dim> f_axis_l1_minus;           // -λ₁ along each axis
    std::array<Real, Dim> f_axis_l2_plus;            // +λ₂ along each axis
    std::array<Real, Dim> f_axis_l2_minus;           // -λ₂ along each axis
    std::array<Real, Dim> f_axis_l3_plus;            // +λ₃ along each axis (deg-9 only)
    std::array<Real, Dim> f_axis_l3_minus;           // -λ₃ along each axis (deg-9 only)
    std::vector<Real> f_pairs;                       // Pair orbits (l1,l1), (l1,l2)
    std::vector<Real> f_triples;                     // Triple orbit (l1,l1,l1) for D≥3
    std::vector<Real> f_diag;                        // Full diagonal (corners)
    
    std::size_t total_evaluations;
    
    orbit_values() : f_center(0), total_evaluations(0) {
      std::fill(f_axis_l1_plus.begin(), f_axis_l1_plus.end(), Real(0));
      std::fill(f_axis_l1_minus.begin(), f_axis_l1_minus.end(), Real(0));
      std::fill(f_axis_l2_plus.begin(), f_axis_l2_plus.end(), Real(0));
      std::fill(f_axis_l2_minus.begin(), f_axis_l2_minus.end(), Real(0));
      std::fill(f_axis_l3_plus.begin(), f_axis_l3_plus.end(), Real(0));
      std::fill(f_axis_l3_minus.begin(), f_axis_l3_minus.end(), Real(0));
    }
  };
  
  /// \brief Evaluate function at orbit representatives
  template <typename F, std::size_t Dim, int Degree>
  static void evaluate_orbits(
      const F& f,
      const region<Real>& reg,
      orbit_values<Dim>& values)
  {
    // Get radii for the rule degree
    const Real lambda1 = static_cast<Real>(rules::detail::gm_lambda1);
    const Real lambda2 = static_cast<Real>(rules::detail::gm_lambda2);
    const Real lambda3 = static_cast<Real>(rules::detail::gm_lambda3);
    
    // Compute center and half-widths
    std::array<Real, Dim> center, half_widths;
    for (std::size_t i = 0; i < Dim; ++i) {
      center[i] = (reg.a[i] + reg.b[i]) / Real(2);
      half_widths[i] = (reg.b[i] - reg.a[i]) / Real(2);
    }
    
    // Workspace for evaluation point
    std::array<Real, Dim> point;
    std::size_t evals = 0;
    
    // 1. Center point
    values.f_center = evaluate_at(f, center);
    evals++;
    
    // 2. Axis points λ₁
    for (std::size_t d = 0; d < Dim; ++d) {
      // +λ₁ along axis d
      point = center;
      point[d] += lambda1 * half_widths[d];
      values.f_axis_l1_plus[d] = evaluate_at(f, point);
      
      // -λ₁ along axis d
      point = center;
      point[d] -= lambda1 * half_widths[d];
      values.f_axis_l1_minus[d] = evaluate_at(f, point);
      
      evals += 2;
    }
    
    // 3. Axis points λ₂
    for (std::size_t d = 0; d < Dim; ++d) {
      // +λ₂ along axis d
      point = center;
      point[d] += lambda2 * half_widths[d];
      values.f_axis_l2_plus[d] = evaluate_at(f, point);
      
      // -λ₂ along axis d
      point = center;
      point[d] -= lambda2 * half_widths[d];
      values.f_axis_l2_minus[d] = evaluate_at(f, point);
      
      evals += 2;
    }
    
    // 4. Axis points λ₃ (degree-9 only)
    if constexpr (Degree >= 9) {
      for (std::size_t d = 0; d < Dim; ++d) {
        // +λ₃ along axis d
        point = center;
        point[d] += lambda3 * half_widths[d];
        values.f_axis_l3_plus[d] = evaluate_at(f, point);
        
        // -λ₃ along axis d
        point = center;
        point[d] -= lambda3 * half_widths[d];
        values.f_axis_l3_minus[d] = evaluate_at(f, point);
        
        evals += 2;
      }
    }
    
    // 5. Pair orbits (λ₁,λ₁)
    values.f_pairs.clear();
    for (std::size_t i = 0; i < Dim; ++i) {
      for (std::size_t j = i + 1; j < Dim; ++j) {
        // All 4 sign combinations for (±λ₁, ±λ₁, 0, ...)
        for (int si : {-1, 1}) {
          for (int sj : {-1, 1}) {
            point = center;
            point[i] += si * lambda1 * half_widths[i];
            point[j] += sj * lambda1 * half_widths[j];
            values.f_pairs.push_back(evaluate_at(f, point));
            evals++;
          }
        }
      }
    }
    
    // 6. Full diagonal (corners with all coordinates non-zero)
    const std::size_t num_corners = 1 << Dim;
    values.f_diag.clear();
    values.f_diag.reserve(num_corners);
    
    // Generate all 2^D corner points
    for (std::size_t mask = 0; mask < num_corners; ++mask) {
      point = center;
      for (std::size_t d = 0; d < Dim; ++d) {
        int sign = (mask & (1 << d)) ? 1 : -1;
        point[d] += sign * half_widths[d]; // Use full width for corners
      }
      values.f_diag.push_back(evaluate_at(f, point));
      evals++;
    }
    
    values.total_evaluations = evals;
  }
  
  /// \brief Compute directional fourth differences for axis selection
  /// \details Implements the formula: f₃ᵢ - 2f₁ - 7*(f₂ᵢ - 2f₁)
  ///          where 7 = (λ₃/λ₂)² as per van Dooren and de Ridder
  template <std::size_t Dim>
  static std::array<Real, Dim> compute_fourth_differences(
      const orbit_values<Dim>& values)
  {
    std::array<Real, Dim> diffs;
    
    // The scaling factor 7 comes from (λ₃/λ₂)² in the original paper
    // For our radii: (λ₃/λ₂)² = (9/10)/(9/70) = 70/10 = 7
    const Real scale_factor = Real(7);
    
    for (std::size_t d = 0; d < Dim; ++d) {
      // f₁ = f_center
      Real f1 = values.f_center;
      
      // f₂ᵢ = f(c + λ₂*eᵢ) + f(c - λ₂*eᵢ)
      Real f2i = values.f_axis_l2_plus[d] + values.f_axis_l2_minus[d];
      
      // f₃ᵢ = f(c + λ₃*eᵢ) + f(c - λ₃*eᵢ) for degree-9
      // For degree-7, use λ₁ as substitute
      Real f3i;
      if (values.f_axis_l3_plus[d] != Real(0) || values.f_axis_l3_minus[d] != Real(0)) {
        f3i = values.f_axis_l3_plus[d] + values.f_axis_l3_minus[d];
      } else {
        // Fallback for degree-7: use λ₁ points
        f3i = values.f_axis_l1_plus[d] + values.f_axis_l1_minus[d];
      }
      
      // Fourth divided difference: |f₃ᵢ - 2f₁ - scale*(f₂ᵢ - 2f₁)|
      Real diff = std::abs(f3i - Real(2)*f1 - scale_factor*(f2i - Real(2)*f1));
      diffs[d] = diff;
    }
    
    return diffs;
  }
  
  /// \brief Select optimal splitting dimension based on fourth differences
  template <std::size_t Dim>
  static std::size_t select_split_dimension(
      const std::array<Real, Dim>& fourth_diffs,
      const region<Real>& reg)
  {
    std::size_t best_dim = 0;
    Real max_weighted_diff = Real(0);
    
    for (std::size_t d = 0; d < Dim; ++d) {
      // Weight by axis width to prefer splitting longer axes
      Real width = reg.b[d] - reg.a[d];
      Real weighted_diff = fourth_diffs[d] * width;
      
      if (weighted_diff > max_weighted_diff) {
        max_weighted_diff = weighted_diff;
        best_dim = d;
      }
    }
    
    return best_dim;
  }
  
  
private:
  /// \brief Helper to evaluate function with proper signature handling
  template <typename F, std::size_t Dim>
  static Real evaluate_at(const F& f, const std::array<Real, Dim>& point) {
    if constexpr (std::is_invocable_v<F, const Real*, std::size_t>) {
      // f(const Real*, std::size_t)
      return f(point.data(), Dim);
    } else if constexpr (std::is_invocable_v<F, const Real*>) {
      // f(const Real*) - legacy support
      return f(point.data());
    } else {
      // f(std::vector<Real>)
      std::vector<Real> vec(point.begin(), point.end());
      return f(vec);
    }
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_ORBIT_EVALUATOR_HPP