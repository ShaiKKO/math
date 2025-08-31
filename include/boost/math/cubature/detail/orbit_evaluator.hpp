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

/// \file orbit_evaluator.hpp
/// \brief Orbit-based function evaluations for Genz-Malik rules
/// \details Implements efficient evaluation of integrands at symmetric orbit
///          groups and computes directional fourth differences for optimal
///          axis selection in adaptive subdivision.
/// 
/// \section orbit_theory Theory
/// Genz-Malik rules exploit the symmetry of the integration region:
/// - Center orbit: 1 point at origin
/// - Axis orbits: 2d points at ±λᵢeⱼ
/// - Pair orbits: 4C(d,2) points at (±λᵢ,±λⱼ,0,...)
/// - Diagonal orbit: 2^d corner points
/// 
/// \section orbit_differences Fourth Differences
/// The fourth divided difference along axis i is:
/// \f$ \Delta^4_i = |f_{3i} - 2f_1 - 7(f_{2i} - 2f_1)| \f$
/// where 7 = (λ₃/λ₂)² is the scaling factor from van Dooren & de Ridder.
/// 
/// \references
/// - A.C. Genz and A.A. Malik, "An adaptive algorithm for numerical integration
///   over an N-dimensional rectangular region", J. Comp. Appl. Math. 6 (1980)
/// - P. van Dooren and L. de Ridder, "An adaptive algorithm for numerical
///   integration over an n-dimensional cube", J. Comp. Appl. Math. 2 (1976)

/// \brief Orbit-based function evaluations for Genz-Malik rules
/// \tparam Real Floating point type
template <typename Real>
class orbit_evaluator {
public:
  
  /// \brief Cached evaluations organized by orbit structure
  /// \tparam Dim Dimension (compile-time for efficiency)
  /// \details Stores function values at orbit representatives to enable
  ///          efficient computation of both integration estimates and
  ///          directional derivatives for axis selection.
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
  /// \tparam F Integrand type
  /// \tparam Dim Dimension
  /// \tparam Degree Rule degree (5, 7, or 9)
  /// \param f Integrand function
  /// \param reg Integration region
  /// \param values Output: function values at orbit points
  /// 
  /// \details Evaluates the integrand at all orbit representatives needed
  ///          for the specified rule degree. Uses symmetry to minimize
  ///          the number of function evaluations.
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
    evaluate_axis_l3_points(f, center, half_widths, lambda3, values, evals,
                            std::integral_constant<bool, (Degree >= 9)>());
    
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
  /// \param values Function values at orbit points
  /// \return Array of fourth differences for each dimension
  /// 
  /// \details Implements the formula: \f$ \Delta^4_i = |f_{3i} - 2f_1 - 7(f_{2i} - 2f_1)| \f$
  /// where:
  /// - \f$ f_1 \f$ = center value
  /// - \f$ f_{2i} \f$ = sum at ±λ₂eᵢ
  /// - \f$ f_{3i} \f$ = sum at ±λ₃eᵢ (or ±λ₁eᵢ for degree-7)
  /// - 7 = (λ₃/λ₂)² = (9/10)/(9/70) = 7
  /// 
  /// The fourth difference estimates the magnitude of the fourth derivative
  /// along each axis, indicating where the integrand varies most rapidly.
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
  /// \param fourth_diffs Fourth differences for each dimension
  /// \param reg Current region
  /// \return Best dimension to split along
  /// 
  /// \details Selects the dimension with largest weighted fourth difference,
  ///          where weighting by axis width prefers splitting longer axes.
  ///          This heuristic tends to equalize error contributions across
  ///          dimensions and maintains aspect ratio.
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
  // Helper function for degree-9 axis points (SFINAE dispatch)
  template <typename F, std::size_t Dim>
  static void evaluate_axis_l3_points(
      const F& f,
      const std::array<Real, Dim>& center,
      const std::array<Real, Dim>& half_widths,
      Real lambda3,
      orbit_values<Dim>& values,
      std::size_t& evals,
      std::true_type) // Degree >= 9
  {
    std::array<Real, Dim> point;
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
  
  template <typename F, std::size_t Dim>
  static void evaluate_axis_l3_points(
      const F&, const std::array<Real, Dim>&,
      const std::array<Real, Dim>&, Real,
      orbit_values<Dim>&, std::size_t&,
      std::false_type) // Degree < 9
  {
    // No-op for degree < 9
  }
  
  // SFINAE helper for function signature detection
  template <typename F, typename Point, typename = void>
  struct function_signature {
    static constexpr int value = 0; // vector
  };
  
  template <typename F, typename Point>
  struct function_signature<F, Point,
      typename std::enable_if<
        std::is_convertible<decltype(std::declval<F>()(std::declval<const Real*>(), std::declval<std::size_t>())), Real>::value
      >::type> {
    static constexpr int value = 1; // pointer + size
  };
  
  template <typename F, typename Point>
  struct function_signature<F, Point,
      typename std::enable_if<
        !std::is_convertible<decltype(std::declval<F>()(std::declval<const Real*>(), std::declval<std::size_t>())), Real>::value &&
        std::is_convertible<decltype(std::declval<F>()(std::declval<const Real*>())), Real>::value
      >::type> {
    static constexpr int value = 2; // pointer only
  };
  
  // Tag dispatch functions
  template <typename F, std::size_t Dim>
  static Real evaluate_at_dispatch(const F& f, const std::array<Real, Dim>& point, std::integral_constant<int, 0>) {
    std::vector<Real> vec(point.begin(), point.end());
    return f(vec);
  }
  
  template <typename F, std::size_t Dim>
  static Real evaluate_at_dispatch(const F& f, const std::array<Real, Dim>& point, std::integral_constant<int, 1>) {
    return f(point.data(), Dim);
  }
  
  template <typename F, std::size_t Dim>
  static Real evaluate_at_dispatch(const F& f, const std::array<Real, Dim>& point, std::integral_constant<int, 2>) {
    return f(point.data());
  }
  
  /// \brief Helper to evaluate function with proper signature handling
  template <typename F, std::size_t Dim>
  static Real evaluate_at(const F& f, const std::array<Real, Dim>& point) {
    return evaluate_at_dispatch(f, point, 
        std::integral_constant<int, function_signature<F, decltype(point)>::value>());
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_ORBIT_EVALUATOR_HPP