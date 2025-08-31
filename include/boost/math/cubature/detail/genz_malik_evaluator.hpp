// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_GENZ_MALIK_EVALUATOR_HPP
#define BOOST_MATH_CUBATURE_DETAIL_GENZ_MALIK_EVALUATOR_HPP

#include <boost/math/cubature/constants.hpp>
#include <boost/math/cubature/detail/adaptivity.hpp>
#include <boost/math/cubature/detail/gm_rules.hpp>
#include <boost/math/cubature/detail/orbit_evaluator.hpp>
#include <boost/math/cubature/detail/function_cache.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <memory>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Result of evaluating embedded Genz-Malik rule pair
/// \details Contains both primary and embedded rule estimates for error estimation,
///          along with directional variation information for adaptive subdivision.
template <typename Real>
struct EmbeddedRulePairResult {
  Real primary_estimate;         ///< Higher-degree rule estimate (degree 7 or 9)
  Real embedded_estimate;        ///< Lower-degree rule estimate (degree 5 or 7)
  Real error_estimate;           ///< Absolute difference between rule estimates
  Real variation_bound;          ///< Max-min spread of partial sums for error safeguarding
  std::size_t function_evaluations; ///< Number of integrand evaluations used
  std::array<Real, 15> directional_fourth_differences; ///< Fourth differences for each axis
  std::size_t optimal_split_axis;   ///< Axis with maximum variation for subdivision
  
  // Performance optimization: cached evaluations for potential reuse
  std::shared_ptr<void> cached_orbit_values; ///< Type-erased storage for orbit points
  
  EmbeddedRulePairResult() 
    : primary_estimate(0), 
      embedded_estimate(0), 
      error_estimate(0),
      variation_bound(0), 
      function_evaluations(0), 
      optimal_split_axis(0) {
    std::fill(directional_fourth_differences.begin(), 
              directional_fourth_differences.end(), Real(0));
  }
  
  // Backward compatibility
  Real estimate_fine() const { return primary_estimate; }
  Real estimate_coarse() const { return embedded_estimate; }
  Real embedded_error() const { return error_estimate; }
  Real spread_estimate() const { return variation_bound; }
  std::size_t evaluations() const { return function_evaluations; }
  std::size_t split_dimension() const { return optimal_split_axis; }
};

// Maintain backward compatibility
template <typename Real>
using embedded_pair_result = EmbeddedRulePairResult<Real>;

// SFINAE helper for detecting integrand signature
template <typename F, typename Point, typename = void>
struct integrand_signature {
    static constexpr int value = 0; // default: vector
};

template <typename F, typename Point>
struct integrand_signature<F, Point,
    typename std::enable_if<
        std::is_convertible<decltype(std::declval<F>()(std::declval<const typename Point::value_type*>(), std::declval<std::size_t>())), typename Point::value_type>::value
    >::type> {
    static constexpr int value = 1; // pointer + size
};

template <typename F, typename Point>
struct integrand_signature<F, Point,
    typename std::enable_if<
        !std::is_convertible<decltype(std::declval<F>()(std::declval<const typename Point::value_type*>(), std::declval<std::size_t>())), typename Point::value_type>::value &&
        std::is_convertible<decltype(std::declval<F>()(std::declval<const typename Point::value_type*>())), typename Point::value_type>::value
    >::type> {
    static constexpr int value = 2; // pointer only
};

/// \brief Genz-Malik rule evaluator with proper table access
/// \details Implementation following Boost.Math patterns
template <typename Real>
class genz_malik_evaluator {
private:
  // Tag dispatch functions for evaluating integrand
  template <typename F, typename Point>
  static typename Point::value_type evaluate_integrand_dispatch(const F& f, const Point& point, std::integral_constant<int, 0>) {
    std::vector<typename Point::value_type> vec(point.begin(), point.end());
    return f(vec);
  }
  
  template <typename F, typename Point>
  static typename Point::value_type evaluate_integrand_dispatch(const F& f, const Point& point, std::integral_constant<int, 1>) {
    return f(point.data(), point.size());
  }
  
  template <typename F, typename Point>
  static typename Point::value_type evaluate_integrand_dispatch(const F& f, const Point& point, std::integral_constant<int, 2>) {
    return f(point.data());
  }
  
public:
  static constexpr int degree_fine = 9;
  static constexpr int degree_coarse = 7;
  
  /// \brief Evaluate embedded rule pair on a region with optional caching
  /// \param f The integrand function
  /// \param reg The region to integrate over
  /// \param result Output structure with estimates and diagnostics
  /// \param parent_cache Optional cached values from parent region
  /// \return true if evaluation succeeded, false if rules unavailable
  template <typename F, std::size_t Dim>
  static bool evaluate_embedded_pair(
      const F& f,
      const region<Real>& reg,
      embedded_pair_result<Real>& result,
      const std::shared_ptr<void>& parent_cache = nullptr) 
  {
    using namespace boost::math::cubature::detail::gm;
    
    // Check rule availability at compile time
    using rule_fine = rule_fam<degree_fine, Dim, family_9_7>;
    using rule_coarse = rule_fam<degree_coarse, Dim, family_9_7>;
    
    // Use SFINAE-based check instead of if constexpr
    return evaluate_embedded_pair_impl<rule_fine, rule_coarse, F, Dim>(f, reg, result, parent_cache,
        std::integral_constant<bool, rule_fine::available && rule_coarse::available>());
  }
  
private:
  // Implementation when rules are available
  template <typename RuleFine, typename RuleCoarse, typename F, std::size_t Dim>
  static bool evaluate_embedded_pair_impl(
      const F& f,
      const region<Real>& reg,
      embedded_pair_result<Real>& result,
      const std::shared_ptr<void>& parent_cache,
      std::true_type) // rules available
  {
    using namespace boost::math::cubature::detail::gm;
    
    // Get nodes and weights from actual tables
    const auto nodes = RuleFine::template nodes<Real>();
    const auto weights_fine = RuleFine::template weights<Real>();
    const auto weights_coarse = RuleCoarse::template weights<Real>();
    
    const std::size_t num_nodes = nodes.size();
    
    // Prepare accumulators with Kahan summation
    kahan_accumulator<Real> sum_fine, sum_coarse;
    
    // For spread estimation (safeguarding)
    Real max_partial = std::numeric_limits<Real>::lowest();
    Real min_partial = std::numeric_limits<Real>::max();
    
    // Store function values for axis selection and caching
    std::vector<Real> function_values;
    function_values.reserve(num_nodes);
    
    // Prepare cache for storing evaluations
    using cache_t = typename function_cache<Real, Dim>::parent_cache;
    auto new_cache = std::make_shared<cache_t>();
    new_cache->evaluations.reserve(num_nodes);
    
    // Extract reusable values from parent cache if available
    std::vector<typename function_cache<Real, Dim>::cached_evaluation> reusable_values;
    if (parent_cache) {
      try {
        auto parent_typed = std::static_pointer_cast<cache_t>(parent_cache);
        if (parent_typed && !parent_typed->evaluations.empty()) {
          // Parent's cached evaluations are in physical space,
          // so they can be directly reused if the nodes match
          reusable_values = parent_typed->evaluations;
        }
      } catch (const std::bad_cast&) {
        // Parent cache is not of expected type, proceed without caching
        reusable_values.clear();
      }
    }
    
    // Transform and evaluate at all nodes using actual rule structure
    std::array<Real, Dim> point;
    std::size_t reused_count = 0;
    
    for (std::size_t i = 0; i < num_nodes; ++i) {
      Real value;
      
      // Transform node from [0,1]^d to physical space
      transform_reference_node_to_region<Dim>(nodes[i], reg, point);
      
      // Check if we can reuse a cached value
      const Real* cached_value = nullptr;
      if (!reusable_values.empty()) {
        cached_value = function_cache<Real, Dim>::find_cached_value(
            reusable_values, point);  // Search using physical coordinates
      }
      
      if (cached_value) {
        // Reuse cached value
        value = *cached_value;
        reused_count++;
      } else {
        // Evaluate integrand using tag dispatch
        value = evaluate_integrand_dispatch(f, point, 
            std::integral_constant<int, integrand_signature<F, decltype(point)>::value>());
      }
      
      function_values.push_back(value);
      new_cache->evaluations.emplace_back(point, value);  // Store physical coordinates
      
      // Accumulate for fine rule
      Real contrib_fine = value * weights_fine[i];
      sum_fine.add(contrib_fine);
      
      // Track spread for safeguarding
      Real partial = sum_fine.sum();
      max_partial = std::max(max_partial, partial);
      min_partial = std::min(min_partial, partial);
      
      // Accumulate for coarse rule (subset of nodes)
      if (i < weights_coarse.size()) {
        sum_coarse.add(value * weights_coarse[i]);
      }
    }
    
    // Compute Jacobian for transformation
    Real jacobian = reg.volume();
    
    // Final estimates
    result.primary_estimate = sum_fine.sum() * jacobian;
    result.embedded_estimate = sum_coarse.sum() * jacobian;
    result.error_estimate = std::abs(result.primary_estimate - result.embedded_estimate);
    result.variation_bound = (max_partial - min_partial) * jacobian;
    result.function_evaluations = num_nodes - reused_count;  // Only count new evaluations
    
    // Compute fourth differences from the evaluated points for axis selection
    compute_directional_fourth_differences<Dim>(nodes, function_values, reg, result);
    
    // Store cache for potential reuse by child regions
    result.cached_orbit_values = new_cache;
    
    // Apply safeguarding
    apply_error_safeguarding(result);
    
    return true;
  }
  
  // Implementation when rules are not available
  template <typename RuleFine, typename RuleCoarse, typename F, std::size_t Dim>
  static bool evaluate_embedded_pair_impl(
      const F&,
      const region<Real>&,
      embedded_pair_result<Real>&,
      const std::shared_ptr<void>&,
      std::false_type) // rules not available
  {
    return false;
  }
  
  /// \brief Get split dimension from result (already computed)
  template <std::size_t Dim>
  static std::size_t select_optimal_split_axis(
      const embedded_pair_result<Real>& result,
      const region<Real>& /*reg*/)
  {
    return result.optimal_split_axis;
  }
  
private:
  /// \brief Transform node from [0,1]^d reference to region
  template <std::size_t Dim>
  static void transform_reference_node_to_region(
      const std::array<Real, Dim>& ref_node,
      const region<Real>& reg,
      std::array<Real, Dim>& point)
  {
    for (std::size_t i = 0; i < Dim; ++i) {
      point[i] = reg.a[i] + ref_node[i] * (reg.b[i] - reg.a[i]);
    }
  }
  
  /// \brief Compute fourth differences and select split dimension
  /// \details Uses orbit structure to compute proper fourth differences
  ///          Following van Dooren and de Ridder formula
  template <std::size_t Dim, typename NodesArray>
  static void compute_directional_fourth_differences(
      const NodesArray& nodes,
      const std::vector<Real>& values,
      const region<Real>& reg,
      embedded_pair_result<Real>& result)
  {
    // Initialize fourth differences
    std::fill(result.directional_fourth_differences.begin(), 
              result.directional_fourth_differences.begin() + Dim, Real(0));
    
    // Extract orbit values from the evaluated nodes
    // The nodes are ordered according to the Genz-Malik rule structure
    
    // Identify key orbit representatives based on rule structure
    // Center point is always first
    Real f_center = values[0];
    
    // Extract axis points (λ₁, λ₂, λ₃) by identifying nodes that differ in only one dimension
    std::array<Real, Dim> f_axis_l1_plus, f_axis_l1_minus;
    std::array<Real, Dim> f_axis_l2_plus, f_axis_l2_minus;
    std::array<Real, Dim> f_axis_l3_plus, f_axis_l3_minus;
    
    // Initialize arrays
    std::fill(f_axis_l1_plus.begin(), f_axis_l1_plus.end(), Real(0));
    std::fill(f_axis_l1_minus.begin(), f_axis_l1_minus.end(), Real(0));
    std::fill(f_axis_l2_plus.begin(), f_axis_l2_plus.end(), Real(0));
    std::fill(f_axis_l2_minus.begin(), f_axis_l2_minus.end(), Real(0));
    std::fill(f_axis_l3_plus.begin(), f_axis_l3_plus.end(), Real(0));
    std::fill(f_axis_l3_minus.begin(), f_axis_l3_minus.end(), Real(0));
    
    // Identify axis points by checking which nodes differ from center in exactly one dimension
    for (std::size_t i = 1; i < nodes.size(); ++i) {
      std::size_t diff_dim = Dim;  // Invalid dimension
      std::size_t diff_count = 0;
      Real diff_val = 0;
      
      for (std::size_t d = 0; d < Dim; ++d) {
        Real delta = nodes[i][d] - nodes[0][d];  // Difference from center
        if (std::abs(delta) > 1e-10) {
          diff_count++;
          diff_dim = d;
          diff_val = delta;
        }
      }
      
      // If exactly one dimension differs, it's an axis point
      if (diff_count == 1 && diff_dim < Dim) {
        // Classify by distance (λ₁, λ₂, or λ₃)
        Real abs_diff = std::abs(diff_val);
        Real half_width = (reg.b[diff_dim] - reg.a[diff_dim]) / 2;
        Real normalized_diff = abs_diff / half_width;
        
        // Approximate λ values from Genz-Malik rules
        const Real lambda1 = std::sqrt(Real(9)/Real(70));  // ≈ 0.358
        const Real lambda2 = std::sqrt(Real(9)/Real(10));  // ≈ 0.949
        const Real lambda3 = Real(1);                       // 1.0
        
        if (std::abs(normalized_diff - lambda1) < 0.01) {
          // λ₁ point
          if (diff_val > 0) {
            f_axis_l1_plus[diff_dim] = values[i];
          } else {
            f_axis_l1_minus[diff_dim] = values[i];
          }
        } else if (std::abs(normalized_diff - lambda2) < 0.01) {
          // λ₂ point
          if (diff_val > 0) {
            f_axis_l2_plus[diff_dim] = values[i];
          } else {
            f_axis_l2_minus[diff_dim] = values[i];
          }
        } else if (std::abs(normalized_diff - lambda3) < 0.01) {
          // λ₃ point (corners)
          if (diff_val > 0) {
            f_axis_l3_plus[diff_dim] = values[i];
          } else {
            f_axis_l3_minus[diff_dim] = values[i];
          }
        }
      }
    }
    
    // Compute fourth differences using the formula
    // f₃ᵢ - 2f₁ - 7*(f₂ᵢ - 2f₁) where 7 = (λ₃/λ₂)²
    const Real scale_factor = Real(7);
    
    for (std::size_t d = 0; d < Dim; ++d) {
      Real f1 = f_center;
      Real f2i = f_axis_l2_plus[d] + f_axis_l2_minus[d];
      
      // Use λ₃ points if available, otherwise use λ₁ as fallback
      Real f3i;
      if (f_axis_l3_plus[d] != Real(0) || f_axis_l3_minus[d] != Real(0)) {
        f3i = f_axis_l3_plus[d] + f_axis_l3_minus[d];
      } else {
        f3i = f_axis_l1_plus[d] + f_axis_l1_minus[d];
      }
      
      // Fourth divided difference
      Real diff = std::abs(f3i - Real(2)*f1 - scale_factor*(f2i - Real(2)*f1));
      result.directional_fourth_differences[d] = diff;
    }
    
    // Select dimension with largest weighted fourth difference
    std::size_t best_dim = 0;
    Real max_weighted_diff = 0;
    
    for (std::size_t d = 0; d < Dim; ++d) {
      Real width = reg.b[d] - reg.a[d];
      // Weight by width^5 to account for transformation scaling
      Real weighted_diff = result.directional_fourth_differences[d] * std::pow(width, Real(5));
      
      if (weighted_diff > max_weighted_diff) {
        max_weighted_diff = weighted_diff;
        best_dim = d;
      }
    }
    
    result.optimal_split_axis = best_dim;
  }
  
  
  /// \brief Apply safeguarding to prevent spuriously small errors
  static void apply_error_safeguarding(embedded_pair_result<Real>& result)
  {
    const Real machine_eps = std::numeric_limits<Real>::epsilon();
    const Real safety_factor = Real(boost::math::cubature::constants::dcuhre_empirical_safety_factor);
    
    // If embedded difference is too small relative to spread, use spread estimate
    if (result.error_estimate < safety_factor * machine_eps * std::abs(result.primary_estimate)) {
      if (result.variation_bound > result.error_estimate) {
        result.error_estimate = result.variation_bound;
      }
    }
  }
  
  /// \brief Compute node-based spread statistics for reliability assessment
  template <std::size_t Dim>
  static Real compute_function_value_spread(const std::vector<Real>& values)
  {
    if (values.size() < 2) {
      return Real(0);
    }
    
    // Compute min, max, and standard deviation
    Real min_val = *std::min_element(values.begin(), values.end());
    Real max_val = *std::max_element(values.begin(), values.end());
    Real spread = max_val - min_val;
    
    // Also compute standard deviation for additional reliability info
    Real mean = std::accumulate(values.begin(), values.end(), Real(0)) / values.size();
    Real variance = Real(0);
    for (const auto& val : values) {
      Real diff = val - mean;
      variance += diff * diff;
    }
    variance /= values.size();
    Real std_dev = std::sqrt(variance);
    
    // Return the maximum of raw spread and 3*std_dev as conservative estimate
    return std::max(spread, Real(3) * std_dev);
  }
};

/// \brief Specialized evaluator for vector integrands
template <typename Real>
class genz_malik_vector_evaluator {
public:
  /// \brief Evaluate embedded pair for vector integrand
  /// \param f Vector integrand f(x, out, m) where m is number of components
  /// \param reg Region to integrate
  /// \param num_components Number of integrand components
  /// \param results Vector of results, one per component
  template <typename F, std::size_t Dim>
  static bool evaluate_embedded_pair_vector(
      const F& f,
      const region<Real>& reg,
      std::size_t num_components,
      std::vector<embedded_pair_result<Real>>& results)
  {
    using namespace boost::math::cubature::detail::gm;
    
    // Check rule availability
    using rule_fine = rule_fam<9, Dim, family_9_7>;
    using rule_coarse = rule_fam<7, Dim, family_9_7>;
    
    // Use runtime check for availability
    if (!rule_fine::available || !rule_coarse::available) {
      return false;
    }
    
    // Get nodes and weights
    const auto nodes = rule_fine::template nodes<Real>();
    const auto weights_fine = rule_fine::template weights<Real>();
    const auto weights_coarse = rule_coarse::template weights<Real>();
    
    // Prepare results vector
    results.resize(num_components);
    
    // Accumulators for each component
    std::vector<kahan_accumulator<Real>> sums_fine(num_components);
    std::vector<kahan_accumulator<Real>> sums_coarse(num_components);
    std::vector<Real> max_partials(num_components, std::numeric_limits<Real>::lowest());
    std::vector<Real> min_partials(num_components, std::numeric_limits<Real>::max());
    
    // Temporary storage for function values
    std::vector<Real> values(num_components);
    
    // Evaluate at all nodes
    std::array<Real, Dim> point;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      // Transform node
      for (std::size_t d = 0; d < Dim; ++d) {
        point[d] = reg.a[d] + nodes[i][d] * (reg.b[d] - reg.a[d]);
      }
      
      // Evaluate all components
      f(point.data(), values.data(), num_components);
      
      // Accumulate for each component
      for (std::size_t c = 0; c < num_components; ++c) {
        Real contrib_fine = values[c] * weights_fine[i];
        sums_fine[c].add(contrib_fine);
        
        // Track spread
        Real partial = sums_fine[c].sum();
        max_partials[c] = std::max(max_partials[c], partial);
        min_partials[c] = std::min(min_partials[c], partial);
        
        // Coarse rule
        if (i < weights_coarse.size()) {
          sums_coarse[c].add(values[c] * weights_coarse[i]);
        }
        
        // Cache values for first component (used for splitting decisions)
        if (c == 0) {
          results[0].function_values.push_back(values[c]);
        }
      }
    }
    
    // Compute final estimates for each component
    Real jacobian = reg.volume();
    for (std::size_t c = 0; c < num_components; ++c) {
      results[c].estimate_fine = sums_fine[c].sum() * jacobian;
      results[c].estimate_coarse = sums_coarse[c].sum() * jacobian;
      results[c].embedded_error = std::abs(results[c].estimate_fine - results[c].estimate_coarse);
      results[c].spread_estimate = (max_partials[c] - min_partials[c]) * jacobian;
      results[c].evaluations = nodes.size();
      
      // Apply safeguarding
      safeguard_error(results[c]);
    }
    
    return true;
  }
  
private:
  static void apply_error_safeguarding(embedded_pair_result<Real>& result) {
    genz_malik_evaluator<Real>::safeguard_error(result);
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_GENZ_MALIK_EVALUATOR_HPP