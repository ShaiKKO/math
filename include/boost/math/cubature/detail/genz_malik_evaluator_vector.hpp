// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_GENZ_MALIK_EVALUATOR_VECTOR_HPP
#define BOOST_MATH_CUBATURE_DETAIL_GENZ_MALIK_EVALUATOR_VECTOR_HPP

#include <boost/math/cubature/detail/gm_rules.hpp>
#include <boost/math/cubature/detail/adaptivity.hpp>
#include <boost/math/cubature/detail/vector_adapter.hpp>
#include <boost/math/cubature/regions.hpp>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

namespace boost { namespace math { namespace cubature { namespace detail {

template <typename Real>
struct embedded_pair_result_vector {
  std::vector<Real> estimates_fine;     
  std::vector<Real> estimates_coarse;   
  std::vector<Real> embedded_errors;    
  std::size_t evaluations;
  std::size_t split_dimension;          
  
  embedded_pair_result_vector() 
    : evaluations(0), split_dimension(0) {}
};

template <typename Real>
class genz_malik_evaluator_vector {
public:
  static constexpr int degree_fine = 9;
  static constexpr int degree_coarse = 7;
  
  template <typename F, std::size_t Dim>
  static bool evaluate_embedded_pair(
      const F& f,
      const region<Real>& reg,
      embedded_pair_result_vector<Real>& result,
      std::size_t num_components,
      vector_workspace<Real>* workspace = nullptr) 
  {
    using namespace boost::math::cubature::detail::gm;
    
    // Use raw rules for proper embedded pairs (includes zero-weight nodes)
    using rule_fine = raw_rule_fam<degree_fine, Dim, family_9_7>;
    using rule_coarse = raw_rule_fam<degree_coarse, Dim, family_9_7>;
    
    if constexpr (!rule_fine::available || !rule_coarse::available) {
      return false;
    }
    
    // Get nodes and weights from raw rules (both have same size)
    const auto nodes = rule_fine::template nodes_with_zero<Real>();
    const auto weights_fine = rule_fine::template weights_with_zero<Real>();
    const auto weights_coarse = rule_coarse::template weights_with_zero<Real>();
    
    const std::size_t num_nodes = nodes.size();
    
    // Initialize result vectors
    result.estimates_fine.assign(num_components, Real(0));
    result.estimates_coarse.assign(num_components, Real(0));
    result.embedded_errors.assign(num_components, Real(0));
    result.evaluations = num_nodes;
    
    // Prepare Kahan accumulators for each component
    std::vector<kahan_accumulator<Real>> sums_fine(num_components);
    std::vector<kahan_accumulator<Real>> sums_coarse(num_components);
    
    // Fourth differences for axis selection (per dimension)
    std::vector<std::vector<Real>> fourth_diffs(Dim, std::vector<Real>(num_components, Real(0)));
    
    // Temporary buffer for function evaluation
    std::vector<Real> temp_values(num_components);
    
    // Track orbit evaluations for fourth difference calculation
    std::vector<std::vector<Real>> orbit_values;
    orbit_values.reserve(5);  // Max 5 orbits in GM rules
    
    std::size_t orbit_idx = 0;
    std::size_t prev_orbit_size = 0;
    
    // Evaluate at all nodes
    for (std::size_t i = 0; i < num_nodes; ++i) {
      // Transform node to region
      std::array<Real, Dim> x;
      for (std::size_t d = 0; d < Dim; ++d) {
        x[d] = reg.a[d] + (reg.b[d] - reg.a[d]) * nodes[i][d];
      }
      
      // Evaluate function (single call for all components)
      f(x.data(), temp_values.data(), num_components);
      
      // Store for orbit analysis
      if (i == 0 || (i > prev_orbit_size && orbit_idx < 5)) {
        orbit_values.push_back(temp_values);
        prev_orbit_size = i;
        ++orbit_idx;
      }
      
      // Accumulate weighted sums for each component
      Real w_fine = weights_fine[i] * reg.volume();
      Real w_coarse = weights_coarse[i] * reg.volume();
      
      for (std::size_t c = 0; c < num_components; ++c) {
        sums_fine[c].add(temp_values[c] * w_fine);
        sums_coarse[c].add(temp_values[c] * w_coarse);
      }
    }
    
    // Extract results from accumulators
    for (std::size_t c = 0; c < num_components; ++c) {
      result.estimates_fine[c] = sums_fine[c].sum();
      result.estimates_coarse[c] = sums_coarse[c].sum();
      result.embedded_errors[c] = std::abs(result.estimates_fine[c] - result.estimates_coarse[c]);
    }
    
    // Calculate fourth differences for axis selection
    if (orbit_values.size() >= 3 && Dim > 1) {
      // Using simplified fourth difference formula
      // Based on orbit structure: center, face points, edge points
      const Real lambda_ratio = Real(7);  // (λ₃/λ₂)² for degree 7
      
      for (std::size_t d = 0; d < Dim; ++d) {
        for (std::size_t c = 0; c < num_components; ++c) {
          if (orbit_values.size() > 2) {
            // Simplified: use difference between orbit evaluations
            Real diff = std::abs(orbit_values[2][c] - orbit_values[0][c]);
            fourth_diffs[d][c] = diff * std::pow(reg.b[d] - reg.a[d], 5);
          }
        }
      }
    }
    
    // Select split dimension based on maximum fourth difference across all components
    Real max_diff = Real(0);
    result.split_dimension = 0;
    
    for (std::size_t d = 0; d < Dim; ++d) {
      Real dim_max = Real(0);
      for (std::size_t c = 0; c < num_components; ++c) {
        dim_max = std::max(dim_max, fourth_diffs[d][c]);
      }
      if (dim_max > max_diff) {
        max_diff = dim_max;
        result.split_dimension = d;
      }
    }
    
    return true;
  }
  
  // Simplified version for dimensions where GM rules are available
  template <typename F>
  static bool evaluate_embedded_pair_runtime(
      const F& f,
      const region<Real>& reg,
      embedded_pair_result_vector<Real>& result,
      std::size_t num_components,
      vector_workspace<Real>* workspace = nullptr) 
  {
    const std::size_t dim = reg.dimension();
    
    // Dispatch based on dimension
    switch (dim) {
      case 2: return evaluate_embedded_pair<F, 2>(f, reg, result, num_components, workspace);
      case 3: return evaluate_embedded_pair<F, 3>(f, reg, result, num_components, workspace);
      case 4: return evaluate_embedded_pair<F, 4>(f, reg, result, num_components, workspace);
      case 5: return evaluate_embedded_pair<F, 5>(f, reg, result, num_components, workspace);
      case 6: return evaluate_embedded_pair<F, 6>(f, reg, result, num_components, workspace);
      case 7: return evaluate_embedded_pair<F, 7>(f, reg, result, num_components, workspace);
      case 8: return evaluate_embedded_pair<F, 8>(f, reg, result, num_components, workspace);
      case 9: return evaluate_embedded_pair<F, 9>(f, reg, result, num_components, workspace);
      case 10: return evaluate_embedded_pair<F, 10>(f, reg, result, num_components, workspace);
      default: return false;  // GM rules not available for this dimension
    }
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_GENZ_MALIK_EVALUATOR_VECTOR_HPP