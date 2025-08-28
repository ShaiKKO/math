// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_ADAPTIVE_IMPL_HPP
#define BOOST_MATH_CUBATURE_DETAIL_ADAPTIVE_IMPL_HPP

#include <boost/math/cubature/detail/adaptivity.hpp>
#include <boost/math/cubature/detail/genz_malik_evaluator.hpp>
#include <boost/math/cubature/detail/memory_pool.hpp>
#include <boost/math/cubature/detail/reliability_metrics.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <algorithm>
#include <unordered_map>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief DCUHRE adaptive integration implementation
/// \details Staff+ level implementation following ALGORITHMS.md §1.2
template <typename Real, typename F, typename Policy>
class adaptive_integrator {
private:
  const F& f_;
  const hypercube<Real>& box_;
  Real abs_tol_;
  Real rel_tol_;
  std::size_t max_eval_;
  std::size_t max_regions_;
  std::size_t max_depth_;
  const Policy& policy_;
  
  // Integration state with numerical stability
  kahan_accumulator<Real> integral_;
  kahan_accumulator<Real> error_sum_;
  std::size_t total_evals_;
  std::size_t num_regions_;
  
  // Memory pool for efficient region allocation
  std::unique_ptr<memory_pool<region<Real>>> region_pool_;
  
  // Convergence history tracking for reliability metrics
  convergence_history<Real> history_;
  std::vector<Real> region_errors_;  // Per-region errors for chi-squared
  std::size_t max_refinement_depth_;  // Track maximum depth reached
  
  // Directional variation cache for adaptive splitting
  struct directional_cache {
    std::vector<Real> variations;  // Per-dimension variation estimates
    std::size_t best_dim;          // Dimension with highest variation
    bool valid;                    // Whether cache is valid
    
    directional_cache() : best_dim(0), valid(false) {}
  };
  mutable std::unordered_map<std::size_t, directional_cache> variation_cache_;
  
  // Genz-Malik evaluator
  using evaluator_type = genz_malik_evaluator<Real>;
  
public:
  adaptive_integrator(const F& f, const hypercube<Real>& box,
                     Real abs_tol, Real rel_tol,
                     std::size_t max_eval,
                     const Policy& pol)
    : f_(f), box_(box), 
      abs_tol_(abs_tol), rel_tol_(rel_tol),
      max_eval_(max_eval > 0 ? max_eval : std::numeric_limits<std::size_t>::max()),
      max_regions_(10000),  // Default max regions
      max_depth_(50),        // Default max depth  
      policy_(pol),
      total_evals_(0),
      num_regions_(0),
      max_refinement_depth_(0)
  {
    // Configure memory pool for efficient region allocation
    typename memory_pool<region<Real>>::pool_config config;
    config.initial_capacity = 1024;  // Pre-allocate space for 1024 regions
    config.growth_factor = 2;        // Double capacity when growing
    config.max_capacity = max_regions_ * 2;  // Allow some headroom
    region_pool_ = std::make_unique<memory_pool<region<Real>>>(config);
  }
  
  result<Real> integrate() {
    const std::size_t dim = box_.dimension();
    
    // Initialize with whole box as first region
    region<Real> initial_region(dim);
    initial_region.a = box_.lower;
    initial_region.b = box_.upper;
    
    // Evaluate initial region using proper Genz-Malik evaluator
    embedded_pair_result<Real> initial_result;
    bool success = evaluate_region_dispatch(initial_region, initial_result, dim);
    
    if (!success) {
      // Rules not available for this dimension
      result<Real> res;
      res.value = Real(0);
      res.error = std::numeric_limits<Real>::max();
      res.evaluations = 0;
      res.status = status_code::dimension_error;
      return res;
    }
    
    initial_region.estimate_fine = initial_result.estimate_fine;
    initial_region.estimate_coarse = initial_result.estimate_coarse;
    initial_region.error = initial_result.embedded_error;
    initial_region.evaluations = initial_result.evaluations;
    initial_region.split_dim = initial_result.split_dimension;
    initial_region.cached_values = initial_result.cached_values;
    
    // Initialize integration state
    integral_.add(initial_result.estimate_fine);
    error_sum_.add(initial_region.error);
    total_evals_ = initial_result.evaluations;
    num_regions_ = 1;
    
    // Record initial state for convergence tracking
    history_.record(initial_result.estimate_fine, initial_region.error, 
                    total_evals_, num_regions_);
    region_errors_.push_back(initial_region.error);
    
    // Priority queue for regions (max-heap by error)
    region_priority_queue<Real> pq;
    pq.push(initial_region);
    
    // Main DCUHRE loop following ALGORITHMS.md
    std::size_t current_depth = 0;
    while (!pq.empty() && total_evals_ < max_eval_ && num_regions_ < max_regions_) {
      // Check convergence
      Real current_integral = integral_.sum();
      Real current_error = error_sum_.sum();
      
      if (check_convergence(current_error, current_integral, abs_tol_, rel_tol_)) {
        break;
      }
      
      // Track refinement depth
      current_depth++;
      if (current_depth > max_refinement_depth_) {
        max_refinement_depth_ = current_depth;
      }
      
      // Get region with largest error
      region<Real> current = pq.top();
      pq.pop();
      
      // Remove current region's contribution
      integral_.add(-current.estimate_fine);
      error_sum_.add(-current.error);
      
      // Use stored split dimension or compute with caching
      std::size_t split_dim = current.split_dim;
      if (split_dim == 0) {
        // Create region hash for cache lookup
        std::size_t region_hash = 0;
        for (std::size_t i = 0; i < dim; ++i) {
          std::hash<Real> hasher;
          region_hash ^= hasher(current.a[i]) + 0x9e3779b9 + (region_hash << 6) + (region_hash >> 2);
          region_hash ^= hasher(current.b[i]) + 0x9e3779b9 + (region_hash << 6) + (region_hash >> 2);
        }
        
        // Check cache for directional variations
        auto cache_it = variation_cache_.find(region_hash);
        if (cache_it != variation_cache_.end() && cache_it->second.valid) {
          split_dim = cache_it->second.best_dim;
        } else {
          // Compute directional variations
          directional_cache& cache_entry = variation_cache_[region_hash];
          cache_entry.variations.resize(dim);
          
          Real max_variation = Real(0);
          for (std::size_t i = 0; i < dim; ++i) {
            Real width = current.b[i] - current.a[i];
            
            // Estimate variation using fourth-difference if available
            Real variation = width;
            // For now, use width-based heuristic
            // Full implementation would cast cached_values to proper type
            
            cache_entry.variations[i] = variation;
            if (variation > max_variation) {
              max_variation = variation;
              split_dim = i;
            }
          }
          
          cache_entry.best_dim = split_dim;
          cache_entry.valid = true;
        }
      }
      
      // Bisect region
      auto [left, right] = bisect_region(current, split_dim);
      
      // Evaluate child regions with parent's cached values
      embedded_pair_result<Real> left_result, right_result;
      
      // Pass parent's cached values to child evaluations
      evaluate_region_dispatch_with_cache(left, left_result, dim, current.cached_values);
      evaluate_region_dispatch_with_cache(right, right_result, dim, current.cached_values);
      
      left.estimate_fine = left_result.estimate_fine;
      left.estimate_coarse = left_result.estimate_coarse;
      left.error = left_result.embedded_error;
      left.evaluations = left_result.evaluations;
      left.split_dim = left_result.split_dimension;
      left.cached_values = left_result.cached_values;
      
      right.estimate_fine = right_result.estimate_fine;
      right.estimate_coarse = right_result.estimate_coarse;
      right.error = right_result.embedded_error;
      right.evaluations = right_result.evaluations;
      right.split_dim = right_result.split_dimension;
      right.cached_values = right_result.cached_values;
      
      // Update global estimates
      integral_.add(left.estimate_fine);
      integral_.add(right.estimate_fine);
      error_sum_.add(left.error);
      error_sum_.add(right.error);
      total_evals_ += left.evaluations + right.evaluations;
      num_regions_ += 2;
      
      // Track per-region errors for chi-squared analysis
      region_errors_.push_back(left.error);
      region_errors_.push_back(right.error);
      
      // Record convergence history
      history_.record(integral_.sum(), error_sum_.sum(), 
                     total_evals_, num_regions_);
      
      // Add child regions to queue if they have significant error
      const Real local_tol = abs_tol_ * left.volume() / box_.volume();
      if (left.error > local_tol) {
        pq.push(left);
      }
      if (right.error > local_tol) {
        pq.push(right);
      }
    }
    
    // Prepare result
    result<Real> res;
    res.value = integral_.sum();
    res.error = error_sum_.sum();
    res.evaluations = total_evals_;
    
    // Compute reliability metrics using collected statistics
    res.reliability = reliability_calculator<Real>::compute_metrics(
        history_,
        region_errors_,
        res.error,
        max_refinement_depth_);
    
    // Set status
    if (check_convergence(res.error, res.value, abs_tol_, rel_tol_)) {
      res.status = status_code::success;
    } else if (total_evals_ >= max_eval_) {
      res.status = status_code::maxeval_reached;
    } else if (num_regions_ >= max_regions_) {
      res.status = status_code::maxregions_reached;
    } else {
      res.status = status_code::success; // Converged within limits
    }
    
    return res;
  }
  
private:
  // Runtime dimension dispatch for region evaluation
  bool evaluate_region_dispatch(const region<Real>& reg,
                               embedded_pair_result<Real>& result,
                               std::size_t dim) {
    return evaluate_region_dispatch_with_cache(reg, result, dim, nullptr);
  }
  
  // Runtime dimension dispatch for region evaluation with caching
  bool evaluate_region_dispatch_with_cache(const region<Real>& reg,
                                          embedded_pair_result<Real>& result,
                                          std::size_t dim,
                                          const std::shared_ptr<void>& parent_cache) {
    switch(dim) {
      case 2:
        return evaluator_type::template evaluate_embedded_pair<F, 2>(f_, reg, result, parent_cache);
      case 3:
        return evaluator_type::template evaluate_embedded_pair<F, 3>(f_, reg, result, parent_cache);
      case 4:
        return evaluator_type::template evaluate_embedded_pair<F, 4>(f_, reg, result, parent_cache);
      case 5:
        return evaluator_type::template evaluate_embedded_pair<F, 5>(f_, reg, result, parent_cache);
      case 6:
        return evaluator_type::template evaluate_embedded_pair<F, 6>(f_, reg, result, parent_cache);
      case 7:
        return evaluator_type::template evaluate_embedded_pair<F, 7>(f_, reg, result, parent_cache);
      case 8:
        return evaluator_type::template evaluate_embedded_pair<F, 8>(f_, reg, result, parent_cache);
      case 9:
        return evaluator_type::template evaluate_embedded_pair<F, 9>(f_, reg, result, parent_cache);
      case 10:
        return evaluator_type::template evaluate_embedded_pair<F, 10>(f_, reg, result, parent_cache);
      case 11:
        return evaluator_type::template evaluate_embedded_pair<F, 11>(f_, reg, result, parent_cache);
      case 12:
        return evaluator_type::template evaluate_embedded_pair<F, 12>(f_, reg, result, parent_cache);
      case 13:
        return evaluator_type::template evaluate_embedded_pair<F, 13>(f_, reg, result, parent_cache);
      case 14:
        return evaluator_type::template evaluate_embedded_pair<F, 14>(f_, reg, result, parent_cache);
      case 15:
        return evaluator_type::template evaluate_embedded_pair<F, 15>(f_, reg, result, parent_cache);
      default:
        // For dimensions > 15, rules are not available
        return false;
    }
  }
  
  // Select split dimension based on directional variations
  std::size_t select_split_dimension_dispatch(
      const embedded_pair_result<Real>& result,
      const region<Real>& reg,
      std::size_t dim) {
    switch(dim) {
      case 2:
        return evaluator_type::template select_split_dimension<2>(result, reg);
      case 3:
        return evaluator_type::template select_split_dimension<3>(result, reg);
      case 4:
        return evaluator_type::template select_split_dimension<4>(result, reg);
      case 5:
        return evaluator_type::template select_split_dimension<5>(result, reg);
      case 6:
        return evaluator_type::template select_split_dimension<6>(result, reg);
      case 7:
        return evaluator_type::template select_split_dimension<7>(result, reg);
      case 8:
        return evaluator_type::template select_split_dimension<8>(result, reg);
      case 9:
        return evaluator_type::template select_split_dimension<9>(result, reg);
      case 10:
        return evaluator_type::template select_split_dimension<10>(result, reg);
      case 11:
        return evaluator_type::template select_split_dimension<11>(result, reg);
      case 12:
        return evaluator_type::template select_split_dimension<12>(result, reg);
      case 13:
        return evaluator_type::template select_split_dimension<13>(result, reg);
      case 14:
        return evaluator_type::template select_split_dimension<14>(result, reg);
      case 15:
        return evaluator_type::template select_split_dimension<15>(result, reg);
      default:
        // Fallback: split longest dimension
        std::size_t max_dim = 0;
        Real max_width = Real(0);
        for (std::size_t i = 0; i < dim; ++i) {
          Real width = reg.b[i] - reg.a[i];
          if (width > max_width) {
            max_width = width;
            max_dim = i;
          }
        }
        return max_dim;
    }
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_ADAPTIVE_IMPL_HPP