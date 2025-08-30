// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_ADAPTIVE_INTEGRATOR_VECTOR_HPP
#define BOOST_MATH_CUBATURE_DETAIL_ADAPTIVE_INTEGRATOR_VECTOR_HPP

#include <boost/math/cubature/detail/vector_adapter.hpp>
#include <boost/math/cubature/detail/genz_malik_evaluator_vector.hpp>
#include <boost/math/cubature/detail/adaptivity.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/policies.hpp>
#include <vector>
#include <queue>
#include <memory>
#include <algorithm>

namespace boost { namespace math { namespace cubature { namespace detail {

template <typename Real, typename F, typename Policy>
class adaptive_integrator_vector {
private:
  const F& f_;
  const hypercube<Real>& box_;
  Real abs_tol_;
  Real rel_tol_;
  std::size_t max_eval_;
  std::size_t num_components_;
  error_norm norm_type_;
  Policy pol_;
  
  // Vector-specific workspace
  std::unique_ptr<vector_workspace<Real>> workspace_;
  std::unique_ptr<vector_result_aggregator<Real>> aggregator_;
  
  // Priority queue of regions
  using vector_region_ptr = std::shared_ptr<vector_region<Real>>;
  struct region_comparator {
    bool operator()(const vector_region_ptr& a, const vector_region_ptr& b) const {
      return a->aggregated_error < b->aggregated_error;
    }
  };
  std::priority_queue<vector_region_ptr, std::vector<vector_region_ptr>, region_comparator> region_queue_;
  
  // Tracking
  std::size_t total_evaluations_;
  reliability_metrics<Real> metrics_;
  
public:
  adaptive_integrator_vector(const F& f, const hypercube<Real>& box,
                            Real abs_tol, Real rel_tol, 
                            std::size_t num_components,
                            std::size_t max_eval = 0,
                            error_norm norm = error_norm::l_infinity,
                            Policy const& pol = Policy{})
    : f_(f), box_(box), abs_tol_(abs_tol), rel_tol_(rel_tol),
      num_components_(num_components),
      max_eval_(max_eval > 0 ? max_eval : 100000),
      norm_type_(norm),
      pol_(pol),
      total_evaluations_(0)
  {
    const std::size_t estimated_nodes = 100;  // Estimate for GM rule
    
    workspace_ = std::make_unique<vector_workspace<Real>>(num_components, estimated_nodes);
    aggregator_ = std::make_unique<vector_result_aggregator<Real>>(num_components, norm);
  }
  
  std::vector<result<Real>> integrate() {
    const std::size_t dim = box_.dimension();
    
    // Initial evaluation on whole box
    auto initial_region = std::make_shared<vector_region<Real>>(dim, num_components_);
    initial_region->a = box_.lower;
    initial_region->b = box_.upper;
    
    // Evaluate initial region
    embedded_pair_result_vector<Real> initial_result;
    bool success = genz_malik_evaluator_vector<Real>::evaluate_embedded_pair_runtime(
      f_, *initial_region, initial_result, num_components_, workspace_.get());
    
    if (!success) {
      // Rules not available for this dimension
      std::vector<result<Real>> results(num_components_);
      for (auto& res : results) {
        res.status = status_code::invalid_input;
        res.error = std::numeric_limits<Real>::max();
      }
      return results;
    }
    
    // Store component estimates and errors
    initial_region->component_estimates_fine = initial_result.estimates_fine;
    initial_region->component_estimates_coarse = initial_result.estimates_coarse;
    initial_region->component_errors = initial_result.embedded_errors;
    initial_region->compute_aggregated_error(norm_type_);
    
    // Add to aggregator
    aggregator_->add_results(initial_result.estimates_fine, initial_result.embedded_errors);
    total_evaluations_ += initial_result.evaluations;
    
    // Check if we're done
    Real total_error = aggregator_->get_total_error();
    std::vector<Real> values = aggregator_->get_values();
    Real max_value = *std::max_element(values.begin(), values.end(),
                                       [](Real a, Real b) { return std::abs(a) < std::abs(b); });
    
    if (total_error <= abs_tol_ || total_error <= rel_tol_ * std::abs(max_value)) {
      return create_results(status_code::success);
    }
    
    // Add to queue for refinement
    region_queue_.push(initial_region);
    metrics_.add_iteration(total_error, total_evaluations_);
    
    // Main adaptive loop
    while (!region_queue_.empty() && total_evaluations_ < max_eval_) {
      // Get region with largest error
      auto parent = region_queue_.top();
      region_queue_.pop();
      
      // Remove parent's contribution
      aggregator_->subtract_results(parent->component_estimates_fine, parent->component_errors);
      
      // Split along best dimension
      std::size_t split_dim = parent->split_dimension;
      Real split_point = (parent->a[split_dim] + parent->b[split_dim]) / 2;
      
      // Create child regions
      auto left_child = std::make_shared<vector_region<Real>>(*parent);
      auto right_child = std::make_shared<vector_region<Real>>(*parent);
      
      left_child->b[split_dim] = split_point;
      right_child->a[split_dim] = split_point;
      
      // Evaluate children
      embedded_pair_result_vector<Real> left_result, right_result;
      
      genz_malik_evaluator_vector<Real>::evaluate_embedded_pair_runtime(
        f_, *left_child, left_result, num_components_, workspace_.get());
      
      genz_malik_evaluator_vector<Real>::evaluate_embedded_pair_runtime(
        f_, *right_child, right_result, num_components_, workspace_.get());
      
      // Store results in children
      left_child->component_estimates_fine = left_result.estimates_fine;
      left_child->component_estimates_coarse = left_result.estimates_coarse;
      left_child->component_errors = left_result.embedded_errors;
      left_child->split_dimension = left_result.split_dimension;
      left_child->compute_aggregated_error(norm_type_);
      
      right_child->component_estimates_fine = right_result.estimates_fine;
      right_child->component_estimates_coarse = right_result.estimates_coarse;
      right_child->component_errors = right_result.embedded_errors;
      right_child->split_dimension = right_result.split_dimension;
      right_child->compute_aggregated_error(norm_type_);
      
      // Add children's contributions
      aggregator_->add_results(left_result.estimates_fine, left_result.embedded_errors);
      aggregator_->add_results(right_result.estimates_fine, right_result.embedded_errors);
      
      total_evaluations_ += left_result.evaluations + right_result.evaluations;
      
      // Check convergence
      total_error = aggregator_->get_total_error();
      values = aggregator_->get_values();
      max_value = *std::max_element(values.begin(), values.end(),
                                    [](Real a, Real b) { return std::abs(a) < std::abs(b); });
      
      metrics_.add_iteration(total_error, total_evaluations_);
      
      if (total_error <= abs_tol_ || total_error <= rel_tol_ * std::abs(max_value)) {
        return create_results(status_code::success);
      }
      
      // Add children to queue
      if (left_child->aggregated_error > 0) {
        region_queue_.push(left_child);
      }
      if (right_child->aggregated_error > 0) {
        region_queue_.push(right_child);
      }
    }
    
    // Reached max evaluations
    return create_results(total_evaluations_ >= max_eval_ ? 
                         status_code::maxeval_reached : status_code::success);
  }
  
private:
  std::vector<result<Real>> create_results(status_code status) {
    std::vector<result<Real>> results(num_components_);
    std::vector<Real> values = aggregator_->get_values();
    std::vector<Real> errors = aggregator_->get_errors();
    
    for (std::size_t i = 0; i < num_components_; ++i) {
      results[i].value = values[i];
      results[i].error = std::abs(errors[i]);
      results[i].evaluations = total_evaluations_;
      results[i].status = status;
    }
    
    return results;
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_ADAPTIVE_INTEGRATOR_VECTOR_HPP