// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_VECTOR_ADAPTER_HPP
#define BOOST_MATH_CUBATURE_DETAIL_VECTOR_ADAPTER_HPP

#include <boost/math/cubature/detail/adaptivity.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Error norm types for vector-valued integrands
enum class error_norm {
  l_infinity,  // Maximum component error (default, conservative)
  l2,         // Euclidean norm of errors
  l1          // Sum of absolute errors
};

/// \brief Vector result aggregator with configurable error norms
template <typename Real>
class vector_result_aggregator {
private:
  std::size_t num_components_;
  error_norm norm_type_;
  std::vector<kahan_accumulator<Real>> value_accumulators_;
  std::vector<kahan_accumulator<Real>> error_accumulators_;
  
public:
  vector_result_aggregator(std::size_t num_components, 
                          error_norm norm = error_norm::l_infinity)
    : num_components_(num_components),
      norm_type_(norm),
      value_accumulators_(num_components),
      error_accumulators_(num_components) {}
  
  /// \brief Add component results
  void add_results(const std::vector<Real>& values, 
                   const std::vector<Real>& errors) {
    for (std::size_t i = 0; i < num_components_; ++i) {
      value_accumulators_[i].add(values[i]);
      error_accumulators_[i].add(errors[i]);
    }
  }
  
  /// \brief Remove component results (for region refinement)
  void subtract_results(const std::vector<Real>& values,
                       const std::vector<Real>& errors) {
    for (std::size_t i = 0; i < num_components_; ++i) {
      value_accumulators_[i].add(-values[i]);
      error_accumulators_[i].add(-errors[i]);
    }
  }
  
  /// \brief Get aggregated error using configured norm
  Real get_total_error() const {
    switch (norm_type_) {
      case error_norm::l_infinity: {
        Real max_err = Real(0);
        for (const auto& acc : error_accumulators_) {
          max_err = std::max(max_err, std::abs(acc.sum()));
        }
        return max_err;
      }
      
      case error_norm::l2: {
        Real sum_sq = Real(0);
        for (const auto& acc : error_accumulators_) {
          Real err = acc.sum();
          sum_sq += err * err;
        }
        return std::sqrt(sum_sq);
      }
      
      case error_norm::l1: {
        Real sum_abs = Real(0);
        for (const auto& acc : error_accumulators_) {
          sum_abs += std::abs(acc.sum());
        }
        return sum_abs;
      }
      
      default:
        return Real(0);
    }
  }
  
  /// \brief Get component values
  std::vector<Real> get_values() const {
    std::vector<Real> values(num_components_);
    for (std::size_t i = 0; i < num_components_; ++i) {
      values[i] = value_accumulators_[i].sum();
    }
    return values;
  }
  
  /// \brief Get component errors
  std::vector<Real> get_errors() const {
    std::vector<Real> errors(num_components_);
    for (std::size_t i = 0; i < num_components_; ++i) {
      errors[i] = error_accumulators_[i].sum();
    }
    return errors;
  }
};

/// \brief Zero-allocation vector workspace for efficient evaluation
template <typename Real>
class vector_workspace {
private:
  std::vector<Real> values_buffer_;
  std::vector<Real> scratch_buffer_;
  std::size_t num_components_;
  std::size_t max_nodes_;
  
public:
  vector_workspace(std::size_t num_components, std::size_t max_nodes)
    : num_components_(num_components),
      max_nodes_(max_nodes) {
    // Pre-allocate buffers
    values_buffer_.reserve(num_components * max_nodes);
    scratch_buffer_.reserve(num_components);
  }
  
  /// \brief Get temporary buffer for function values
  Real* get_values_buffer(std::size_t node_index) {
    return values_buffer_.data() + node_index * num_components_;
  }
  
  /// \brief Get scratch buffer for single evaluation
  Real* get_scratch_buffer() {
    return scratch_buffer_.data();
  }
  
  /// \brief Reset workspace for new evaluation
  void reset() {
    // No allocation, just reset usage counters if needed
  }
};

/// \brief Adapter to convert scalar integrand to vector interface
template <typename Real, typename ScalarF>
class scalar_to_vector_adapter {
private:
  const ScalarF& f_;
  
public:
  explicit scalar_to_vector_adapter(const ScalarF& f) : f_(f) {}
  
  void operator()(const Real* x, Real* out, std::size_t /*m*/) const {
    out[0] = f_(x);
  }
};

/// \brief Adapter to extract single component from vector integrand
template <typename Real, typename VectorF>
class vector_component_adapter {
private:
  const VectorF& f_;
  std::size_t component_;
  mutable std::vector<Real> buffer_;
  
public:
  vector_component_adapter(const VectorF& f, std::size_t component, std::size_t num_components)
    : f_(f), component_(component), buffer_(num_components) {}
  
  Real operator()(const Real* x) const {
    f_(x, buffer_.data(), buffer_.size());
    return buffer_[component_];
  }
};

/// \brief SIMD-friendly evaluation pattern for vector integrands
template <typename Real>
class vectorized_evaluator {
public:
  /// \brief Evaluate vector integrand at multiple points efficiently
  /// \details Uses cache-friendly access patterns and allows compiler vectorization
  template <typename F, std::size_t Dim>
  static void evaluate_batch(
      const F& f,
      const std::vector<std::array<Real, Dim>>& points,
      std::size_t num_components,
      std::vector<std::vector<Real>>& values_per_component)
  {
    const std::size_t num_points = points.size();
    values_per_component.resize(num_components);
    
    // Pre-allocate component vectors
    for (auto& component_values : values_per_component) {
      component_values.resize(num_points);
    }
    
    // Temporary buffer for single evaluation
    std::vector<Real> temp_values(num_components);
    
    // Evaluate at each point
    for (std::size_t p = 0; p < num_points; ++p) {
      f(points[p].data(), temp_values.data(), num_components);
      
      // Transpose: store by component for better cache locality
      for (std::size_t c = 0; c < num_components; ++c) {
        values_per_component[c][p] = temp_values[c];
      }
    }
  }
  
  /// \brief Compute weighted sums efficiently with potential vectorization
  static Real compute_weighted_sum(
      const std::vector<Real>& values,
      const std::vector<Real>& weights)
  {
    kahan_accumulator<Real> sum;
    
    // This loop pattern enables auto-vectorization
    const std::size_t n = std::min(values.size(), weights.size());
    for (std::size_t i = 0; i < n; ++i) {
      sum.add(values[i] * weights[i]);
    }
    
    return sum.sum();
  }
};

/// \brief Vector region for adaptive refinement
template <typename Real>
struct vector_region : public region<Real> {
  std::vector<Real> component_estimates_fine;
  std::vector<Real> component_estimates_coarse;
  std::vector<Real> component_errors;
  Real aggregated_error;  // Error using chosen norm
  
  vector_region(std::size_t dim, std::size_t num_components) 
    : region<Real>(dim),
      component_estimates_fine(num_components),
      component_estimates_coarse(num_components),
      component_errors(num_components),
      aggregated_error(0) {}
  
  /// \brief Compute aggregated error using specified norm
  void compute_aggregated_error(error_norm norm) {
    switch (norm) {
      case error_norm::l_infinity:
        aggregated_error = *std::max_element(component_errors.begin(), 
                                            component_errors.end());
        break;
        
      case error_norm::l2: {
        Real sum_sq = std::inner_product(component_errors.begin(),
                                        component_errors.end(),
                                        component_errors.begin(),
                                        Real(0));
        aggregated_error = std::sqrt(sum_sq);
        break;
      }
      
      case error_norm::l1:
        aggregated_error = std::accumulate(component_errors.begin(),
                                          component_errors.end(),
                                          Real(0));
        break;
    }
    
    // Store in base class error field for priority queue
    this->error = aggregated_error;
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_VECTOR_ADAPTER_HPP