// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_WORKSPACE_HPP
#define BOOST_MATH_CUBATURE_WORKSPACE_HPP

#include <boost/math/cubature/policies.hpp>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace boost { namespace math { namespace cubature {

struct execution_options {
  std::size_t max_threads = 0; // 0 -> use hardware_concurrency at call-site
  std::uint64_t max_eval  = 0;
  bool deterministic      = true;
};

/// \brief Thread-local workspace for integration algorithms
/// \details Pre-allocated buffers to avoid dynamic allocation in hot paths
template <typename Real>
class workspace {
private:
  // Buffer sizes (can be configured via constructor)
  static constexpr std::size_t default_node_capacity = 1024;
  static constexpr std::size_t default_value_capacity = 1024;
  static constexpr std::size_t default_region_capacity = 256;
  
public:
  // Node and weight buffers for quadrature rules
  std::vector<Real> abscissas;
  std::vector<Real> weights;
  
  // Function evaluation buffers
  std::vector<Real> values;
  std::vector<Real> scratch;
  
  // Region subdivision buffers (for adaptive algorithms)
  std::vector<Real> region_centers;
  std::vector<Real> region_widths;
  std::vector<Real> region_errors;
  
  // QMC-specific state
  struct qmc_state {
    std::vector<Real> sobol_points;
    std::vector<std::uint32_t> scramble_seeds;
    std::size_t sequence_index;
    std::uint64_t random_state;
    
    qmc_state() : sequence_index(0), random_state(12345) {}
  } qmc;
  
  // Priority queue workspace for adaptive integration
  struct pq_workspace {
    std::vector<std::size_t> indices;
    std::vector<Real> keys;
    std::size_t capacity;
    
    pq_workspace() : capacity(default_region_capacity) {
      indices.reserve(capacity);
      keys.reserve(capacity);
    }
  } priority_queue;
  
  // Sparse grid workspace
  struct sparse_workspace {
    std::vector<std::size_t> multi_indices;
    std::vector<int> coefficients;
    std::unordered_map<std::size_t, Real> node_weights;
    
    void clear() {
      multi_indices.clear();
      coefficients.clear();
      node_weights.clear();
    }
  } sparse;
  
  /// \brief Default constructor with default buffer sizes
  workspace() {
    reserve_buffers(default_node_capacity, default_value_capacity);
  }
  
  /// \brief Constructor with custom buffer sizes
  workspace(std::size_t node_capacity, std::size_t value_capacity) {
    reserve_buffers(node_capacity, value_capacity);
  }
  
  /// \brief Reserve buffer capacity
  void reserve_buffers(std::size_t node_capacity, std::size_t value_capacity) {
    abscissas.reserve(node_capacity);
    weights.reserve(node_capacity);
    values.reserve(value_capacity);
    scratch.reserve(value_capacity);
    
    region_centers.reserve(default_region_capacity * 10);  // Dimension * regions
    region_widths.reserve(default_region_capacity * 10);
    region_errors.reserve(default_region_capacity);
    
    qmc.sobol_points.reserve(node_capacity);
    qmc.scramble_seeds.reserve(64);  // Maximum dimensions
  }
  
  /// \brief Clear all buffers (keeps capacity)
  void clear() {
    abscissas.clear();
    weights.clear();
    values.clear();
    scratch.clear();
    region_centers.clear();
    region_widths.clear();
    region_errors.clear();
    qmc.sobol_points.clear();
    priority_queue.indices.clear();
    priority_queue.keys.clear();
    sparse.clear();
  }
  
  /// \brief Reset QMC state for new integration
  void reset_qmc() {
    qmc.sequence_index = 0;
    qmc.random_state = 12345;
    qmc.sobol_points.clear();
  }
};

// Type-erased workspace for policy-based dispatch
using default_workspace = workspace<double>;

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_WORKSPACE_HPP

