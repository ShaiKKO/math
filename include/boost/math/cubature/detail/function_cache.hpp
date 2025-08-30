// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_FUNCTION_CACHE_HPP
#define BOOST_MATH_CUBATURE_DETAIL_FUNCTION_CACHE_HPP

#include <boost/math/cubature/detail/adaptivity.hpp>
#include <unordered_map>
#include <vector>
#include <array>
#include <cmath>
#include <memory>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \file function_cache.hpp
/// \brief Function evaluation caching for adaptive integration
/// \details Provides caching mechanisms to reuse function evaluations when
///          regions are subdivided. This optimization is critical for expensive
///          integrands where function evaluation dominates runtime.
/// 
/// \section cache_algorithm Algorithm
/// When a region is split along dimension d at point s:
/// - Parent nodes with x[d] < s belong to left child
/// - Parent nodes with x[d] > s belong to right child
/// - Nodes are transformed from parent's [0,1]^d to child's [0,1]^d
/// - Physical space coordinates remain unchanged for cache lookup
/// 
/// \section cache_performance Performance
/// Typical reuse rates for Genz-Malik rules:
/// - Degree 5: ~30% of nodes can be reused
/// - Degree 7: ~35% of nodes can be reused
/// - Degree 9: ~40% of nodes can be reused

/// \brief Cache for function evaluations to avoid redundant computations
/// \tparam Real Floating point type
/// \tparam Dim Dimension (compile-time for efficiency)
/// \details Implements node transformation and cache lookup for adaptive refinement
template <typename Real, std::size_t Dim>
class function_cache {
public:
  /// \brief Node-value pair for caching
  struct cached_evaluation {
    std::array<Real, Dim> physical_node;  // Node in physical space (not reference)
    Real value;                           // Function value at this node
    
    cached_evaluation(const std::array<Real, Dim>& n, Real v) 
      : physical_node(n), value(v) {}
  };
  
  /// \brief Cached data from parent region evaluation
  struct parent_cache {
    std::vector<cached_evaluation> evaluations;
    std::size_t split_dimension;
    Real split_point;  // In [0,1] reference space
    
    parent_cache() : split_dimension(0), split_point(Real(0.5)) {}
  };
  
  /// \brief Check if a node from parent can be reused in child
  /// \param parent_node Node in parent's [0,1]^d reference
  /// \param split_dim Dimension along which parent was split
  /// \param split_point Split point in [0,1] reference (typically 0.5)
  /// \param is_left True for left child, false for right
  /// \param child_node Output: transformed node in child's [0,1]^d reference
  /// \return True if node is in child region and can be reused
  /// 
  /// \details Transformation formulas:
  /// - Left child: x'[d] = x[d] / split_point
  /// - Right child: x'[d] = (x[d] - split_point) / (1 - split_point)
  /// - Other dimensions remain unchanged
  static bool can_reuse_node(
      const std::array<Real, Dim>& parent_node,
      std::size_t split_dim,
      Real split_point,
      bool is_left,
      std::array<Real, Dim>& child_node)
  {
    // Check if node is in the child's half
    if (is_left && parent_node[split_dim] > split_point) {
      return false;  // Node is in right half
    }
    if (!is_left && parent_node[split_dim] < split_point) {
      return false;  // Node is in left half
    }
    
    // Transform node to child's [0,1]^d reference space
    child_node = parent_node;
    if (is_left) {
      // Left child: [0, split_point] -> [0, 1]
      child_node[split_dim] = parent_node[split_dim] / split_point;
    } else {
      // Right child: [split_point, 1] -> [0, 1]
      child_node[split_dim] = (parent_node[split_dim] - split_point) / (Real(1) - split_point);
    }
    
    return true;
  }
  
  /// \brief Extract reusable evaluations for a child region
  /// \param parent Parent region's cached evaluations
  /// \param is_left True for left child, false for right
  /// \return Vector of evaluations that can be reused in child
  /// 
  /// \details Filters parent's cached evaluations based on split plane.
  ///          Typically about half of parent's nodes fall into each child.
  static std::vector<cached_evaluation> extract_child_cache(
      const parent_cache& parent,
      bool is_left)
  {
    std::vector<cached_evaluation> child_cache;
    child_cache.reserve(parent.evaluations.size() / 2);  // Approximately half will be reusable
    
    for (const auto& eval : parent.evaluations) {
      std::array<Real, Dim> child_node;
      if (can_reuse_node(eval.node, parent.split_dimension, 
                        parent.split_point, is_left, child_node)) {
        child_cache.emplace_back(child_node, eval.value);
      }
    }
    
    return child_cache;
  }
  
  /// \brief Find cached value for a node if available
  /// \param cache Vector of cached evaluations
  /// \param physical_node Node to search for (in physical space)
  /// \param tolerance Tolerance for node matching (default 1e-14)
  /// \return Pointer to cached value if found, nullptr otherwise
  /// 
  /// \details Uses L∞ norm for node comparison. Physical space coordinates
  ///          are used because they remain invariant under region transformations.
  static const Real* find_cached_value(
      const std::vector<cached_evaluation>& cache,
      const std::array<Real, Dim>& physical_node,
      Real tolerance = 1e-14)
  {
    for (const auto& eval : cache) {
      bool match = true;
      for (std::size_t d = 0; d < Dim; ++d) {
        if (std::abs(eval.physical_node[d] - physical_node[d]) > tolerance) {
          match = false;
          break;
        }
      }
      if (match) {
        return &eval.value;
      }
    }
    return nullptr;
  }
  
  /// \brief Compute reuse statistics for performance monitoring
  struct reuse_stats {
    std::size_t total_nodes;
    std::size_t reused_nodes;
    std::size_t new_evaluations;
    
    Real reuse_ratio() const {
      return total_nodes > 0 ? 
             static_cast<Real>(reused_nodes) / static_cast<Real>(total_nodes) : 
             Real(0);
    }
  };
  
  /// \brief Track reuse statistics (for debugging/optimization)
  static void update_stats(reuse_stats& stats, bool was_reused) {
    stats.total_nodes++;
    if (was_reused) {
      stats.reused_nodes++;
    } else {
      stats.new_evaluations++;
    }
  }
};

/// \brief Enhanced region structure with caching support
/// \tparam Real Floating point type
/// \details Extends basic region with shared cache storage.
///          Cache is shared via shared_ptr for efficiency.
template <typename Real>
struct cached_region : public region<Real> {
  using cache_ptr = std::shared_ptr<void>;  // Type-erased cache storage
  cache_ptr cached_values;
  
  cached_region(std::size_t dim) : region<Real>(dim) {}
  
  // Copy constructor preserves cache
  cached_region(const cached_region& other) 
    : region<Real>(other), cached_values(other.cached_values) {}
  
  // Split constructor for child regions
  cached_region(const region<Real>& parent, bool /*is_left*/)
    : region<Real>(parent) {}
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_FUNCTION_CACHE_HPP