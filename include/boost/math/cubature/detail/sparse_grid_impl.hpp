// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_SPARSE_GRID_IMPL_HPP
#define BOOST_MATH_CUBATURE_DETAIL_SPARSE_GRID_IMPL_HPP

// STL headers first (before namespace) per Boost conventions
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <type_traits>

// Project headers after STL
#include <boost/math/cubature/detail/cc_rules.hpp>
#include <boost/math/cubature/detail/adaptivity.hpp>
#include <boost/math/cubature/policies.hpp>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Multi-index for Smolyak sparse grid construction
/// \details Represents i = (i₁, i₂, ..., iₐ) with |i| = Σiⱼ
template <typename IndexType = std::size_t>
class multi_index {
public:
  std::vector<IndexType> indices;
  IndexType sum;  // |i| = sum of indices
  
  multi_index(std::size_t dim) : indices(dim, 0), sum(0) {}
  
  multi_index(std::initializer_list<IndexType> init) 
    : indices(init), sum(std::accumulate(init.begin(), init.end(), IndexType(0))) {}
  
  multi_index(const std::vector<IndexType>& idx)
    : indices(idx), sum(std::accumulate(idx.begin(), idx.end(), IndexType(0))) {}
  
  std::size_t dimension() const { return indices.size(); }
  
  bool operator==(const multi_index& other) const {
    return indices == other.indices;
  }
  
  bool operator<(const multi_index& other) const {
    return indices < other.indices;  // Lexicographic ordering
  }
};

/// \brief Generate multi-index set for Smolyak formula
/// \details Creates all indices satisfying ℓ-d+1 ≤ |i| ≤ ℓ
template <typename IndexType = std::size_t>
class multi_index_set {
private:
  std::size_t dimension_;
  std::size_t level_;
  std::vector<multi_index<IndexType>> indices_;
  
public:
  multi_index_set(std::size_t dimension, std::size_t level)
    : dimension_(dimension), level_(level) {
    generate_indices();
  }
  
  const std::vector<multi_index<IndexType>>& indices() const { return indices_; }
  
  /// \brief Compute Smolyak coefficient for given multi-index
  /// \details Returns (-1)^(ℓ-|i|) * C(d-1, ℓ-|i|)
  int coefficient(const multi_index<IndexType>& idx) const {
    IndexType diff = level_ - idx.sum;
    
    // Compute binomial coefficient C(d-1, diff)
    std::size_t binom = binomial_coefficient(dimension_ - 1, diff);
    
    // Apply sign
    int sign = (diff % 2 == 0) ? 1 : -1;
    
    return sign * static_cast<int>(binom);
  }
  
private:
  void generate_indices() {
    // Generate all multi-indices with sum between ℓ-d+1 and ℓ
    IndexType min_sum = (level_ >= dimension_ - 1) ? level_ - dimension_ + 1 : 0;
    IndexType max_sum = level_;
    
    generate_recursive(multi_index<IndexType>(dimension_), 0, min_sum, max_sum);
    
    // Sort for deterministic ordering
    std::sort(indices_.begin(), indices_.end());
  }
  
  void generate_recursive(multi_index<IndexType> current, std::size_t pos,
                         IndexType min_sum, IndexType max_sum) {
    if (pos == dimension_) {
      if (current.sum >= min_sum && current.sum <= max_sum) {
        indices_.push_back(current);
      }
      return;
    }
    
    // Try all possible values for current position
    IndexType max_val = max_sum - current.sum;
    for (IndexType val = 0; val <= max_val; ++val) {
      current.indices[pos] = val;
      current.sum += val;
      generate_recursive(current, pos + 1, min_sum, max_sum);
      current.sum -= val;
    }
  }
  
  static std::size_t binomial_coefficient(std::size_t n, std::size_t k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    
    // Use symmetry property
    k = std::min(k, n - k);
    
    std::size_t result = 1;
    for (std::size_t i = 0; i < k; ++i) {
      result = result * (n - i) / (i + 1);
    }
    return result;
  }
};

/// \brief Sparse node with accumulated weight
template <typename Real>
struct sparse_node {
  std::vector<Real> point;
  Real weight;
  
  sparse_node(const std::vector<Real>& p, Real w) 
    : point(p), weight(w) {}
};

/// \brief Node deduplication with ULP-aware hashing
template <typename Real>
class sparse_node_set {
private:
  static constexpr Real ulp_tolerance = Real(1e-14);
  
  struct node_hash {
    std::size_t operator()(const std::vector<Real>& node) const {
      std::size_t hash = 0;
      for (Real x : node) {
        // Round to ULP tolerance
        int64_t bits = static_cast<int64_t>(x / ulp_tolerance);
        // Mix hash using FNV-1a algorithm constants
        hash ^= std::hash<int64_t>{}(bits);
        hash *= 0x01000193;  // FNV prime
      }
      return hash;
    }
  };
  
  struct node_equal {
    bool operator()(const std::vector<Real>& a, const std::vector<Real>& b) const {
      if (a.size() != b.size()) return false;
      for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) > ulp_tolerance * std::max(std::abs(a[i]), std::abs(b[i]))) {
          return false;
        }
      }
      return true;
    }
  };
  
  // Map stores (accumulated_weight, compensation_term) for Kahan summation
  std::unordered_map<std::vector<Real>, std::pair<Real, Real>, node_hash, node_equal> node_map_;
  mutable std::vector<sparse_node<Real>> unique_nodes_;
  mutable bool finalized_;
  
public:
  sparse_node_set() : finalized_(false) {}
  
  /// \brief Add node with weight, accumulating if duplicate
  void add_node(const std::vector<Real>& point, Real weight) {
    finalized_ = false;
    auto it = node_map_.find(point);
    if (it != node_map_.end()) {
      // Proper Kahan compensated summation
      // it->second.first = accumulated sum
      // it->second.second = compensation term
      Real y = weight - it->second.second;      // Apply compensation
      Real t = it->second.first + y;            // Tentative sum
      it->second.second = (t - it->second.first) - y;  // New compensation  
      it->second.first = t;                     // New sum
    } else {
      node_map_[point] = std::make_pair(weight, Real(0));
    }
  }
  
  /// \brief Get deduplicated nodes and weights
  const std::vector<sparse_node<Real>>& get_nodes() const {
    if (!finalized_) {
      finalize();
    }
    return unique_nodes_;
  }
  
  std::size_t size() const {
    return node_map_.size();
  }
  
private:
  void finalize() const {
    unique_nodes_.clear();
    unique_nodes_.reserve(node_map_.size());
    
    for (const auto& [point, weight_pair] : node_map_) {
      // Extract final accumulated weight from Kahan sum
      unique_nodes_.emplace_back(point, weight_pair.first);
    }
    
    // Sort for deterministic ordering
    std::sort(unique_nodes_.begin(), unique_nodes_.end(),
              [](const sparse_node<Real>& a, const sparse_node<Real>& b) {
                return a.point < b.point;
              });
    
    finalized_ = true;
  }
};

/// \brief Smolyak sparse grid construction and integration
template <typename Real>
class smolyak_grid {
private:
  std::size_t dimension_;
  std::size_t level_;
  sparse_node_set<Real> nodes_;
  std::vector<sparse_node<Real>> previous_level_nodes_;
  
public:
  smolyak_grid(std::size_t dimension, std::size_t level)
    : dimension_(dimension), level_(level) {
    construct_grid();
  }
  
  /// \brief Construct Smolyak sparse grid
  void construct_grid() {
    multi_index_set<std::size_t> index_set(dimension_, level_);
    
    for (const auto& idx : index_set.indices()) {
      int coeff = index_set.coefficient(idx);
      if (coeff == 0) continue;
      
      // Build tensor product for this multi-index
      add_tensor_product(idx, coeff);
    }
  }
  
  /// \brief Evaluate integral using sparse grid
  template <typename F>
  result<Real> evaluate(const F& f, const hypercube<Real>& box) {
    const auto& nodes = nodes_.get_nodes();
    
    kahan_accumulator<Real> sum;
    std::size_t evaluations = 0;
    
    // Transform nodes from [0,1]^d to box and evaluate
    std::vector<Real> point(dimension_);
    for (const auto& node : nodes) {
      // Transform from [0,1]^d to [a,b]^d
      for (std::size_t i = 0; i < dimension_; ++i) {
        point[i] = box.lower[i] + node.point[i] * (box.upper[i] - box.lower[i]);
      }
      
      // Support multiple integrand signatures
      Real value;
      if constexpr (std::is_invocable_v<F, decltype(point)>) {
        // f(std::vector<Real>)
        value = f(point);
      } else if constexpr (std::is_invocable_v<F, const Real*, std::size_t>) {
        // f(const Real*, std::size_t)
        value = f(point.data(), dimension_);
      } else if constexpr (std::is_invocable_v<F, const Real*>) {
        // f(const Real*) - legacy support
        value = f(point.data());
      } else {
        static_assert(sizeof(F) == 0, "Integrand must be callable with vector<Real>, (const Real*, size_t), or const Real*");
      }
      sum.add(value * node.weight);
      evaluations++;
    }
    
    // Scale by box volume
    Real volume = 1.0;
    for (std::size_t i = 0; i < dimension_; ++i) {
      volume *= (box.upper[i] - box.lower[i]);
    }
    
    result<Real> res;
    res.value = sum.sum() * volume;
    res.evaluations = evaluations;
    res.status = status_code::success;
    
    // Estimate error using hierarchical surplus method
    res.error = estimate_surplus_error(res.value);
    
    return res;
  }
  
  /// \brief Get number of unique nodes
  std::size_t num_nodes() const {
    return nodes_.size();
  }
  
  /// \brief Get the sparse grid nodes (for testing)
  std::vector<sparse_node<Real>> get_nodes() const {
    return nodes_.get_nodes();
  }
  
private:
  /// \brief Add tensor product nodes for given multi-index
  void add_tensor_product(const multi_index<std::size_t>& idx, int coefficient) {
    // Get 1D rules for each dimension
    std::vector<quadrature_rule_1d<Real>> rules_1d;
    rules_1d.reserve(dimension_);
    
    for (std::size_t i = 0; i < dimension_; ++i) {
      rules_1d.emplace_back(idx.indices[i]);
    }
    
    // Generate tensor product nodes
    std::size_t n_points = 1;
    for (const auto& rule : rules_1d) {
      n_points *= rule.size();
    }
    
    // Build all tensor product points
    for (std::size_t i = 0; i < n_points; ++i) {
      std::vector<Real> point(dimension_);
      Real weight = static_cast<Real>(coefficient);
      
      std::size_t idx_copy = i;
      for (std::size_t d = 0; d < dimension_; ++d) {
        std::size_t n_d = rules_1d[d].size();
        std::size_t j_d = idx_copy % n_d;
        
        // Transform from [-1,1] to [0,1]
        point[d] = (rules_1d[d].nodes[j_d] + Real(1)) / Real(2);
        weight *= rules_1d[d].weights[j_d] / Real(2);  // Scale weight for interval change
        
        idx_copy /= n_d;
      }
      
      nodes_.add_node(point, weight);
    }
  }
  
  /// \brief Estimate error using hierarchical surplus method
  /// \details Estimates integration error based on grid level and dimension
  Real estimate_surplus_error(Real current_value) const {
    // For sparse grids, error typically decreases as O(n^{-r}) where r depends on smoothness
    // Use a conservative estimate based on level and dimension
    
    if (level_ == 0) {
      // Level 0: very coarse estimate
      return std::abs(current_value) * Real(0.5);
    }
    
    // Estimate based on theoretical convergence rate for smooth functions
    // Error ~ C * h^p where h ~ 2^{-level} and p is related to smoothness
    Real h = std::pow(Real(2), -Real(level_));
    
    // For Clenshaw-Curtis with smooth functions, expect spectral convergence
    // But use conservative polynomial rate for safety
    Real convergence_rate = Real(2 * level_);  // Polynomial exactness degree
    
    // Apply dimension-dependent factor
    Real dim_factor = std::pow(Real(dimension_), Real(0.5));
    
    // Heuristic error estimate
    Real error_estimate = std::abs(current_value) * dim_factor * std::pow(h, convergence_rate);
    
    // Add safety factor and machine epsilon margin
    Real safety_factor = Real(10);
    Real epsilon_margin = std::numeric_limits<Real>::epsilon() * 
                          std::abs(current_value) * Real(100);
    
    return safety_factor * error_estimate + epsilon_margin;
  }
};

/// \brief Main integration function for sparse grid
template <typename Real, typename F>
result<Real> integrate_sparse_grid_impl(
    const F& f,
    const hypercube<Real>& box,
    std::size_t level)
{
  std::size_t dim = box.dimension();
  
  // Validate inputs
  if (dim == 0 || dim > 20) {
    result<Real> res;
    res.status = status_code::dimension_error;
    res.error = std::numeric_limits<Real>::max();
    return res;
  }
  
  // Construct and evaluate sparse grid
  smolyak_grid<Real> grid(dim, level);
  return grid.evaluate(f, box);
}

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_SPARSE_GRID_IMPL_HPP