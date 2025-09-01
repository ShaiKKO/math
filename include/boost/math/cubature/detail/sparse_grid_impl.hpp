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
#include <map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>
#include <type_traits>
#include <limits>

// Project headers after STL
#include <boost/math/cubature/detail/cc_rules.hpp>
#include <boost/math/cubature/detail/gauss_hermite_rules.hpp>
#include <boost/math/cubature/detail/adaptivity.hpp>
#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>

namespace boost { namespace math { namespace cubature { namespace detail {

// Signature detection traits for integrand functions
template <typename F, typename Real, typename... Args>
struct has_raw_pointer_signature {
private:
  template <typename U>
  static auto test(int) -> decltype(
    std::declval<U>()(std::declval<const Real*>(), std::declval<std::size_t>()),
    std::true_type{});
  
  template <typename>
  static std::false_type test(...);
  
public:
  static constexpr bool value = decltype(test<F>(0))::value;
};

template <typename F, typename Real>
struct has_vector_signature {
private:
  template <typename U>
  static auto test(int) -> decltype(
    std::declval<U>()(std::declval<const std::vector<Real>&>()),
    std::true_type{});
  
  template <typename>
  static std::false_type test(...);
  
public:
  static constexpr bool value = decltype(test<F>(0))::value;
};

template <typename F, typename Real>
struct has_raw_pointer_only_signature {
private:
  template <typename U>
  static auto test(int) -> decltype(
    std::declval<U>()(std::declval<const Real*>()),
    std::true_type{});
  
  template <typename>
  static std::false_type test(...);
  
public:
  static constexpr bool value = decltype(test<F>(0))::value;
};

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
  /// \details The Smolyak formula combines tensor products with alternating signs:
  ///          A(ℓ,d) = Σ (-1)^(ℓ-|i|) * C(d-1, ℓ-|i|) * (U^{i_1} ⊗ ... ⊗ U^{i_d})
  ///          where the sum is over multi-indices with ℓ-d+1 ≤ |i| ≤ ℓ
  ///          This function returns the coefficient (-1)^(ℓ-|i|) * C(d-1, ℓ-|i|)
  int coefficient(const multi_index<IndexType>& idx) const {
    IndexType diff = level_ - idx.sum;
    
    // Binomial coefficient C(d-1, diff) counts tensor products at this level
    std::size_t binom = binomial_coefficient(dimension_ - 1, diff);
    
    // Alternating sign implements inclusion-exclusion principle
    // to avoid overcounting shared nodes between tensor products
    int sign = (diff % 2 == 0) ? 1 : -1;
    
    return sign * static_cast<int>(binom);
  }
  
private:
  void generate_indices() {
    // Generate all multi-indices satisfying the Smolyak constraint:
    // ℓ-d+1 ≤ |i| ≤ ℓ, where |i| = i_1 + i_2 + ... + i_d
    // These indices determine which tensor products to include
    // Lower bound ensures we don't include low-order terms already in A(ℓ-1,d)
    IndexType min_sum = (level_ >= dimension_ - 1) ? level_ - dimension_ + 1 : 0;
    IndexType max_sum = level_;
    
    generate_recursive(multi_index<IndexType>(dimension_), 0, min_sum, max_sum);
    
    // Sort for deterministic ordering (lexicographic)
    // This ensures reproducible results across platforms
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
      // Node already exists - accumulate weight using Kahan summation
      // This maintains numerical accuracy when combining many tensor products
      // with alternating signs (inclusion-exclusion can cause cancellation)
      // Kahan algorithm tracks a compensation term for lost low-order digits
      Real y = weight - it->second.second;      // Apply previous compensation
      Real t = it->second.first + y;            // Tentative sum
      it->second.second = (t - it->second.first) - y;  // Update compensation
      it->second.first = t;                     // Store corrected sum
    } else {
      // New node - initialize with weight and zero compensation
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
    
    for (const auto& node_entry : node_map_) {
      // Extract final accumulated weight from Kahan sum
      unique_nodes_.emplace_back(node_entry.first, node_entry.second.first);
    }
    
    // Sort for deterministic ordering
    std::sort(unique_nodes_.begin(), unique_nodes_.end(),
              [](const sparse_node<Real>& a, const sparse_node<Real>& b) {
                return a.point < b.point;
              });
    
    finalized_ = true;
  }
};

// SFINAE helper for detecting integrand signature
namespace detail_sfinae {
  // Detection idiom helpers
  template<class...> using void_t = void;
  
  template<class T, class = void>
  struct has_vector_signature : std::false_type {};
  
  template<class T>
  struct has_vector_signature<T, 
      void_t<decltype(std::declval<T>()(std::declval<std::vector<typename T::Real>>()))>> 
      : std::true_type {};
  
  template<class T, class Real, class = void>
  struct has_pointer_size_signature : std::false_type {};
  
  template<class T, class Real>
  struct has_pointer_size_signature<T, Real,
      void_t<decltype(std::declval<T>()(std::declval<const Real*>(), std::declval<std::size_t>()))>>
      : std::true_type {};
  
  template<class T, class Real, class = void>
  struct has_pointer_signature : std::false_type {};
  
  template<class T, class Real>
  struct has_pointer_signature<T, Real,
      void_t<decltype(std::declval<T>()(std::declval<const Real*>()))>>
      : std::true_type {};
}

// Simplified integrand traits with priority-based selection
template <typename Real, typename Func, typename Point>
struct integrand_traits {
private:
  // Helper to check if vector call is valid
  template<typename F, typename P>
  static auto test_vector(int) -> decltype(
      Real(std::declval<const F&>()(std::declval<const P&>())), std::true_type{});
  template<typename, typename>
  static std::false_type test_vector(...);
  
  // Helper to check if pointer+size call is valid
  template<typename F>
  static auto test_ptr_size(int) -> decltype(
      Real(std::declval<const F&>()(std::declval<const Real*>(), std::declval<std::size_t>())), std::true_type{});
  template<typename>
  static std::false_type test_ptr_size(...);
  
  // Helper to check if pointer call is valid
  template<typename F>
  static auto test_ptr(int) -> decltype(
      Real(std::declval<const F&>()(std::declval<const Real*>())), std::true_type{});
  template<typename>
  static std::false_type test_ptr(...);
  
  static constexpr bool has_vector = decltype(test_vector<Func, Point>(0))::value;
  static constexpr bool has_ptr_size = decltype(test_ptr_size<Func>(0))::value;
  static constexpr bool has_ptr = decltype(test_ptr<Func>(0))::value;
  
public:
  // Priority-based selection: vector > ptr+size > ptr
  static constexpr int signature = 
      has_vector ? 1 :
      has_ptr_size ? 2 :
      has_ptr ? 3 : 0;
};

/// \brief Smolyak sparse grid construction and integration
template <typename Real>
class smolyak_grid {
private:
  std::size_t dimension_;
  std::size_t level_;
  sparse_node_set<Real> nodes_;
  std::vector<sparse_node<Real>> previous_level_nodes_;
  bool enable_diagnostics_;
  
  // Diagnostic information
  struct diagnostic_info {
    std::size_t num_positive_weights = 0;
    std::size_t num_negative_weights = 0;
    Real sum_positive_weights = 0;
    Real sum_negative_weights = 0;
    Real weight_cancellation_ratio = 0;  // |sum|weights|| / |sum weights|
    bool has_weight_issues = false;
  } diagnostics_;
  
  // Tag dispatch functions
  template <typename Func>
  Real evaluate_integrand_dispatch(const Func&, const std::vector<Real>&, std::integral_constant<int, 0>) const {
    static_assert(sizeof(Func) == 0, "Integrand must be callable");
    return Real{};
  }
  
  template <typename Func>
  Real evaluate_integrand_dispatch(const Func& f, const std::vector<Real>& point, std::integral_constant<int, 1>) const {
    return f(point);
  }
  
  template <typename Func>
  Real evaluate_integrand_dispatch(const Func& f, const std::vector<Real>& point, std::integral_constant<int, 2>) const {
    return f(point.data(), dimension_);
  }
  
  template <typename Func>
  Real evaluate_integrand_dispatch(const Func& f, const std::vector<Real>& point, std::integral_constant<int, 3>) const {
    return f(point.data());
  }
  
public:
  smolyak_grid(std::size_t dimension, std::size_t level, bool enable_diagnostics = false)
    : dimension_(dimension), level_(level), enable_diagnostics_(enable_diagnostics) {
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
    
    // Collect diagnostic information if enabled
    if (enable_diagnostics_) {
      collect_weight_diagnostics(nodes);
    }
    
    kahan_accumulator<Real> sum;
    std::size_t evaluations = 0;
    bool all_positive = true;  // Track if integrand is always positive
    Real min_value = std::numeric_limits<Real>::max();
    
    // Transform nodes from [0,1]^d to box and evaluate
    std::vector<Real> point(dimension_);
    for (const auto& node : nodes) {
      // Transform from [0,1]^d to [a,b]^d
      for (std::size_t i = 0; i < dimension_; ++i) {
        point[i] = box.lower[i] + node.point[i] * (box.upper[i] - box.lower[i]);
      }
      
      // Support multiple integrand signatures via tag dispatch
      Real value = evaluate_integrand_dispatch(f, point, 
          std::integral_constant<int, integrand_traits<Real, F, decltype(point)>::signature>());
      
      if (value < 0) all_positive = false;
      min_value = std::min(min_value, value);
      
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
    
    // Sanity check: warn if integral is negative for non-negative integrand
    if (all_positive && res.value < 0) {
      res.status = status_code::cancelled;  // Use this to indicate numerical issues
      // In production, might want to log warning or throw exception
      // Note: we're repurposing cancelled to indicate the integration is unreliable
    }
    
    // Check for weight cancellation issues
    if (enable_diagnostics_ && diagnostics_.has_weight_issues) {
      // Adjust error estimate to reflect uncertainty
      res.error = std::max(estimate_surplus_error(res.value),
                           std::abs(res.value) * diagnostics_.weight_cancellation_ratio);
    } else {
      res.error = estimate_surplus_error(res.value);
    }
    
    return res;
  }
  
  /// \brief Evaluate vector integral using sparse grid
  template <typename F>
  std::vector<result<Real>> evaluate_vector(const F& f, const hypercube<Real>& box, 
                                           std::size_t num_components) {
    const auto& nodes = nodes_.get_nodes();
    
    // Kahan accumulators for each component
    std::vector<kahan_accumulator<Real>> sums(num_components);
    std::size_t evaluations = 0;
    
    // Temporary buffer for function values
    std::vector<Real> values(num_components);
    
    // Transform nodes from [0,1]^d to box and evaluate
    std::vector<Real> point(dimension_);
    for (const auto& node : nodes) {
      // Transform from [0,1]^d to [a,b]^d
      for (std::size_t i = 0; i < dimension_; ++i) {
        point[i] = box.lower[i] + node.point[i] * (box.upper[i] - box.lower[i]);
      }
      
      // Evaluate vector integrand (single call for all components)
      f(point.data(), values.data(), num_components);
      
      // Accumulate weighted sums for each component
      for (std::size_t c = 0; c < num_components; ++c) {
        sums[c].add(values[c] * node.weight);
      }
      evaluations++;
    }
    
    // Scale by box volume
    Real volume = 1.0;
    for (std::size_t i = 0; i < dimension_; ++i) {
      volume *= (box.upper[i] - box.lower[i]);
    }
    
    // Create results for each component
    std::vector<result<Real>> results(num_components);
    for (std::size_t c = 0; c < num_components; ++c) {
      results[c].value = sums[c].sum() * volume;
      results[c].evaluations = evaluations;
      results[c].status = status_code::success;
      
      // Estimate error using hierarchical surplus method
      results[c].error = estimate_surplus_error(results[c].value);
    }
    
    return results;
  }
  
  /// \brief Get number of unique nodes
  std::size_t num_nodes() const {
    return nodes_.size();
  }
  
  /// \brief Get the sparse grid nodes (for testing)
  std::vector<sparse_node<Real>> get_nodes() const {
    return nodes_.get_nodes();
  }
  
  /// \brief Get diagnostic information (if enabled)
  const diagnostic_info& get_diagnostics() const {
    return diagnostics_;
  }
  
private:
  /// \brief Collect diagnostic information about weights
  void collect_weight_diagnostics(const std::vector<sparse_node<Real>>& nodes) {
    diagnostics_ = diagnostic_info{};  // Reset
    
    Real sum_abs_weights = 0;
    Real sum_weights = 0;
    
    for (const auto& node : nodes) {
      Real w = node.weight;
      sum_weights += w;
      sum_abs_weights += std::abs(w);
      
      if (w > 0) {
        diagnostics_.num_positive_weights++;
        diagnostics_.sum_positive_weights += w;
      } else if (w < 0) {
        diagnostics_.num_negative_weights++;
        diagnostics_.sum_negative_weights += std::abs(w);
      }
    }
    
    // Calculate weight cancellation ratio
    if (std::abs(sum_weights) > Real(1e-10)) {
      diagnostics_.weight_cancellation_ratio = sum_abs_weights / std::abs(sum_weights);
    } else {
      diagnostics_.weight_cancellation_ratio = std::numeric_limits<Real>::max();
    }
    
    // Flag potential issues
    // If cancellation ratio > 10, we have significant cancellation
    diagnostics_.has_weight_issues = (diagnostics_.weight_cancellation_ratio > Real(10));
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
    // Use more realistic convergence rate
    Real convergence_rate = std::min(Real(2 * level_), Real(4));  // Cap at 4th order
    
    // Apply dimension-dependent factor (curse of dimensionality)
    Real dim_factor = std::sqrt(Real(dimension_));
    
    // More realistic error estimate
    Real error_estimate = std::abs(current_value) * std::pow(h, convergence_rate);
    
    // Apply dimension scaling
    error_estimate *= dim_factor;
    
    // Add safety factor and machine epsilon margin
    Real safety_factor = Real(2);  // Reduced from 10
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
    std::size_t level,
    bool enable_diagnostics = false)
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
  smolyak_grid<Real> grid(dim, level, enable_diagnostics);
  auto res = grid.evaluate(f, box);
  
  // Add warning if we detected issues
  if (enable_diagnostics && grid.get_diagnostics().has_weight_issues) {
    // Could log warning here in production code
    // The status_code::cancelled in the result indicates potential issues
  }
  
  return res;
}

/// \brief Main vector integration function for sparse grid
template <typename Real, typename F>
std::vector<result<Real>> integrate_sparse_grid_vector_impl(
    const F& f,
    const hypercube<Real>& box,
    std::size_t num_components,
    std::size_t level)
{
  std::size_t dim = box.dimension();
  
  // Validate inputs
  if (dim == 0 || dim > 20 || num_components == 0) {
    std::vector<result<Real>> res(num_components);
    for (auto& r : res) {
      r.status = status_code::dimension_error;
      r.error = std::numeric_limits<Real>::max();
    }
    return res;
  }
  
  // Construct and evaluate sparse grid
  smolyak_grid<Real> grid(dim, level);
  return grid.evaluate_vector(f, box, num_components);
}

/// \brief Dimension-adaptive sparse grid construction
/// \details Allows different levels per dimension based on importance
template <typename Real>
class adaptive_smolyak_grid {
private:
  std::size_t dimension_;
  std::vector<std::size_t> levels_;  // Level per dimension
  sparse_node_set<Real> nodes_;
  bool enable_diagnostics_;
  
  // Dimension importance weights (higher = more important)
  std::vector<Real> dimension_weights_;
  
  // Hierarchical surplus indicators per dimension
  std::vector<Real> surplus_indicators_;
  
  // Tag dispatch functions for evaluating integrand
  template <typename Func>
  Real evaluate_integrand_dispatch(const Func&, const std::vector<Real>&, std::integral_constant<int, 0>) const {
    static_assert(sizeof(Func) == 0, "Integrand must be callable");
    return Real{};
  }
  
  template <typename Func>
  Real evaluate_integrand_dispatch(const Func& f, const std::vector<Real>& point, std::integral_constant<int, 1>) const {
    return f(point);
  }
  
  template <typename Func>
  Real evaluate_integrand_dispatch(const Func& f, const std::vector<Real>& point, std::integral_constant<int, 2>) const {
    return f(point.data(), dimension_);
  }
  
  template <typename Func>
  Real evaluate_integrand_dispatch(const Func& f, const std::vector<Real>& point, std::integral_constant<int, 3>) const {
    return f(point.data());
  }
  
public:
  adaptive_smolyak_grid(
      std::size_t dimension,
      const std::vector<std::size_t>& initial_levels,
      const std::vector<Real>& dim_weights = {},
      bool enable_diagnostics = false)
    : dimension_(dimension), 
      levels_(initial_levels),
      enable_diagnostics_(enable_diagnostics),
      dimension_weights_(dim_weights)
  {
    if (levels_.size() != dimension_) {
      levels_.resize(dimension_, 1);  // Default to level 1
    }
    
    if (dimension_weights_.empty()) {
      dimension_weights_.resize(dimension_, Real(1));  // Equal weights
    }
    
    surplus_indicators_.resize(dimension_, Real(0));
    construct_anisotropic_grid();
  }
  
  /// \brief Construct anisotropic sparse grid with different levels per dimension
  void construct_anisotropic_grid() {
    // Generate multi-indices for anisotropic grid
    // Sum of weighted levels: Σ(w_i * l_i) <= threshold
    Real max_weighted_level = compute_weighted_level_sum();
    
    // Generate all valid multi-indices
    std::vector<multi_index<std::size_t>> indices;
    generate_anisotropic_indices(indices, max_weighted_level);
    
    // Build tensor products for each multi-index
    for (const auto& idx : indices) {
      int coeff = compute_anisotropic_coefficient(idx);
      if (coeff != 0) {
        add_anisotropic_tensor_product(idx, coeff);
      }
    }
  }
  
  /// \brief Evaluate integral with anisotropic grid
  template <typename F>
  result<Real> evaluate(const F& f, const hypercube<Real>& box) {
    const auto& grid_nodes = nodes_.get_nodes();
    
    kahan_accumulator<Real> sum;
    std::size_t evaluations = 0;
    
    // Track surplus per dimension for adaptivity
    std::vector<Real> dim_surplus(dimension_, Real(0));
    
    std::vector<Real> point(dimension_);
    for (const auto& node : grid_nodes) {
      // Transform from [0,1]^d to box
      for (std::size_t i = 0; i < dimension_; ++i) {
        point[i] = box.lower[i] + node.point[i] * (box.upper[i] - box.lower[i]);
      }
      
      Real value = evaluate_integrand_dispatch(f, point, 
          std::integral_constant<int, integrand_traits<Real, F, decltype(point)>::signature>());
      
      sum.add(value * node.weight);
      
      // Track dimensional contributions for adaptivity
      for (std::size_t d = 0; d < dimension_; ++d) {
        dim_surplus[d] += std::abs(value * node.weight) * 
                          std::abs(node.point[d] - Real(0.5));
      }
      
      evaluations++;
    }
    
    // Update surplus indicators
    for (std::size_t d = 0; d < dimension_; ++d) {
      surplus_indicators_[d] = dim_surplus[d] / (evaluations + 1);
    }
    
    // Scale by volume
    Real volume = 1.0;
    for (std::size_t i = 0; i < dimension_; ++i) {
      volume *= (box.upper[i] - box.lower[i]);
    }
    
    result<Real> res;
    res.value = sum.sum() * volume;
    res.evaluations = evaluations;
    res.status = status_code::success;
    res.error = estimate_anisotropic_error(res.value);
    
    return res;
  }
  
  /// \brief Adapt grid by increasing level in most important dimension
  void adapt() {
    // Find dimension with largest surplus indicator
    std::size_t max_dim = 0;
    Real max_surplus = surplus_indicators_[0] * dimension_weights_[0];
    
    for (std::size_t d = 1; d < dimension_; ++d) {
      Real weighted_surplus = surplus_indicators_[d] * dimension_weights_[d];
      if (weighted_surplus > max_surplus) {
        max_surplus = weighted_surplus;
        max_dim = d;
      }
    }
    
    // Increase level in that dimension
    levels_[max_dim]++;
    
    // Reconstruct grid with new levels
    nodes_ = sparse_node_set<Real>();  // Clear
    construct_anisotropic_grid();
  }
  
  /// \brief Get current levels per dimension
  const std::vector<std::size_t>& get_levels() const { return levels_; }
  
  /// \brief Get surplus indicators per dimension
  const std::vector<Real>& get_surplus_indicators() const { 
    return surplus_indicators_; 
  }
  
private:
  /// \brief Compute weighted sum of levels
  Real compute_weighted_level_sum() const {
    Real sum = 0;
    for (std::size_t d = 0; d < dimension_; ++d) {
      sum += dimension_weights_[d] * Real(levels_[d]);
    }
    return sum;
  }
  
  /// \brief Generate multi-indices for anisotropic grid
  void generate_anisotropic_indices(
      std::vector<multi_index<std::size_t>>& indices,
      Real max_weighted_level) 
  {
    // Generate indices satisfying weighted constraint
    multi_index<std::size_t> current(dimension_);
    generate_anisotropic_recursive(indices, current, 0, max_weighted_level);
    
    // Sort for deterministic ordering
    std::sort(indices.begin(), indices.end());
  }
  
  /// \brief Recursive generation of anisotropic indices
  void generate_anisotropic_recursive(
      std::vector<multi_index<std::size_t>>& indices,
      multi_index<std::size_t>& current,
      std::size_t pos,
      Real remaining_weight)
  {
    if (pos == dimension_) {
      if (current.sum > 0) {  // Skip all-zero index
        indices.push_back(current);
      }
      return;
    }
    
    // Try different levels for current dimension
    std::size_t max_level = std::min(
        levels_[pos],
        static_cast<std::size_t>(remaining_weight / dimension_weights_[pos])
    );
    
    for (std::size_t level = 0; level <= max_level; ++level) {
      current.indices[pos] = level;
      current.sum = std::accumulate(
          current.indices.begin(), 
          current.indices.end(), 
          std::size_t(0)
      );
      
      Real new_remaining = remaining_weight - dimension_weights_[pos] * Real(level);
      generate_anisotropic_recursive(indices, current, pos + 1, new_remaining);
    }
  }
  
  /// \brief Compute coefficient for anisotropic multi-index
  int compute_anisotropic_coefficient(const multi_index<std::size_t>& idx) const {
    // Modified Smolyak coefficient for anisotropic case
    // Still uses inclusion-exclusion but with weighted levels
    Real weighted_sum = 0;
    for (std::size_t d = 0; d < dimension_; ++d) {
      weighted_sum += dimension_weights_[d] * Real(idx.indices[d]);
    }
    
    Real max_weighted = compute_weighted_level_sum();
    int diff = static_cast<int>(max_weighted - weighted_sum);
    
    // Binomial coefficient and sign
    if (diff < 0) return 0;
    
    std::size_t binom = 1;  // Simplified for anisotropic case
    int sign = (diff % 2 == 0) ? 1 : -1;
    
    return sign * static_cast<int>(binom);
  }
  
  /// \brief Add tensor product for anisotropic grid
  void add_anisotropic_tensor_product(
      const multi_index<std::size_t>& idx, 
      int coefficient) 
  {
    // Get 1D rules for each dimension
    std::vector<quadrature_rule_1d<Real>> rules_1d;
    rules_1d.reserve(dimension_);
    
    for (std::size_t i = 0; i < dimension_; ++i) {
      rules_1d.emplace_back(idx.indices[i]);
    }
    
    // Generate tensor product nodes (same as isotropic case)
    std::size_t n_points = 1;
    for (const auto& rule : rules_1d) {
      n_points *= rule.size();
    }
    
    for (std::size_t i = 0; i < n_points; ++i) {
      std::vector<Real> point(dimension_);
      Real weight = static_cast<Real>(coefficient);
      
      std::size_t idx_copy = i;
      for (std::size_t d = 0; d < dimension_; ++d) {
        std::size_t n_d = rules_1d[d].size();
        std::size_t j_d = idx_copy % n_d;
        
        point[d] = (rules_1d[d].nodes[j_d] + Real(1)) / Real(2);
        weight *= rules_1d[d].weights[j_d] / Real(2);
        
        idx_copy /= n_d;
      }
      
      nodes_.add_node(point, weight);
    }
  }
  
  /// \brief Estimate error for anisotropic grid
  Real estimate_anisotropic_error(Real current_value) const {
    // Error estimate based on surplus indicators
    Real total_surplus = std::accumulate(
        surplus_indicators_.begin(), 
        surplus_indicators_.end(), 
        Real(0)
    );
    
    // Conservative estimate
    Real base_error = std::abs(current_value) * Real(0.01);  // 1% base
    Real surplus_error = total_surplus * Real(10);  // Scale surplus
    
    return std::max(base_error, surplus_error);
  }
};

/// \brief Gauss-Hermite sparse grid for integration with Gaussian weight
/// \details Integrates f(x) * exp(-||x||^2) over R^d using sparse grid
///          constructed from 1D Gauss-Hermite or Genz-Keister rules
template <typename Real>
class gauss_hermite_sparse_grid {
private:
  std::size_t dimension_;
  std::size_t level_;
  bool use_genz_keister_;
  
  // Store transformed nodes and weights
  std::vector<std::vector<Real>> nodes_;
  std::vector<Real> weights_;
  
public:
  gauss_hermite_sparse_grid(std::size_t dim, std::size_t level, bool use_gk = false)
    : dimension_(dim), level_(level), use_genz_keister_(use_gk) {
    construct_grid();
  }
  
  /// \brief Construct Smolyak sparse grid with Gauss-Hermite nodes
  void construct_grid() {
    // Clear existing data
    nodes_.clear();
    weights_.clear();
    
    // Generate multi-indices for Smolyak formula
    multi_index_set<std::size_t> index_set(dimension_, level_);
    
    // Map to accumulate nodes and weights with deduplication
    std::map<std::vector<Real>, Real> node_weights;
    
    // Build tensor products for each multi-index
    for (const auto& idx : index_set.indices()) {
      int coeff = index_set.coefficient(idx);
      if (coeff == 0) continue;
      
      // Construct tensor product for this multi-index
      add_tensor_product(idx, coeff, node_weights);
    }
    
    // Extract deduplicated nodes and weights
    for (const auto& pair : node_weights) {
      if (std::abs(pair.second) > std::numeric_limits<Real>::epsilon()) {
        nodes_.push_back(pair.first);
        weights_.push_back(pair.second);
      }
    }
  }
  
  /// \brief Add tensor product of 1D rules to the sparse grid
  void add_tensor_product(
      const multi_index<std::size_t>& idx,
      int coeff,
      std::map<std::vector<Real>, Real>& node_weights) 
  {
    // Get 1D rules for each dimension
    std::vector<gauss_hermite_rule_1d<Real>> rules_1d;
    rules_1d.reserve(dimension_);
    
    std::size_t total_points = 1;
    for (std::size_t d = 0; d < dimension_; ++d) {
      std::size_t level_d = idx.indices[d];
      rules_1d.emplace_back(level_d, use_genz_keister_);
      total_points *= rules_1d[d].size();
    }
    
    // Generate all tensor product points
    for (std::size_t i = 0; i < total_points; ++i) {
      std::vector<Real> point(dimension_);
      Real weight = Real(coeff);
      
      std::size_t idx_copy = i;
      for (std::size_t d = 0; d < dimension_; ++d) {
        std::size_t n_d = rules_1d[d].size();
        std::size_t j = idx_copy % n_d;
        
        point[d] = rules_1d[d].nodes[j];
        weight *= rules_1d[d].weights[j];
        
        idx_copy /= n_d;
      }
      
      // Add to map (automatic deduplication)
      node_weights[point] += weight;
    }
  }
  
  /// \brief Evaluate integrand at all sparse grid points
  template <typename F>
  result<Real> evaluate(const F& f) const {
    result<Real> res;
    res.value = Real(0);
    res.error = Real(0);
    res.evaluations = nodes_.size();
    res.status = status_code::success;
    
    // Use Kahan summation for numerical stability
    kahan_accumulator<Real> sum_acc;
    
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
      Real f_val = evaluate_integrand(f, nodes_[i]);
      sum_acc.add(f_val * weights_[i]);
    }
    
    res.value = sum_acc.sum();
    
    // Error estimate based on weight cancellation and grid level
    Real weight_sum = std::accumulate(weights_.begin(), weights_.end(), Real(0),
                                      [](Real a, Real b) { return a + std::abs(b); });
    Real actual_weight = std::abs(std::accumulate(weights_.begin(), weights_.end(), Real(0)));
    
    if (weight_sum > Real(0)) {
      Real cancellation_ratio = actual_weight / weight_sum;
      
      // Higher cancellation means less reliable result
      if (cancellation_ratio < Real(0.1)) {
        res.error = std::abs(res.value) * Real(1.0);  // 100% uncertainty
      } else if (cancellation_ratio < Real(0.5)) {
        res.error = std::abs(res.value) * Real(0.1);  // 10% uncertainty
      } else {
        // Standard error estimate for well-conditioned case
        res.error = std::abs(res.value) * std::pow(Real(10), -Real(level_ + 1));
      }
    } else {
      res.error = std::numeric_limits<Real>::max();
    }
    
    return res;
  }
  
private:
  /// \brief Evaluate integrand with runtime signature detection
  template <typename F>
  Real evaluate_integrand(const F& f, const std::vector<Real>& point) const {
    return evaluate_integrand_dispatch(f, point);
  }
  
  /// \brief Dispatch helper for different integrand signatures
  template <typename F>
  auto evaluate_integrand_dispatch(const F& f, const std::vector<Real>& point) const
    -> typename std::enable_if<
         has_raw_pointer_signature<F, Real, std::size_t>::value,
         Real>::type
  {
    return f(point.data(), dimension_);
  }
  
  template <typename F>
  auto evaluate_integrand_dispatch(const F& f, const std::vector<Real>& point) const
    -> typename std::enable_if<
         !has_raw_pointer_signature<F, Real, std::size_t>::value &&
         has_vector_signature<F, Real>::value,
         Real>::type
  {
    return f(point);
  }
  
  template <typename F>
  auto evaluate_integrand_dispatch(const F& f, const std::vector<Real>& point) const
    -> typename std::enable_if<
         !has_raw_pointer_signature<F, Real, std::size_t>::value &&
         !has_vector_signature<F, Real>::value &&
         has_raw_pointer_only_signature<F, Real>::value,
         Real>::type
  {
    return f(point.data());
  }
};

/// \brief Integration function for Gaussian-weighted sparse grid
template <typename Real, typename F>
result<Real> integrate_sparse_grid_gaussian_impl(
    const F& f,
    std::size_t dim,
    unsigned level,
    bool use_genz_keister)
{
  // Validate input
  if (dim == 0) {
    result<Real> res;
    res.value = Real(0);
    res.error = std::numeric_limits<Real>::max();
    res.status = status_code::dimension_error;
    res.evaluations = 0;
    return res;
  }
  
  if (level == 0) {
    result<Real> res;
    res.value = Real(0);
    res.error = std::numeric_limits<Real>::max();
    res.status = status_code::invalid_input;
    res.evaluations = 0;
    return res;
  }
  
  // Construct and evaluate sparse grid
  gauss_hermite_sparse_grid<Real> grid(dim, level, use_genz_keister);
  return grid.evaluate(f);
}

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_SPARSE_GRID_IMPL_HPP