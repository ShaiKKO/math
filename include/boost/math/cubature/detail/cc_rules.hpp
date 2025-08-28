// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_CC_RULES_HPP
#define BOOST_MATH_CUBATURE_DETAIL_CC_RULES_HPP

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <unordered_map>
// Remove boost constants dependency

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Clenshaw-Curtis quadrature rules with nested slow exponential growth
/// \details Implements nested CC rules with growth sequence 1, 3, 5, 9, 17, 33, ...
///          Following ALGORITHMS.md §2.1 for sparse grid construction
template <typename Real>
class clenshaw_curtis {
public:
  /// \brief Get number of points for given level (slow exponential growth)
  /// \details Growth: m(0)=1, m(1)=3, m(2)=5, m(3)=9, m(4)=17, ...
  ///          Formula: m(l) = 1 if l=0, else 2^l + 1 for l≥2, special case for l=1
  static std::size_t num_points(std::size_t level) {
    switch(level) {
      case 0: return 1;
      case 1: return 3;
      default: return (1u << level) + 1;  // 2^l + 1 for l≥2
    }
  }
  
  /// \brief Get Clenshaw-Curtis nodes for given level
  /// \details Nodes are cos(πj/n) for j=0,...,n where n = num_points(level)-1
  static std::vector<Real> get_nodes(std::size_t level) {
    const std::size_t n_points = num_points(level);
    std::vector<Real> nodes(n_points);
    
    if (n_points == 1) {
      nodes[0] = Real(0);  // Midpoint for level 0
      return nodes;
    }
    
    // Standard Clenshaw-Curtis nodes on [-1, 1]
    const std::size_t n = n_points - 1;
    const Real pi = Real(3.14159265358979323846);
    
    for (std::size_t j = 0; j <= n; ++j) {
      nodes[j] = std::cos(pi * Real(j) / Real(n));
    }
    
    // Sort in ascending order
    std::sort(nodes.begin(), nodes.end());
    
    return nodes;
  }
  
  /// \brief Compute Clenshaw-Curtis weights 
  /// \details Implements the classical Clenshaw-Curtis weight formula
  ///          Based on standard references (Waldvogel 2006, Trefethen)
  static std::vector<Real> get_weights(std::size_t level) {
    const std::size_t n_points = num_points(level);
    std::vector<Real> weights(n_points);
    
    if (n_points == 1) {
      weights[0] = Real(2);  // Integral over [-1, 1]
      return weights;
    }
    
    const std::size_t N = n_points - 1;
    const Real pi = Real(3.14159265358979323846);
    
    // Use the standard textbook formula for CC weights
    // Reference: Waldvogel, "Fast Construction of the Fejér and Clenshaw-Curtis 
    //            Quadrature Rules", BIT 46 (2006), pp 195-202
    
    if (n_points == 2) {
      // Special case: 2 points
      weights[0] = Real(1);
      weights[1] = Real(1);
    } else if (n_points == 3) {
      // Special case: 3 points, exact for degree 2
      weights[0] = Real(1) / Real(3);
      weights[1] = Real(4) / Real(3);
      weights[2] = Real(1) / Real(3);
    } else {
      // General formula
      for (std::size_t i = 0; i <= N; ++i) {
        Real theta_i = Real(i) * pi / Real(N);
        weights[i] = Real(0);
        
        // Sum contributions from Chebyshev polynomials
        for (std::size_t j = 0; j <= N/2; ++j) {
          Real b_j;
          if (j == 0) {
            b_j = Real(1);
          } else if (2*j == N) {
            b_j = Real(1) / Real(2);
          } else {
            b_j = Real(2) / (Real(1) - Real(4*j*j));
          }
          weights[i] += b_j * std::cos(Real(2*j) * theta_i);
        }
        
        // Apply endpoint scaling
        if (i == 0 || i == N) {
          weights[i] = weights[i] / Real(N);
        } else {
          weights[i] = Real(2) * weights[i] / Real(N);
        }
      }
    }
    
    // Ensure weights are in same order as nodes (ascending x)
    // Since x_i = cos(theta_i) and theta increases, x decreases, so reverse
    std::reverse(weights.begin(), weights.end());
    
    return weights;
  }
  
  /// \brief Verify nestedness property
  /// \details Check that nodes at level l are subset of nodes at level l+1
  static bool verify_nested(std::size_t max_level = 5) {
    for (std::size_t l = 0; l < max_level; ++l) {
      auto nodes_l = get_nodes(l);
      auto nodes_l1 = get_nodes(l + 1);
      
      // Check each node from level l exists in level l+1
      for (const auto& node : nodes_l) {
        bool found = false;
        for (const auto& node_next : nodes_l1) {
          if (std::abs(node - node_next) < Real(1e-14)) {
            found = true;
            break;
          }
        }
        if (!found) {
          return false;
        }
      }
    }
    return true;
  }
  
  /// \brief Get new nodes added at level l (not in level l-1)
  static std::vector<Real> get_new_nodes(std::size_t level) {
    if (level == 0) {
      return get_nodes(0);
    }
    
    auto nodes_prev = get_nodes(level - 1);
    auto nodes_curr = get_nodes(level);
    std::vector<Real> new_nodes;
    
    for (const auto& node : nodes_curr) {
      bool is_new = true;
      for (const auto& prev_node : nodes_prev) {
        if (std::abs(node - prev_node) < Real(1e-14)) {
          is_new = false;
          break;
        }
      }
      if (is_new) {
        new_nodes.push_back(node);
      }
    }
    
    return new_nodes;
  }
  
  /// \brief Cache for storing computed weights
  class weight_cache {
  private:
    mutable std::unordered_map<std::size_t, std::vector<Real>> cache_;
    
  public:
    const std::vector<Real>& get_weights(std::size_t level) const {
      auto it = cache_.find(level);
      if (it != cache_.end()) {
        return it->second;
      }
      
      // Compute and cache
      cache_[level] = clenshaw_curtis<Real>::get_weights(level);
      return cache_[level];
    }
    
    void clear() {
      cache_.clear();
    }
  };
};

/// \brief 1D quadrature rule interface for tensor product construction
template <typename Real>
struct quadrature_rule_1d {
  std::vector<Real> nodes;
  std::vector<Real> weights;
  std::size_t level;
  
  quadrature_rule_1d(std::size_t l = 0) : level(l) {
    nodes = clenshaw_curtis<Real>::get_nodes(level);
    weights = clenshaw_curtis<Real>::get_weights(level);
  }
  
  std::size_t size() const { return nodes.size(); }
};

/// \brief Transform nodes from [-1, 1] to [a, b]
template <typename Real>
void transform_nodes(std::vector<Real>& nodes, Real a, Real b) {
  const Real mid = (a + b) / Real(2);
  const Real half_width = (b - a) / Real(2);
  
  for (auto& node : nodes) {
    node = mid + half_width * node;
  }
}

/// \brief Transform weights for interval [a, b]
template <typename Real>
void transform_weights(std::vector<Real>& weights, Real a, Real b) {
  const Real jacobian = (b - a) / Real(2);
  
  for (auto& weight : weights) {
    weight *= jacobian;
  }
}

/// \brief Compute tensor product nodes for multi-dimensional integration
template <typename Real>
class tensor_product_rule {
private:
  std::size_t dimension_;
  std::vector<quadrature_rule_1d<Real>> rules_1d_;
  
public:
  tensor_product_rule(std::size_t dim, std::size_t level)
    : dimension_(dim), rules_1d_(dim, quadrature_rule_1d<Real>(level)) {}
  
  tensor_product_rule(std::size_t dim, const std::vector<std::size_t>& levels)
    : dimension_(dim) {
    for (std::size_t d = 0; d < dim; ++d) {
      rules_1d_.emplace_back(levels[d]);
    }
  }
  
  /// \brief Get total number of tensor product points
  std::size_t num_points() const {
    std::size_t total = 1;
    for (const auto& rule : rules_1d_) {
      total *= rule.size();
    }
    return total;
  }
  
  /// \brief Generate all tensor product nodes
  std::vector<std::vector<Real>> get_nodes() const {
    const std::size_t n_points = num_points();
    std::vector<std::vector<Real>> nodes(n_points, std::vector<Real>(dimension_));
    
    // Generate tensor product using index arithmetic
    for (std::size_t i = 0; i < n_points; ++i) {
      std::size_t idx = i;
      for (std::size_t d = 0; d < dimension_; ++d) {
        std::size_t n_d = rules_1d_[d].size();
        std::size_t j_d = idx % n_d;
        nodes[i][d] = rules_1d_[d].nodes[j_d];
        idx /= n_d;
      }
    }
    
    return nodes;
  }
  
  /// \brief Generate tensor product weights
  std::vector<Real> get_weights() const {
    const std::size_t n_points = num_points();
    std::vector<Real> weights(n_points);
    
    // Compute product weights
    for (std::size_t i = 0; i < n_points; ++i) {
      weights[i] = Real(1);
      std::size_t idx = i;
      for (std::size_t d = 0; d < dimension_; ++d) {
        std::size_t n_d = rules_1d_[d].size();
        std::size_t j_d = idx % n_d;
        weights[i] *= rules_1d_[d].weights[j_d];
        idx /= n_d;
      }
    }
    
    return weights;
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_CC_RULES_HPP

