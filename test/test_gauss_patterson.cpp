// Test Gauss-Patterson quadrature rules for sparse grids
// Copyright 2025 Boost
// Distributed under the Boost Software License, Version 1.0.

#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/detail/cc_rules.hpp>
#include <boost/math/cubature/detail/gauss_patterson_nodes.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace boost::math::cubature;

// Test 1D Gauss-Patterson rules
void test_1d_rules() {
  std::cout << "=== Testing 1D Gauss-Patterson Rules ===\n\n";
  
  using Real = double;
  
  // Test polynomial exactness for each level
  for (std::size_t level = 0; level <= 4; ++level) {
    auto rule = detail::gauss_patterson_nodes<Real>::get_nodes_weights(level);
    const auto& nodes = rule.first;
    const auto& weights = rule.second;
    
    std::cout << "Level " << level << ": " << nodes.size() << " nodes\n";
    std::cout << "Expected polynomial degree: " 
              << detail::gauss_patterson_nodes<Real>::polynomial_degree(level) << "\n";
    
    // Test integral of 1 (should be 2 for [-1,1])
    Real sum_weights = 0;
    for (const auto& w : weights) {
      sum_weights += w;
    }
    std::cout << "  Sum of weights: " << sum_weights << " (should be 2.0)\n";
    
    // Test integral of x^2 (should be 2/3)
    Real integral_x2 = 0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      integral_x2 += weights[i] * nodes[i] * nodes[i];
    }
    Real exact_x2 = Real(2)/Real(3);
    std::cout << "  Integral of x^2: " << integral_x2 
              << " (exact: " << exact_x2 << ", error: " 
              << std::abs(integral_x2 - exact_x2) << ")\n";
    
    // Test integral of x^4 (should be 2/5)
    Real integral_x4 = 0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      Real x2 = nodes[i] * nodes[i];
      integral_x4 += weights[i] * x2 * x2;
    }
    Real exact_x4 = Real(2)/Real(5);
    std::cout << "  Integral of x^4: " << integral_x4 
              << " (exact: " << exact_x4 << ", error: " 
              << std::abs(integral_x4 - exact_x4) << ")\n";
    
    std::cout << "\n";
  }
}

// Compare Clenshaw-Curtis vs Gauss-Patterson for smooth function
void compare_quadratures() {
  std::cout << "=== Comparing Clenshaw-Curtis vs Gauss-Patterson ===\n\n";
  
  using Real = double;
  
  // Test function: smooth exponential-cosine combination
  auto f = [](const std::vector<Real>& x) -> Real {
    Real sum = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      sum += std::exp(-x[i]*x[i]) * std::cos(2*x[i]);
    }
    return sum;
  };
  
  // Test in 2D
  hypercube<Real> box(2);
  box.lower = {-1, -1};
  box.upper = {1, 1};
  
  std::cout << "2D integration of exp(-x^2-y^2) * (cos(2x) + cos(2y))\n\n";
  
  // Test different levels
  for (unsigned level = 1; level <= 4; ++level) {
    std::cout << "Level " << level << ":\n";
    
    // Clenshaw-Curtis
    auto result_cc = integrate_sparse_grid<Real>(
        f, box, level, detail::quadrature_type::clenshaw_curtis);
    std::cout << "  Clenshaw-Curtis: " << std::scientific << std::setprecision(10) 
              << result_cc.value << " (" << result_cc.evaluations << " evals)\n";
    
    // Gauss-Patterson
    auto result_gp = integrate_sparse_grid<Real>(
        f, box, level, detail::quadrature_type::gauss_patterson);
    std::cout << "  Gauss-Patterson: " << result_gp.value 
              << " (" << result_gp.evaluations << " evals)\n";
    
    // Difference
    std::cout << "  Difference: " << std::abs(result_cc.value - result_gp.value) << "\n";
    
    // For smooth functions, Gauss-Patterson should generally be more accurate
    // with fewer points due to higher polynomial exactness
    if (result_gp.evaluations < result_cc.evaluations) {
      std::cout << "  ✓ Gauss-Patterson uses fewer points\n";
    }
    
    std::cout << "\n";
  }
}

// Test high-dimensional integration
void test_high_dimensional() {
  std::cout << "=== Testing High-Dimensional Integration ===\n\n";
  
  using Real = double;
  
  // Product of cosines - separable function
  auto f = [](const std::vector<Real>& x) -> Real {
    Real prod = 1;
    for (const auto& xi : x) {
      prod *= std::cos(xi);
    }
    return prod;
  };
  
  // Test in increasing dimensions
  for (std::size_t dim : {3, 4, 5}) {
    hypercube<Real> box(dim);
    std::fill(box.lower.begin(), box.lower.end(), Real(0));
    std::fill(box.upper.begin(), box.upper.end(), Real(1));
    
    std::cout << "Dimension " << dim << ":\n";
    
    // Exact value for product of cos on [0,1]^d is (sin(1))^d
    Real exact = std::pow(std::sin(Real(1)), dim);
    std::cout << "  Exact value: " << exact << "\n";
    
    // Level 2 sparse grid
    auto result_cc = integrate_sparse_grid<Real>(
        f, box, 2, detail::quadrature_type::clenshaw_curtis);
    auto result_gp = integrate_sparse_grid<Real>(
        f, box, 2, detail::quadrature_type::gauss_patterson);
    
    std::cout << "  Clenshaw-Curtis: " << result_cc.value 
              << " (error: " << std::abs(result_cc.value - exact) 
              << ", evals: " << result_cc.evaluations << ")\n";
    std::cout << "  Gauss-Patterson: " << result_gp.value 
              << " (error: " << std::abs(result_gp.value - exact) 
              << ", evals: " << result_gp.evaluations << ")\n";
    
    std::cout << "\n";
  }
}

int main() {
  test_1d_rules();
  compare_quadratures();
  test_high_dimensional();
  
  std::cout << "=== Test Complete ===\n";
  
  return 0;
}