#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "./include/boost/math/cubature/detail/cc_rules.hpp"

using Real = double;
using namespace boost::math::cubature::detail;

void test_cc_integration() {
  std::cout << std::fixed << std::setprecision(15);
  
  // Test CC rules on [-1,1] first
  for (std::size_t level = 0; level <= 3; ++level) {
    std::cout << "\n=== CC Level " << level << " ===" << std::endl;
    
    auto nodes = clenshaw_curtis<Real>::get_nodes(level);
    auto weights = clenshaw_curtis<Real>::get_weights(level);
    
    std::cout << "Nodes and weights on [-1,1]:" << std::endl;
    Real weight_sum = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      std::cout << "  x=" << nodes[i] << " w=" << weights[i] << std::endl;
      weight_sum += weights[i];
    }
    std::cout << "Weight sum: " << weight_sum << " (should be 2 for [-1,1])" << std::endl;
    
    // Test integration of x^2 on [-1,1]
    Real integral = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      integral += weights[i] * nodes[i] * nodes[i];
    }
    Real exact = 2.0/3.0;  // ∫₋₁¹ x² dx = 2/3
    std::cout << "∫₋₁¹ x² dx = " << integral << " (exact: " << exact << ")" << std::endl;
    std::cout << "Error: " << std::abs(integral - exact) << std::endl;
    
    // Now test on [0,1]
    std::cout << "\nTransformed to [0,1]:" << std::endl;
    Real integral_01 = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      Real x = (nodes[i] + 1) / 2;  // Transform to [0,1]
      Real w = weights[i] / 2;       // Scale weight
      std::cout << "  x=" << x << " w=" << w << std::endl;
      integral_01 += w * x * x;
    }
    Real exact_01 = 1.0/3.0;  // ∫₀¹ x² dx = 1/3
    std::cout << "∫₀¹ x² dx = " << integral_01 << " (exact: " << exact_01 << ")" << std::endl;
    std::cout << "Error: " << std::abs(integral_01 - exact_01) << std::endl;
  }
}

int main() {
  test_cc_integration();
  return 0;
}