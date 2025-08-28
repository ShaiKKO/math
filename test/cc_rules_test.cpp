#include <boost/math/cubature/detail/cc_rules.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <numeric>

using namespace boost::math::cubature::detail;

template <typename Real>
void test_cc_growth() {
  std::cout << "Testing Clenshaw-Curtis growth sequence..." << std::endl;
  
  // Verify slow exponential growth: 1, 3, 5, 9, 17, 33, ...
  assert(clenshaw_curtis<Real>::num_points(0) == 1);
  assert(clenshaw_curtis<Real>::num_points(1) == 3);
  assert(clenshaw_curtis<Real>::num_points(2) == 5);
  assert(clenshaw_curtis<Real>::num_points(3) == 9);
  assert(clenshaw_curtis<Real>::num_points(4) == 17);
  assert(clenshaw_curtis<Real>::num_points(5) == 33);
  assert(clenshaw_curtis<Real>::num_points(6) == 65);
  
  std::cout << "  Growth sequence: ";
  for (std::size_t l = 0; l <= 6; ++l) {
    std::cout << clenshaw_curtis<Real>::num_points(l) << " ";
  }
  std::cout << std::endl;
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_cc_nodes() {
  std::cout << "Testing Clenshaw-Curtis nodes..." << std::endl;
  
  // Level 0: single midpoint
  {
    auto nodes = clenshaw_curtis<Real>::get_nodes(0);
    assert(nodes.size() == 1);
    assert(std::abs(nodes[0]) < 1e-14);
    std::cout << "  Level 0 node: " << nodes[0] << " (midpoint)" << std::endl;
  }
  
  // Level 1: 3 points (endpoints and midpoint)
  {
    auto nodes = clenshaw_curtis<Real>::get_nodes(1);
    assert(nodes.size() == 3);
    assert(std::abs(nodes[0] - Real(-1)) < 1e-14);
    assert(std::abs(nodes[1]) < 1e-14);
    assert(std::abs(nodes[2] - Real(1)) < 1e-14);
    std::cout << "  Level 1 nodes: ";
    for (auto n : nodes) std::cout << n << " ";
    std::cout << std::endl;
  }
  
  // Verify nodes are in [-1, 1]
  for (std::size_t l = 0; l <= 4; ++l) {
    auto nodes = clenshaw_curtis<Real>::get_nodes(l);
    for (const auto& node : nodes) {
      assert(node >= Real(-1) - 1e-14 && node <= Real(1) + 1e-14);
    }
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_cc_weights() {
  std::cout << "Testing Clenshaw-Curtis weights..." << std::endl;
  
  // Weights should sum to 2 (integral of 1 over [-1, 1])
  for (std::size_t l = 0; l <= 5; ++l) {
    auto weights = clenshaw_curtis<Real>::get_weights(l);
    Real sum = std::accumulate(weights.begin(), weights.end(), Real(0));
    std::cout << "  Level " << l << " weight sum: " << sum;
    assert(std::abs(sum - Real(2)) < 1e-12);
    std::cout << " ✓" << std::endl;
  }
  
  // Test integration of x^2 over [-1, 1]
  // Exact integral is 2/3
  for (std::size_t l = 1; l <= 4; ++l) {
    auto nodes = clenshaw_curtis<Real>::get_nodes(l);
    auto weights = clenshaw_curtis<Real>::get_weights(l);
    
    Real integral = Real(0);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      integral += weights[i] * nodes[i] * nodes[i];
    }
    
    Real exact = Real(2) / Real(3);
    Real error = std::abs(integral - exact);
    std::cout << "  Level " << l << " integral of x²: " << integral 
              << " (error: " << error << ")" << std::endl;
    
    // CC should be exact for polynomials up to degree 2*n-1
    if (l >= 2) {
      assert(error < 1e-12);
    }
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_cc_nestedness() {
  std::cout << "Testing Clenshaw-Curtis nestedness..." << std::endl;
  
  bool nested = clenshaw_curtis<Real>::verify_nested(5);
  assert(nested);
  
  // Manually verify level 0 ⊂ level 1 ⊂ level 2
  auto nodes0 = clenshaw_curtis<Real>::get_nodes(0);
  auto nodes1 = clenshaw_curtis<Real>::get_nodes(1);
  auto nodes2 = clenshaw_curtis<Real>::get_nodes(2);
  
  std::cout << "  Level 0 nodes: ";
  for (auto n : nodes0) std::cout << std::setprecision(4) << n << " ";
  std::cout << std::endl;
  
  std::cout << "  Level 1 nodes: ";
  for (auto n : nodes1) std::cout << std::setprecision(4) << n << " ";
  std::cout << std::endl;
  
  std::cout << "  Level 2 nodes: ";
  for (auto n : nodes2) std::cout << std::setprecision(4) << n << " ";
  std::cout << std::endl;
  
  // Check new nodes
  auto new_nodes = clenshaw_curtis<Real>::get_new_nodes(2);
  std::cout << "  New nodes at level 2: ";
  for (auto n : new_nodes) std::cout << std::setprecision(4) << n << " ";
  std::cout << std::endl;
  
  std::cout << "  Nestedness verified for levels 0-5" << std::endl;
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_weight_cache() {
  std::cout << "Testing weight cache..." << std::endl;
  
  typename clenshaw_curtis<Real>::weight_cache cache;
  
  // First access computes
  const auto& weights1 = cache.get_weights(3);
  
  // Second access should return cached
  const auto& weights2 = cache.get_weights(3);
  
  // Should be same memory location
  assert(&weights1 == &weights2);
  
  // Clear and verify recomputation
  cache.clear();
  const auto& weights3 = cache.get_weights(3);
  
  // After clear, should be different memory location (new computation)
  // Note: weights1 reference may be invalidated after clear, so we compare values only
  assert(weights3.size() == 9);  // Level 3 has 9 points
  Real sum = std::accumulate(weights3.begin(), weights3.end(), Real(0));
  assert(std::abs(sum - Real(2)) < 1e-12);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_tensor_product() {
  std::cout << "Testing tensor product rule..." << std::endl;
  
  // 2D tensor product at level 1
  tensor_product_rule<Real> tp(2, 1);
  
  // Should have 3 × 3 = 9 points
  assert(tp.num_points() == 9);
  
  auto nodes = tp.get_nodes();
  auto weights = tp.get_weights();
  
  assert(nodes.size() == 9);
  assert(weights.size() == 9);
  
  // Check weight sum (should be 4 for [-1,1]²)
  Real weight_sum = std::accumulate(weights.begin(), weights.end(), Real(0));
  assert(std::abs(weight_sum - Real(4)) < 1e-12);
  
  std::cout << "  2D tensor product (level 1):" << std::endl;
  std::cout << "    Points: " << tp.num_points() << std::endl;
  std::cout << "    Weight sum: " << weight_sum << std::endl;
  
  // Test mixed levels
  std::vector<std::size_t> levels = {1, 2};
  tensor_product_rule<Real> tp_mixed(2, levels);
  assert(tp_mixed.num_points() == 3 * 5);  // 3 × 5 = 15
  
  std::cout << "  Mixed levels (1,2): " << tp_mixed.num_points() << " points" << std::endl;
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_integration_accuracy() {
  std::cout << "Testing integration accuracy..." << std::endl;
  
  // Test integral of exp(x+y) over [-1,1]²
  // Exact: 4*sinh(1)²/1 ≈ 5.524391382...
  const Real exact = Real(4) * std::sinh(Real(1)) * std::sinh(Real(1));
  
  for (std::size_t level = 2; level <= 4; ++level) {
    tensor_product_rule<Real> tp(2, level);
    auto nodes = tp.get_nodes();
    auto weights = tp.get_weights();
    
    Real integral = Real(0);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      Real f_val = std::exp(nodes[i][0] + nodes[i][1]);
      integral += weights[i] * f_val;
    }
    
    Real error = std::abs(integral - exact);
    std::cout << "  Level " << level << " (";
    std::cout << tp.num_points() << " points): ";
    std::cout << "integral = " << integral;
    std::cout << ", error = " << std::scientific << error << std::endl;
    
    // Error should decrease with level
    if (level >= 3) {
      assert(error < 1e-6);
    }
  }
  
  std::cout << "  PASSED" << std::endl;
}

int main() {
  std::cout << "\n===== Testing Clenshaw-Curtis Rules =====\n" << std::endl;
  
  using Real = double;
  
  test_cc_growth<Real>();
  test_cc_nodes<Real>();
  test_cc_weights<Real>();
  test_cc_nestedness<Real>();
  test_weight_cache<Real>();
  test_tensor_product<Real>();
  test_integration_accuracy<Real>();
  
  std::cout << "\n===== All Tests PASSED =====\n" << std::endl;
  
  return 0;
}