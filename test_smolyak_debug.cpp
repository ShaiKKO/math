#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "./include/boost/math/cubature/regions.hpp"
#include "./include/boost/math/cubature/policies.hpp"
#include "./include/boost/math/cubature/detail/cc_rules.hpp"
#include "./include/boost/math/cubature/detail/sparse_grid_impl.hpp"

using Real = double;
using namespace boost::math::cubature::detail;

void debug_1d_smolyak() {
  std::cout << "\n=== 1D Smolyak Debug (Level 1) ===" << std::endl;
  
  // For 1D, level 1, Smolyak should give us CC level 1 (3 points)
  // Multi-indices for d=1, l=1: |i| must be between 1 and 1
  // So only i=[1] with coefficient = (-1)^(1-1) * C(0,0) = 1
  
  multi_index_set<std::size_t> index_set(1, 1);
  const auto& indices = index_set.indices();
  
  std::cout << "Multi-indices:" << std::endl;
  for (const auto& idx : indices) {
    int coeff = index_set.coefficient(idx);
    std::cout << "  [" << idx.indices[0] << "] |i|=" << idx.sum 
              << " coeff=" << coeff << std::endl;
  }
  
  // Build the sparse grid manually
  smolyak_grid<Real> grid(1, 1);
  const auto& nodes = grid.get_nodes();
  
  std::cout << "\nSparse grid nodes and weights:" << std::endl;
  Real weight_sum = 0.0;
  for (const auto& node : nodes) {
    std::cout << "  x=" << node.point[0] << " weight=" << node.weight << std::endl;
    weight_sum += node.weight;
  }
  std::cout << "Total weight: " << weight_sum << " (should be 1)" << std::endl;
  
  // Test quadratic x^2
  boost::math::cubature::hypercube<Real> box(1);
  box.lower = {0.0};
  box.upper = {1.0};
  
  auto quad_func = [](const Real* x) { return x[0] * x[0]; };
  auto result = grid.evaluate(quad_func, box);
  
  std::cout << "\n∫₀¹ x² dx = " << result.value << " (exact: 1/3 = " << (1.0/3.0) << ")" << std::endl;
  
  // Manual calculation
  Real manual_sum = 0.0;
  for (const auto& node : nodes) {
    Real x = node.point[0];  // Already in [0,1]
    Real fx = x * x;
    manual_sum += fx * node.weight;
  }
  std::cout << "Manual calculation: " << manual_sum << std::endl;
  
  // Check CC rule directly
  std::cout << "\nDirect CC rule (level 1):" << std::endl;
  quadrature_rule_1d<Real> cc_rule(1);
  for (std::size_t i = 0; i < cc_rule.size(); ++i) {
    Real x_cc = (cc_rule.nodes[i] + 1) / 2;  // Transform to [0,1]
    Real w_cc = cc_rule.weights[i] / 2;      // Scale weight
    std::cout << "  x=" << x_cc << " weight=" << w_cc << std::endl;
  }
}

void debug_smolyak_formula() {
  std::cout << "\n=== Smolyak Formula Debug ===" << std::endl;
  
  // The Smolyak formula for 2D, level 2:
  // A(2,2) = Σ (-1)^(2-|i|) * C(1, 2-|i|) * (U^i₁ ⊗ U^i₂)
  // where |i| ranges from 1 to 2
  
  std::cout << "2D, Level 2 multi-indices:" << std::endl;
  multi_index_set<std::size_t> index_set(2, 2);
  for (const auto& idx : index_set.indices()) {
    int coeff = index_set.coefficient(idx);
    if (coeff != 0) {
      std::cout << "  [" << idx.indices[0] << "," << idx.indices[1] 
                << "] |i|=" << idx.sum 
                << " coeff=" << coeff;
      
      // Show tensor product size
      quadrature_rule_1d<Real> rule_x(idx.indices[0]);
      quadrature_rule_1d<Real> rule_y(idx.indices[1]);
      std::cout << " (tensor: " << rule_x.size() << "×" << rule_y.size() << ")" << std::endl;
    }
  }
  
  // Build and check nodes
  smolyak_grid<Real> grid(2, 2);
  std::cout << "\nTotal unique nodes: " << grid.num_nodes() << std::endl;
  
  // Check weight sum
  const auto& nodes = grid.get_nodes();
  Real weight_sum = 0.0;
  for (const auto& node : nodes) {
    weight_sum += node.weight;
  }
  std::cout << "Total weight sum: " << weight_sum << " (should be 1)" << std::endl;
}

int main() {
  std::cout << std::fixed << std::setprecision(10);
  
  debug_1d_smolyak();
  debug_smolyak_formula();
  
  return 0;
}