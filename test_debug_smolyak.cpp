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

void debug_tensor_product() {
  std::cout << "\n=== DEBUG: Tensor Product Construction ===" << std::endl;
  
  // Test simple 2D tensor product for multi-index [1,1]
  multi_index<std::size_t> idx({1, 1});
  
  // Get 1D rules
  quadrature_rule_1d<Real> rule_x(1);  // Level 1: 3 points
  quadrature_rule_1d<Real> rule_y(1);  // Level 1: 3 points
  
  std::cout << "Level 1 CC nodes: ";
  for (auto n : rule_x.nodes) std::cout << n << " ";
  std::cout << std::endl;
  
  std::cout << "Level 1 CC weights: ";
  for (auto w : rule_x.weights) std::cout << w << " ";
  std::cout << std::endl;
  
  // Build tensor product manually
  std::cout << "\nTensor product [1,1] points:" << std::endl;
  Real weight_sum = 0.0;
  for (std::size_t i = 0; i < rule_x.size(); ++i) {
    for (std::size_t j = 0; j < rule_y.size(); ++j) {
      Real x = (rule_x.nodes[i] + 1) / 2;  // Transform to [0,1]
      Real y = (rule_y.nodes[j] + 1) / 2;
      Real w = rule_x.weights[i] * rule_y.weights[j] / 4;  // Scale for [0,1]²
      
      std::cout << "  (" << x << ", " << y << ") weight=" << w << std::endl;
      weight_sum += w;
    }
  }
  std::cout << "Weight sum: " << weight_sum << " (should be 1)" << std::endl;
}

void debug_multi_index_coefficients() {
  std::cout << "\n=== DEBUG: Multi-Index Coefficients ===" << std::endl;
  
  multi_index_set<std::size_t> index_set(2, 3);
  const auto& indices = index_set.indices();
  
  Real total_weight = 0.0;
  for (const auto& idx : indices) {
    int coeff = index_set.coefficient(idx);
    if (coeff != 0) {
      std::cout << "[" << idx.indices[0] << "," << idx.indices[1] 
                << "] |i|=" << idx.sum 
                << " coeff=" << coeff << std::endl;
      
      // Count tensor product points
      quadrature_rule_1d<Real> rule_x(idx.indices[0]);
      quadrature_rule_1d<Real> rule_y(idx.indices[1]);
      std::size_t n_points = rule_x.size() * rule_y.size();
      
      std::cout << "  Tensor size: " << rule_x.size() << " x " << rule_y.size() 
                << " = " << n_points << " points" << std::endl;
    }
  }
}

void debug_node_accumulation() {
  std::cout << "\n=== DEBUG: Node Weight Accumulation ===" << std::endl;
  
  // Test the problematic Kahan summation
  sparse_node_set<Real> nodes;
  
  // Add a node
  std::vector<Real> point = {0.5, 0.5};
  nodes.add_node(point, 1.0);
  std::cout << "Added weight 1.0 to (0.5, 0.5)" << std::endl;
  
  // Add same node again with different weight
  nodes.add_node(point, 2.0);
  std::cout << "Added weight 2.0 to (0.5, 0.5)" << std::endl;
  
  // Check accumulated weight
  const auto& final_nodes = nodes.get_nodes();
  for (const auto& node : final_nodes) {
    std::cout << "Final node (" << node.point[0] << ", " << node.point[1] 
              << ") has weight " << node.weight << std::endl;
  }
  
  if (std::abs(final_nodes[0].weight - 3.0) > 1e-10) {
    std::cout << "ERROR: Weight should be 3.0, got " << final_nodes[0].weight << std::endl;
    std::cout << "Bug is in line 175 of sparse_grid_impl.hpp!" << std::endl;
  }
}

void debug_simple_integral() {
  std::cout << "\n=== DEBUG: Simple Integral Test ===" << std::endl;
  
  // Create box [0,1]²
  boost::math::cubature::hypercube<Real> box(2);
  box.lower = {0.0, 0.0};
  box.upper = {1.0, 1.0};
  
  // Test constant function = 1
  auto constant = [](const Real* x) { return Real(1); };
  
  // Level 1 should be simple
  smolyak_grid<Real> grid(2, 1);
  auto result = grid.evaluate(constant, box);
  
  std::cout << "Level 1, constant function:" << std::endl;
  std::cout << "  Result: " << result.value << " (should be 1)" << std::endl;
  std::cout << "  Nodes: " << result.evaluations << std::endl;
  
  // Check the actual nodes and weights
  const auto& nodes = grid.get_nodes();
  Real weight_sum = 0.0;
  std::cout << "  Node details:" << std::endl;
  for (const auto& node : nodes) {
    std::cout << "    (" << node.point[0] << ", " << node.point[1] 
              << ") weight=" << node.weight << std::endl;
    weight_sum += node.weight;
  }
  std::cout << "  Total weight sum: " << weight_sum << " (should be 1)" << std::endl;
}

int main() {
  std::cout << std::fixed << std::setprecision(6);
  
  debug_tensor_product();
  debug_multi_index_coefficients();
  debug_node_accumulation();
  debug_simple_integral();
  
  return 0;
}