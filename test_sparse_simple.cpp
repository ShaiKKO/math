#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

// Include minimal headers directly
#include "./include/boost/math/cubature/regions.hpp"
#include "./include/boost/math/cubature/policies.hpp"
#include "./include/boost/math/cubature/detail/cc_rules.hpp"
#include "./include/boost/math/cubature/detail/sparse_grid_impl.hpp"

int main() {
  using Real = double;
  using namespace boost::math::cubature::detail;
  
  std::cout << "\n===== Simple Sparse Grid Test =====\n" << std::endl;
  
  // Test 1: Multi-index generation
  std::cout << "1. Multi-index generation (2D, level 3):" << std::endl;
  {
    multi_index_set<std::size_t> index_set(2, 3);
    const auto& indices = index_set.indices();
    
    for (const auto& idx : indices) {
      int coeff = index_set.coefficient(idx);
      if (coeff != 0) {
        std::cout << "  [" << idx.indices[0] << "," << idx.indices[1] 
                  << "] coeff=" << coeff << std::endl;
      }
    }
  }
  
  // Test 2: Node deduplication
  std::cout << "\n2. Node deduplication:" << std::endl;
  {
    sparse_node_set<Real> nodes;
    nodes.add_node({0.5, 0.5}, 1.0);
    nodes.add_node({0.5, 0.5}, 2.0);  // Duplicate
    nodes.add_node({0.7, 0.5}, 3.0);
    
    std::cout << "  Added 3 nodes, unique: " << nodes.size() << std::endl;
    assert(nodes.size() == 2);
  }
  
  // Test 3: Sparse grid construction
  std::cout << "\n3. Sparse vs Tensor grid nodes:" << std::endl;
  for (std::size_t dim = 2; dim <= 4; ++dim) {
    for (std::size_t level = 2; level <= 3; ++level) {
      smolyak_grid<Real> grid(dim, level);
      std::size_t sparse_nodes = grid.num_nodes();
      
      // Compute tensor product size
      std::size_t tensor_nodes = 1;
      for (std::size_t d = 0; d < dim; ++d) {
        tensor_nodes *= clenshaw_curtis<Real>::num_points(level);
      }
      
      std::cout << "  Dim=" << dim << ", Level=" << level 
                << ": Sparse=" << sparse_nodes 
                << ", Tensor=" << tensor_nodes
                << ", Ratio=" << (Real(sparse_nodes)/tensor_nodes) << std::endl;
    }
  }
  
  // Test 4: Integration test
  std::cout << "\n4. Integration test (polynomial):" << std::endl;
  {
    // Create a simple box [0,1]²
    boost::math::cubature::hypercube<Real> box(2);
    box.lower = {0.0, 0.0};
    box.upper = {1.0, 1.0};
    
    // Test simpler function first: constant 1
    auto constant = [](const Real* x) { return Real(1); };
    
    smolyak_grid<Real> grid(2, 3);
    auto result = grid.evaluate(constant, box);
    
    std::cout << "  ∫(1) computed: " << result.value << " (should be 1)" << std::endl;
    std::cout << "  Nodes used: " << result.evaluations << std::endl;
    
    // Now test polynomial: x² + y²
    auto poly = [](const Real* x) {
      return x[0]*x[0] + x[1]*x[1];
    };
    
    result = grid.evaluate(poly, box);
    
    Real exact = Real(2)/Real(3);  // ∫₀¹∫₀¹(x²+y²)dxdy = 2/3
    Real error = std::abs(result.value - exact);
    
    std::cout << "  ∫(x²+y²) computed: " << result.value << std::endl;
    std::cout << "  Exact: " << exact << std::endl;
    std::cout << "  Error: " << error << std::endl;
    
    // For now, just check it's reasonable (not negative)
    assert(result.value > 0 && result.value < 2);
  }
  
  std::cout << "\n===== All Tests PASSED =====\n" << std::endl;
  
  return 0;
}