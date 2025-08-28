#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <cassert>

// Include sparse grid directly without the full framework
#include "boost/math/cubature/regions.hpp"
#include "boost/math/cubature/policies.hpp"
#include "boost/math/cubature/detail/cc_rules.hpp"
#include "boost/math/cubature/detail/sparse_grid_impl.hpp"

using namespace boost::math::cubature;
using namespace boost::math::cubature::detail;

// Test functions from the full test suite
template <typename Real>
void test_multi_index_generation() {
  std::cout << "\n1. Testing multi-index generation..." << std::endl;
  
  // Test 2D, level 3
  {
    multi_index_set<std::size_t> index_set(2, 3);
    const auto& indices = index_set.indices();
    
    // Should have indices with sum between 2 and 3
    std::size_t expected_count = 7;  // [0,2], [0,3], [1,1], [1,2], [2,0], [2,1], [3,0]
    assert(indices.size() == expected_count);
    
    // Check coefficients sum to something sensible
    int coeff_sum = 0;
    for (const auto& idx : indices) {
      coeff_sum += index_set.coefficient(idx);
    }
    assert(coeff_sum != 0);  // Should be non-zero for proper combination
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_node_deduplication() {
  std::cout << "\n2. Testing node deduplication..." << std::endl;
  
  sparse_node_set<Real> nodes;
  
  // Add same node multiple times
  nodes.add_node({0.5, 0.5}, 1.0);
  nodes.add_node({0.5, 0.5}, 2.0);
  nodes.add_node({0.5, 0.5}, 3.0);
  
  assert(nodes.size() == 1);
  
  // Check accumulated weight
  const auto& final_nodes = nodes.get_nodes();
  assert(final_nodes.size() == 1);
  assert(std::abs(final_nodes[0].weight - 6.0) < 1e-14);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_sparse_vs_tensor() {
  std::cout << "\n3. Testing sparse vs tensor grid sizes..." << std::endl;
  
  for (std::size_t dim = 2; dim <= 4; ++dim) {
    for (std::size_t level = 2; level <= 3; ++level) {
      smolyak_grid<Real> grid(dim, level);
      std::size_t sparse_nodes = grid.num_nodes();
      
      // Compute tensor product size
      std::size_t tensor_nodes = 1;
      for (std::size_t d = 0; d < dim; ++d) {
        tensor_nodes *= clenshaw_curtis<Real>::num_points(level);
      }
      
      // Sparse should be significantly smaller
      assert(sparse_nodes < tensor_nodes);
      
      Real ratio = Real(sparse_nodes) / Real(tensor_nodes);
      assert(ratio < 1.0);  // Always fewer nodes
      
      if (dim >= 3 && level >= 3) {
        assert(ratio < 0.2);  // Should be very efficient in high dimensions
      }
    }
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_polynomial_exactness() {
  std::cout << "\n4. Testing polynomial exactness..." << std::endl;
  
  // Test 2D quadratic
  hypercube<Real> box(2);
  box.lower = {0.0, 0.0};
  box.upper = {1.0, 1.0};
  
  // f(x,y) = x^2 + y^2
  auto quadratic = [](const Real* x) {
    return x[0]*x[0] + x[1]*x[1];
  };
  
  smolyak_grid<Real> grid(2, 2);  // Level 2 should be exact for degree 2
  auto result = grid.evaluate(quadratic, box);
  
  Real exact = Real(2) / Real(3);  // integral of x^2 + y^2 over [0,1]^2
  Real error = std::abs(result.value - exact);
  
  assert(error < 1e-14);
  
  std::cout << "  Quadratic exact to " << error << std::endl;
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_smooth_function() {
  std::cout << "\n5. Testing smooth function convergence..." << std::endl;
  
  hypercube<Real> box(2);
  box.lower = {0.0, 0.0};
  box.upper = {1.0, 1.0};
  
  // Gaussian-like smooth function
  auto smooth = [](const Real* x) {
    Real r2 = (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[1]-0.5);
    return std::exp(-4*r2);
  };
  
  Real prev_error = 1.0;
  for (std::size_t level = 2; level <= 4; ++level) {
    smolyak_grid<Real> grid(2, level);
    auto result = grid.evaluate(smooth, box);
    
    // Compare with higher accuracy reference
    smolyak_grid<Real> ref_grid(2, 6);
    auto ref_result = ref_grid.evaluate(smooth, box);
    
    Real error = std::abs(result.value - ref_result.value);
    
    // Should converge
    if (level > 2) {
      assert(error < prev_error);
    }
    prev_error = error;
    
    std::cout << "  Level " << level << ": error = " << error << std::endl;
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_performance() {
  std::cout << "\n6. Testing performance scaling..." << std::endl;
  
  hypercube<Real> box(4);
  for (std::size_t i = 0; i < 4; ++i) {
    box.lower[i] = -1.0;
    box.upper[i] = 1.0;
  }
  
  auto test_func = [](const Real* x) {
    return x[0] + x[1]*x[1] + x[2]*x[2]*x[2] + x[3]*x[3]*x[3]*x[3];
  };
  
  for (std::size_t level = 3; level <= 5; ++level) {
    auto start = std::chrono::high_resolution_clock::now();
    
    smolyak_grid<Real> grid(4, level);
    auto res = grid.evaluate(test_func, box);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    std::cout << "  Level " << level << ": " 
              << res.evaluations << " nodes, "
              << duration.count() << " μs" << std::endl;
    
    // Should complete quickly even for high dimensions
    assert(duration.count() < 100000);  // < 100ms
  }
  
  std::cout << "  PASSED" << std::endl;
}

int main() {
  std::cout << "\n===== Full Smolyak Sparse Grid Test Suite =====\n";
  std::cout << "Testing complete integration with Boost.Math::cubature" << std::endl;
  
  using Real = double;
  
  test_multi_index_generation<Real>();
  test_node_deduplication<Real>();
  test_sparse_vs_tensor<Real>();
  test_polynomial_exactness<Real>();
  test_smooth_function<Real>();
  test_performance<Real>();
  
  std::cout << "\n===== All Full Tests PASSED =====\n" << std::endl;
  
  return 0;
}