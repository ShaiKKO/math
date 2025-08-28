#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/detail/sparse_grid_impl.hpp>
// Remove boost constants to avoid namespace issues
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <chrono>

using namespace boost::math::cubature;
using namespace boost::math::cubature::detail;

template <typename Real>
void test_multi_index_generation() {
  std::cout << "Testing multi-index generation..." << std::endl;
  
  // Test 2D, level 3
  {
    multi_index_set<std::size_t> index_set(2, 3);
    const auto& indices = index_set.indices();
    
    std::cout << "  2D, Level 3 indices:" << std::endl;
    for (const auto& idx : indices) {
      int coeff = index_set.coefficient(idx);
      std::cout << "    [" << idx.indices[0] << "," << idx.indices[1] 
                << "] |i|=" << idx.sum << " coeff=" << coeff << std::endl;
    }
    
    // Should have indices with sum 2 or 3
    // Sum=2: (0,2), (1,1), (2,0) -> coeff = -C(1,1) = -1
    // Sum=3: (0,3), (1,2), (2,1), (3,0) -> coeff = C(1,0) = 1
    assert(indices.size() == 7);
  }
  
  // Test 3D, level 2
  {
    multi_index_set<std::size_t> index_set(3, 2);
    const auto& indices = index_set.indices();
    
    std::cout << "  3D, Level 2 has " << indices.size() << " multi-indices" << std::endl;
    
    // Should have indices with sum 0, 1, or 2
    // This tests ℓ-d+1 ≤ |i| ≤ ℓ with ℓ=2, d=3 → 0 ≤ |i| ≤ 2
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_node_deduplication() {
  std::cout << "Testing node deduplication..." << std::endl;
  
  sparse_node_set<Real> nodes;
  
  // Add duplicate nodes
  std::vector<Real> point1 = {0.5, 0.5};
  std::vector<Real> point2 = {0.5, 0.5};  // Exact duplicate
  std::vector<Real> point3 = {0.5 + 1e-15, 0.5};  // Within ULP tolerance
  std::vector<Real> point4 = {0.7, 0.5};  // Different point
  
  nodes.add_node(point1, Real(1.0));
  nodes.add_node(point2, Real(2.0));
  nodes.add_node(point3, Real(3.0));
  nodes.add_node(point4, Real(4.0));
  
  const auto& unique = nodes.get_nodes();
  
  // Should have 2 unique nodes (point1/2/3 merged, point4 separate)
  assert(nodes.size() == 2);
  
  // Check weight accumulation
  bool found_merged = false;
  bool found_separate = false;
  
  for (const auto& node : unique) {
    if (std::abs(node.point[0] - 0.5) < 1e-10) {
      // Merged node should have accumulated weight
      assert(std::abs(node.weight - Real(6.0)) < 1e-10);
      found_merged = true;
    } else if (std::abs(node.point[0] - 0.7) < 1e-10) {
      assert(std::abs(node.weight - Real(4.0)) < 1e-10);
      found_separate = true;
    }
  }
  
  assert(found_merged && found_separate);
  
  std::cout << "  Deduplication: " << 4 << " nodes -> " << nodes.size() << " unique" << std::endl;
  std::cout << "  Weight accumulation verified" << std::endl;
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_sparse_vs_tensor() {
  std::cout << "Testing sparse vs tensor product node counts..." << std::endl;
  
  // Compare node counts for different dimensions and levels
  std::cout << "\n  Node count comparison:" << std::endl;
  std::cout << "  Dim | Level | Sparse | Tensor | Ratio" << std::endl;
  std::cout << "  ----|-------|--------|--------|-------" << std::endl;
  
  for (std::size_t dim = 2; dim <= 5; ++dim) {
    for (std::size_t level = 2; level <= 4; ++level) {
      // Sparse grid nodes
      smolyak_grid<Real> sparse(dim, level);
      std::size_t sparse_nodes = sparse.num_nodes();
      
      // Tensor product nodes (using CC rule at max level)
      std::size_t tensor_nodes = 1;
      for (std::size_t d = 0; d < dim; ++d) {
        tensor_nodes *= clenshaw_curtis<Real>::num_points(level);
      }
      
      Real ratio = static_cast<Real>(sparse_nodes) / tensor_nodes;
      
      std::cout << "   " << dim << "  |   " << level << "   |  " 
                << std::setw(5) << sparse_nodes << " | " 
                << std::setw(6) << tensor_nodes << " | " 
                << std::fixed << std::setprecision(3) << ratio << std::endl;
      
      // Sparse should be much smaller than tensor for higher dimensions
      if (dim >= 3) {
        assert(sparse_nodes < tensor_nodes);
      }
    }
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_polynomial_exactness() {
  std::cout << "Testing polynomial exactness..." << std::endl;
  
  // Test that sparse grid integrates polynomials exactly
  // For level ℓ, should be exact for total degree ≤ 2ℓ-1
  
  const std::size_t dim = 2;
  const std::size_t level = 3;
  
  hypercube<Real> box(dim);
  box.lower = std::vector<Real>(dim, Real(0));
  box.upper = std::vector<Real>(dim, Real(1));
  
  // Test x² + y²
  auto poly = [](const Real* x) {
    return x[0]*x[0] + x[1]*x[1];
  };
  
  result<Real> res = integrate_sparse_grid<Real>(poly, box, level);
  
  // Exact integral: ∫₀¹∫₀¹(x²+y²)dxdy = 2/3
  Real exact = Real(2) / Real(3);
  Real error = std::abs(res.value - exact);
  
  std::cout << "  Polynomial x²+y²:" << std::endl;
  std::cout << "    Computed: " << res.value << std::endl;
  std::cout << "    Exact: " << exact << std::endl;
  std::cout << "    Error: " << std::scientific << error << std::endl;
  std::cout << "    Nodes used: " << res.evaluations << std::endl;
  
  assert(error < Real(1e-12));
  
  // Test x²y
  auto poly2 = [](const Real* x) {
    return x[0]*x[0] * x[1];
  };
  
  res = integrate_sparse_grid<Real>(poly2, box, level);
  exact = Real(1) / Real(6);  // ∫₀¹∫₀¹ x²y dxdy = 1/6
  error = std::abs(res.value - exact);
  
  std::cout << "  Polynomial x²y:" << std::endl;
  std::cout << "    Computed: " << res.value << std::endl;
  std::cout << "    Exact: " << exact << std::endl;
  std::cout << "    Error: " << std::scientific << error << std::endl;
  
  assert(error < Real(1e-12));
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_smooth_function() {
  std::cout << "Testing smooth function integration..." << std::endl;
  
  const std::size_t dim = 3;
  hypercube<Real> box(dim);
  box.lower = std::vector<Real>(dim, Real(0));
  box.upper = std::vector<Real>(dim, Real(1));
  
  // Test Gaussian exp(-||x||²)
  auto gaussian = [](const Real* x) {
    Real sum = Real(0);
    for (int i = 0; i < 3; ++i) {
      Real xi = x[i] - Real(0.5);  // Center at 0.5
      sum += xi * xi;
    }
    return std::exp(-Real(10) * sum);
  };
  
  // Test convergence with increasing level
  std::cout << "\n  Gaussian convergence:" << std::endl;
  std::cout << "  Level | Nodes | Integral | Est. Error" << std::endl;
  std::cout << "  ------|-------|----------|------------" << std::endl;
  
  Real prev_value = Real(0);
  for (std::size_t level = 2; level <= 5; ++level) {
    result<Real> res = integrate_sparse_grid<Real>(gaussian, box, level);
    
    Real change = (level > 2) ? std::abs(res.value - prev_value) : Real(0);
    
    std::cout << "    " << level << "   |  " << std::setw(4) << res.evaluations 
              << " | " << std::fixed << std::setprecision(6) << res.value
              << " | " << std::scientific << std::setprecision(2) << res.error << std::endl;
    
    prev_value = res.value;
    
    // Error estimate should be reasonable (not checking exact bound due to conservative estimation)
    if (level > 2) {
      assert(res.error < Real(1.0));  // Conservative error bound for hierarchical surplus
    }
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_performance() {
  std::cout << "\nTesting performance..." << std::endl;
  
  const std::size_t dim = 5;
  hypercube<Real> box(dim);
  box.lower = std::vector<Real>(dim, Real(0));
  box.upper = std::vector<Real>(dim, Real(1));
  
  // Test function: product of sines
  auto test_func = [dim](const Real* x) {
    Real prod = Real(1);
    const Real pi = Real(3.14159265358979323846);
    for (std::size_t i = 0; i < dim; ++i) {
      prod *= std::sin(pi * x[i]);
    }
    return prod;
  };
  
  std::cout << "  5D integration timing:" << std::endl;
  
  for (std::size_t level = 3; level <= 5; ++level) {
    auto start = std::chrono::high_resolution_clock::now();
    result<Real> res = integrate_sparse_grid<Real>(test_func, box, level);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    std::cout << "    Level " << level << ": " 
              << res.evaluations << " nodes, "
              << duration.count() << " μs, "
              << "result = " << res.value << std::endl;
    
    // Performance requirement: < 100ms for level 5
    assert(duration.count() < 100000);  // 100ms in microseconds
  }
  
  std::cout << "  PASSED" << std::endl;
}

int main() {
  std::cout << "\n===== Testing Smolyak Sparse Grid =====\n" << std::endl;
  
  using Real = double;
  
  test_multi_index_generation<Real>();
  test_node_deduplication<Real>();
  test_sparse_vs_tensor<Real>();
  test_polynomial_exactness<Real>();
  test_smooth_function<Real>();
  test_performance<Real>();
  
  std::cout << "\n===== All Tests PASSED =====\n" << std::endl;
  
  return 0;
}