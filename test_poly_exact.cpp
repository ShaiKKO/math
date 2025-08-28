#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>

// Include minimal headers directly
#include "./include/boost/math/cubature/regions.hpp"
#include "./include/boost/math/cubature/policies.hpp" 
#include "./include/boost/math/cubature/detail/cc_rules.hpp"
#include "./include/boost/math/cubature/detail/sparse_grid_impl.hpp"

using Real = double;
using namespace boost::math::cubature::detail;

// Test polynomial exactness for degree d
void test_polynomial_exactness(std::size_t dim, std::size_t level) {
  std::cout << "\n=== Testing Polynomial Exactness ===" << std::endl;
  std::cout << "Dimension: " << dim << ", Level: " << level << std::endl;
  
  // Create box [0,1]^dim
  boost::math::cubature::hypercube<Real> box(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    box.lower[i] = 0.0;
    box.upper[i] = 1.0;
  }
  
  // Clenshaw-Curtis rules should integrate polynomials exactly
  // Level l integrates polynomials up to degree 2*n-1 where n = num_points(l)
  // For slow growth: level 0: 1 point (degree 0), level 1: 3 points (degree 2), 
  //                  level 2: 5 points (degree 4), level 3: 9 points (degree 8)
  
  // Test 1: Constant function (degree 0)
  {
    auto f = [](const Real* x) { return Real(1); };
    smolyak_grid<Real> grid(dim, level);
    auto result = grid.evaluate(f, box);
    Real exact = 1.0;
    Real error = std::abs(result.value - exact);
    std::cout << "  Constant: result=" << result.value << ", exact=" << exact 
              << ", error=" << error << " (nodes=" << result.evaluations << ")" << std::endl;
    assert(error < 1e-14);
  }
  
  // Test 2: Linear function (degree 1)
  if (dim <= 3) {
    auto f = [dim](const Real* x) { 
      Real sum = 0;
      for (std::size_t i = 0; i < dim; ++i) sum += x[i];
      return sum;
    };
    smolyak_grid<Real> grid(dim, level);
    auto result = grid.evaluate(f, box);
    Real exact = dim * 0.5;  // ∫₀¹ x dx = 1/2 for each dimension
    Real error = std::abs(result.value - exact);
    std::cout << "  Linear: result=" << result.value << ", exact=" << exact 
              << ", error=" << error << std::endl;
    assert(error < 1e-14 || error < 1e-10 * std::abs(exact));
  }
  
  // Test 3: Quadratic function (degree 2)
  if (dim <= 3 && level >= 1) {
    auto f = [dim](const Real* x) { 
      Real sum = 0;
      for (std::size_t i = 0; i < dim; ++i) sum += x[i] * x[i];
      return sum;
    };
    smolyak_grid<Real> grid(dim, level);
    auto result = grid.evaluate(f, box);
    Real exact = dim / 3.0;  // ∫₀¹ x² dx = 1/3 for each dimension
    Real error = std::abs(result.value - exact);
    std::cout << "  Quadratic: result=" << result.value << ", exact=" << exact 
              << ", error=" << error << std::endl;
    assert(error < 1e-14 || error < 1e-10 * std::abs(exact));
  }
  
  // Test 4: Mixed degree polynomial x₀*x₁ (degree 2)
  if (dim >= 2 && level >= 1) {
    auto f = [](const Real* x) { 
      return x[0] * x[1];
    };
    smolyak_grid<Real> grid(dim, level);
    auto result = grid.evaluate(f, box);
    Real exact = 0.25;  // ∫₀¹∫₀¹ x*y dxdy = 1/4
    Real error = std::abs(result.value - exact);
    std::cout << "  Mixed x₀*x₁: result=" << result.value << ", exact=" << exact 
              << ", error=" << error << std::endl;
    assert(error < 1e-14 || error < 1e-10 * std::abs(exact));
  }
  
  // Test 5: Higher degree polynomial x³ (degree 3) - should be exact for level >= 2
  if (dim == 1 && level >= 2) {
    auto f = [](const Real* x) { 
      return x[0] * x[0] * x[0];
    };
    smolyak_grid<Real> grid(dim, level);
    auto result = grid.evaluate(f, box);
    Real exact = 0.25;  // ∫₀¹ x³ dx = 1/4
    Real error = std::abs(result.value - exact);
    std::cout << "  Cubic (1D): result=" << result.value << ", exact=" << exact 
              << ", error=" << error << std::endl;
    if (level >= 2) {
      assert(error < 1e-14 || error < 1e-10 * std::abs(exact));
    }
  }
}

int main() {
  std::cout << std::fixed << std::setprecision(15);
  
  std::cout << "\n===== Testing Smolyak Sparse Grid Polynomial Exactness =====\n";
  
  // Test different dimensions and levels
  test_polynomial_exactness(1, 1);
  test_polynomial_exactness(1, 2);
  test_polynomial_exactness(1, 3);
  
  test_polynomial_exactness(2, 1);
  test_polynomial_exactness(2, 2);
  test_polynomial_exactness(2, 3);
  
  test_polynomial_exactness(3, 1);
  test_polynomial_exactness(3, 2);
  
  std::cout << "\n===== All Polynomial Exactness Tests PASSED =====\n" << std::endl;
  
  return 0;
}