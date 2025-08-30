// Test sparse grid vector integration

#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace boost::math::cubature;
using namespace boost::math::constants;

// Vector test integrand: [cos(πx₀/2), sin(πx₁/2), exp(-||x||²)]
template <typename Real>
void vector_integrand(const Real* x, Real* out, std::size_t num_components) {
  std::size_t dim = 3; // We'll test in 3D
  
  if (num_components >= 1) {
    // Component 0: cos(π*x₀/2)
    out[0] = std::cos(pi<Real>() * x[0] / 2);
  }
  
  if (num_components >= 2) {
    // Component 1: sin(π*x₁/2)
    out[1] = std::sin(pi<Real>() * x[1] / 2);
  }
  
  if (num_components >= 3) {
    // Component 2: exp(-||x||²) where x is mapped from [0,1] to [-1,1]
    Real sum = 0;
    for (std::size_t i = 0; i < dim; ++i) {
      Real xi = 2*x[i] - 1;
      sum += xi * xi;
    }
    out[2] = std::exp(-sum);
  }
}

// Another vector integrand for testing: polynomial components
template <typename Real>
void polynomial_vector(const Real* x, Real* out, std::size_t num_components) {
  if (num_components >= 1) out[0] = x[0] * x[0];           // x₀²
  if (num_components >= 2) out[1] = x[0] * x[1];           // x₀x₁
  if (num_components >= 3) out[2] = x[0] + x[1] + x[2];    // x₀ + x₁ + x₂
  if (num_components >= 4) out[3] = Real(1);               // constant
}

template <typename Real>
void test_sparse_grid_vector(std::size_t dim, std::size_t level) {
  std::cout << "\n=== Testing Sparse Grid Vector Integration ===" << std::endl;
  std::cout << "Dimension: " << dim << ", Level: " << level << std::endl;
  
  // Create integration domain [0,1]^dim
  hypercube<Real> box(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    box.lower[i] = 0;
    box.upper[i] = 1;
  }
  
  // Test vector integrand with 3 components
  {
    std::cout << "\nVector integrand test (3 components):" << std::endl;
    auto results = integrate_sparse_grid_vector<Real>(
        vector_integrand<Real>, box, 3, level);
    
    // Expected values (computed analytically or with high precision)
    // Component 0: integral of cos(πx/2) from 0 to 1 = 2/π ≈ 0.6366197724
    // Component 1: integral of sin(πx/2) from 0 to 1 = 2/π ≈ 0.6366197724
    // Component 2: integral of exp(-||x||²) over [-1,1]³ mapped to [0,1]³
    Real exact0 = 2 / pi<Real>();
    Real exact1 = 2 / pi<Real>();
    Real erf1 = std::erf(Real(1));
    Real exact2 = std::pow(root_pi<Real>() * erf1 / 2, dim);
    
    std::cout << "Component 0 (cos):" << std::endl;
    std::cout << "  Result: " << std::setprecision(10) << results[0].value << std::endl;
    std::cout << "  Exact:  " << exact0 << std::endl;
    std::cout << "  Error:  " << std::scientific << std::abs(results[0].value - exact0) << std::endl;
    std::cout << "  Est err:" << results[0].error << std::endl;
    
    std::cout << "Component 1 (sin):" << std::endl;
    std::cout << "  Result: " << std::setprecision(10) << results[1].value << std::endl;
    std::cout << "  Exact:  " << exact1 << std::endl;
    std::cout << "  Error:  " << std::scientific << std::abs(results[1].value - exact1) << std::endl;
    std::cout << "  Est err:" << results[1].error << std::endl;
    
    std::cout << "Component 2 (gaussian):" << std::endl;
    std::cout << "  Result: " << std::setprecision(10) << results[2].value << std::endl;
    std::cout << "  Exact:  " << exact2 << std::endl;
    std::cout << "  Error:  " << std::scientific << std::abs(results[2].value - exact2) << std::endl;
    std::cout << "  Est err:" << results[2].error << std::endl;
    
    std::cout << "  Evaluations: " << results[0].evaluations << std::endl;
  }
  
  // Test polynomial vector with 4 components
  {
    std::cout << "\nPolynomial vector test (4 components):" << std::endl;
    auto results = integrate_sparse_grid_vector<Real>(
        polynomial_vector<Real>, box, 4, level);
    
    // Exact values for polynomials over [0,1]³
    // Component 0: ∫x₀² = 1/3
    // Component 1: ∫x₀x₁ = 1/4
    // Component 2: ∫(x₀ + x₁ + x₂) = 3/2
    // Component 3: ∫1 = 1
    Real exact[4] = {Real(1)/3, Real(1)/4, Real(3)/2, Real(1)};
    
    for (std::size_t i = 0; i < 4; ++i) {
      std::cout << "Component " << i << ":" << std::endl;
      std::cout << "  Result: " << std::setprecision(10) << results[i].value << std::endl;
      std::cout << "  Exact:  " << exact[i] << std::endl;
      std::cout << "  Error:  " << std::scientific 
                << std::abs(results[i].value - exact[i]) << std::endl;
    }
    std::cout << "  Evaluations: " << results[0].evaluations << std::endl;
  }
}

template <typename Real>
void test_different_levels() {
  std::cout << "\n=== Testing Different Sparse Grid Levels ===" << std::endl;
  
  hypercube<Real> box(2);  // 2D for simplicity
  box.lower[0] = box.lower[1] = 0;
  box.upper[0] = box.upper[1] = 1;
  
  // Simple vector integrand: [x₀, x₁]
  auto f = [](const Real* x, Real* out, std::size_t m) {
    if (m >= 1) out[0] = x[0];
    if (m >= 2) out[1] = x[1];
  };
  
  std::cout << "\nExact values: [0.5, 0.5]" << std::endl;
  
  for (std::size_t level = 0; level <= 4; ++level) {
    auto results = integrate_sparse_grid_vector<Real>(f, box, 2, level);
    
    std::cout << "\nLevel " << level << ":" << std::endl;
    std::cout << "  Results: [" << std::setprecision(10) 
              << results[0].value << ", " << results[1].value << "]" << std::endl;
    std::cout << "  Errors:  [" << std::scientific 
              << std::abs(results[0].value - 0.5) << ", " 
              << std::abs(results[1].value - 0.5) << "]" << std::endl;
    std::cout << "  Evaluations: " << results[0].evaluations << std::endl;
  }
}

int main() {
  std::cout << "===== Sparse Grid Vector Integration Tests =====\n" << std::endl;
  
  using Real = double;
  
  // Test 3D integration at different levels
  test_sparse_grid_vector<Real>(3, 2);
  test_sparse_grid_vector<Real>(3, 3);
  
  // Test convergence with increasing levels
  test_different_levels<Real>();
  
  std::cout << "\n===== Tests Complete =====\n" << std::endl;
  
  return 0;
}