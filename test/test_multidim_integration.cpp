// Test multi-dimensional integration with fixed GM rules

#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace boost::math::cubature;
using namespace boost::math::constants;

// Test integrand: product of cosines
template <typename Real, std::size_t Dim>
Real product_cosine(const Real* x, std::size_t) {
  Real prod = 1;
  for (std::size_t i = 0; i < Dim; ++i) {
    prod *= std::cos(pi<Real>() * x[i] / 2);
  }
  return prod;
}

// Test integrand: exponential of negative sum of squares
template <typename Real, std::size_t Dim>
Real gaussian(const Real* x, std::size_t) {
  Real sum = 0;
  for (std::size_t i = 0; i < Dim; ++i) {
    Real xi = 2*x[i] - 1; // Map from [0,1] to [-1,1]
    sum += xi * xi;
  }
  return std::exp(-sum);
}

template <typename Real, std::size_t Dim>
void test_dimension() {
  std::cout << "\n=== Testing " << Dim << "D integration ===" << std::endl;
  
  hypercube<Real> box(Dim);
  for (std::size_t i = 0; i < Dim; ++i) {
    box.lower[i] = 0;
    box.upper[i] = 1;
  }
  
  // Test 1: Product of cosines
  {
    auto f = product_cosine<Real, Dim>;
    auto result = integrate_adaptive<Real>(f, box, 1e-6, 1e-6, 100000);
    
    // Exact value: (2/π)^Dim
    Real exact = std::pow(2/pi<Real>(), Dim);
    Real error = std::abs(result.value - exact);
    
    std::cout << "Product of cosines:" << std::endl;
    std::cout << "  Result: " << std::setprecision(10) << result.value << std::endl;
    std::cout << "  Exact:  " << exact << std::endl;
    std::cout << "  Error:  " << std::scientific << error << std::endl;
    std::cout << "  Est err:" << result.error << std::endl;
    std::cout << "  Evals:  " << result.evaluations << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.status);
    switch(result.status) {
      case status_code::success: std::cout << " (success)"; break;
      case status_code::maxeval_reached: std::cout << " (maxeval)"; break;
      default: std::cout << " (other)";
    }
    std::cout << std::endl;
  }
  
  // Test 2: Gaussian (mapped to [-1,1]^d)
  {
    auto f = gaussian<Real, Dim>;
    auto result = integrate_adaptive<Real>(f, box, 1e-6, 1e-6, 100000);
    
    // Exact value: (√π * erf(1))^Dim / 2^Dim
    Real erf1 = std::erf(Real(1));
    Real exact = std::pow(root_pi<Real>() * erf1 / 2, Dim);
    Real error = std::abs(result.value - exact);
    
    std::cout << "Gaussian exp(-||x||²):" << std::endl;
    std::cout << "  Result: " << std::setprecision(10) << result.value << std::endl;
    std::cout << "  Exact:  " << exact << std::endl;
    std::cout << "  Error:  " << std::scientific << error << std::endl;
    std::cout << "  Est err:" << result.error << std::endl;
    std::cout << "  Evals:  " << result.evaluations << std::endl;
    std::cout << "  Status: " << static_cast<int>(result.status);
    switch(result.status) {
      case status_code::success: std::cout << " (success)"; break;
      case status_code::maxeval_reached: std::cout << " (maxeval)"; break;
      default: std::cout << " (other)";
    }
    std::cout << std::endl;
  }
}

// Test vector integration in multiple dimensions
template <typename Real, std::size_t Dim>
void test_vector_dimension() {
  std::cout << "\n=== Testing " << Dim << "D vector integration ===" << std::endl;
  
  // Vector integrand: [cos(πx₀/2), sin(πx₁/2), exp(-||x||²)]
  auto f_vec = [](const Real* x, Real* out, std::size_t m) {
    if (m >= 1) out[0] = std::cos(pi<Real>() * x[0] / 2);
    if (m >= 2) out[1] = (Dim >= 2) ? std::sin(pi<Real>() * x[1] / 2) : Real(1);
    if (m >= 3) {
      Real sum = 0;
      for (std::size_t i = 0; i < Dim; ++i) {
        Real xi = 2*x[i] - 1;
        sum += xi * xi;
      }
      out[2] = std::exp(-sum);
    }
  };
  
  hypercube<Real> box(Dim);
  for (std::size_t i = 0; i < Dim; ++i) {
    box.lower[i] = 0;
    box.upper[i] = 1;
  }
  
  auto results = integrate_adaptive_vector<Real>(f_vec, box, 3, 1e-6, 1e-6, 100000);
  
  std::cout << "Vector components:" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    std::cout << "  Component " << i << ": value=" << std::setprecision(6) 
              << results[i].value << ", error=" << std::scientific 
              << results[i].error << ", evals=" << results[i].evaluations << std::endl;
  }
}

int main() {
  std::cout << "===== Multi-dimensional Integration Tests =====" << std::endl;
  std::cout << "Testing integration in various dimensions\n" << std::endl;
  
  using Real = double;
  
  // Test scalar integration
  test_dimension<Real, 2>();
  test_dimension<Real, 3>();
  test_dimension<Real, 4>();
  test_dimension<Real, 5>();
  test_dimension<Real, 6>();
  
  // Test vector integration
  test_vector_dimension<Real, 2>();
  test_vector_dimension<Real, 3>();
  
  std::cout << "\n===== Tests Complete =====" << std::endl;
  
  return 0;
}