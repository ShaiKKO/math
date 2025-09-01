// Test to verify that Genz-Malik rules are properly connected for all dimensions
// Copyright 2025 Boost
// Distributed under the Boost Software License, Version 1.0.

#include <boost/math/cubature/detail/genz_malik_evaluator.hpp>
#include <boost/math/cubature/detail/adaptivity.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace boost::math::cubature;

// Simple test integrand: f(x) = product of all coordinates
template <typename Real>
struct product_integrand {
  Real operator()(const Real* x, std::size_t dim) const {
    Real prod = 1;
    for (std::size_t i = 0; i < dim; ++i) {
      prod *= x[i];
    }
    return prod;
  }
};

// Test that rules are available and give reasonable results
template <std::size_t Dim>
void test_dimension() {
  using Real = double;
  
  std::cout << "\nTesting dimension " << Dim << ":\n";
  
  // Create unit hypercube
  detail::region<Real> reg;
  reg.a.resize(Dim);
  reg.b.resize(Dim);
  for (std::size_t i = 0; i < Dim; ++i) {
    reg.a[i] = 0;
    reg.b[i] = 1;
  }
  
  // Create integrand
  product_integrand<Real> f;
  
  // Evaluate using Genz-Malik evaluator
  detail::embedded_pair_result<Real> result;
  bool success = detail::genz_malik_evaluator<Real>::evaluate_embedded_pair<product_integrand<Real>, Dim>(
      f, reg, result);
  
  if (!success) {
    std::cout << "  ERROR: Rules not available for dimension " << Dim << "\n";
    return;
  }
  
  // The exact integral of x1*x2*...*xd over [0,1]^d is 1/2^d
  Real exact = std::pow(0.5, Dim);
  
  std::cout << "  Primary estimate (deg 9): " << std::scientific << std::setprecision(10) 
            << result.primary_estimate << "\n";
  std::cout << "  Embedded estimate (deg 7): " << result.embedded_estimate << "\n";
  std::cout << "  Error estimate: " << result.error_estimate << "\n";
  std::cout << "  Exact value: " << exact << "\n";
  
  Real primary_error = std::abs(result.primary_estimate - exact);
  Real embedded_error = std::abs(result.embedded_estimate - exact);
  
  std::cout << "  Primary error: " << primary_error << "\n";
  std::cout << "  Embedded error: " << embedded_error << "\n";
  std::cout << "  Function evaluations: " << result.function_evaluations << "\n";
  
  // Verify reasonable accuracy (should be exact for polynomial of low degree)
  // Product of d variables has degree d, so:
  // - Degree 7 rule should be exact for d <= 7
  // - Degree 9 rule should be exact for d <= 9
  
  if (Dim <= 7) {
    // Embedded rule (degree 7) should be exact within roundoff
    if (embedded_error > 1e-10) {
      std::cout << "  WARNING: Degree-7 rule not exact for degree-" << Dim 
                << " polynomial (error = " << embedded_error << ")\n";
    } else {
      std::cout << "  ✓ Degree-7 rule exact for degree-" << Dim << " polynomial\n";
    }
  }
  
  if (Dim <= 9) {
    // Primary rule (degree 9) should be exact within roundoff
    if (primary_error > 1e-10) {
      std::cout << "  WARNING: Degree-9 rule not exact for degree-" << Dim 
                << " polynomial (error = " << primary_error << ")\n";
    } else {
      std::cout << "  ✓ Degree-9 rule exact for degree-" << Dim << " polynomial\n";
    }
  }
  
  // Also test a constant function (should be exact for any rule)
  struct constant_integrand {
    Real operator()(const Real*, std::size_t) const { return 1.0; }
  };
  
  constant_integrand g;
  detail::embedded_pair_result<Real> const_result;
  success = detail::genz_malik_evaluator<Real>::evaluate_embedded_pair<constant_integrand, Dim>(
      g, reg, const_result);
  
  if (success) {
    Real const_exact = 1.0; // Integral of 1 over unit hypercube
    Real const_error = std::abs(const_result.primary_estimate - const_exact);
    if (const_error > 1e-12) {
      std::cout << "  ERROR: Constant function not integrated exactly (error = " 
                << const_error << ")\n";
    } else {
      std::cout << "  ✓ Constant function integrated exactly\n";
    }
  }
}

int main() {
  std::cout << "=== Testing Genz-Malik Rule Connectivity ===\n";
  std::cout << "Testing that rules are available and accurate for dimensions 2-15\n";
  
  // Test dimensions 2-15 (the range claimed to be supported)
  test_dimension<2>();
  test_dimension<3>();
  test_dimension<4>();
  test_dimension<5>();
  test_dimension<6>();
  test_dimension<7>();
  test_dimension<8>();
  test_dimension<9>();
  test_dimension<10>();
  test_dimension<11>();
  test_dimension<12>();
  test_dimension<13>();
  test_dimension<14>();
  test_dimension<15>();
  
  // Also test dimension 16 (should fail gracefully)
  // Note: Can't test at compile time due to template instantiation
  // std::cout << "\n--- Testing unsupported dimension ---\n";
  // test_dimension<16>();
  
  std::cout << "\n=== Test Complete ===\n";
  
  return 0;
}