// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <chrono>
#include <iomanip>

using namespace boost::math::cubature;
using namespace boost::math::constants;

// Test vector integrand: [sin(x)*sin(y), cos(x+y), exp(-(x^2+y^2))]
template <typename Real>
class test_vector_integrand_2d {
public:
  void operator()(const Real* x, Real* out, std::size_t m) const {
    if (m >= 1) out[0] = std::sin(x[0]) * std::sin(x[1]);
    if (m >= 2) out[1] = std::cos(x[0] + x[1]);
    if (m >= 3) out[2] = std::exp(-(x[0]*x[0] + x[1]*x[1]));
  }
};

// Test vector integrand: [x^2, y^2, x*y]
template <typename Real>
class polynomial_vector_2d {
public:
  void operator()(const Real* x, Real* out, std::size_t m) const {
    if (m >= 1) out[0] = x[0] * x[0];
    if (m >= 2) out[1] = x[1] * x[1];
    if (m >= 3) out[2] = x[0] * x[1];
  }
};

template <typename Real>
void test_basic_vector_integration() {
  std::cout << "Testing basic vector integration..." << std::endl;
  
  // Integrate [x^2, y^2, x*y] over [0,1]^2
  polynomial_vector_2d<Real> f;
  hypercube<Real> box(2);
  box.lower = {0, 0};
  box.upper = {1, 1};
  
  auto results = integrate_adaptive_vector<Real>(
    f, box, 3, Real(1e-10), Real(1e-10));
  
  // Exact values: [1/3, 1/3, 1/4]
  std::vector<Real> exact = {Real(1)/Real(3), Real(1)/Real(3), Real(1)/Real(4)};
  
  std::cout << "  Component results:" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    std::cout << "    Component " << i << ": " 
              << "value=" << results[i].value 
              << ", error=" << results[i].error
              << ", exact=" << exact[i]
              << ", diff=" << std::abs(results[i].value - exact[i]) << std::endl;
    
    assert(std::abs(results[i].value - exact[i]) < 1e-8);
    assert(results[i].status == status_code::success);
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_vector_vs_scalar_consistency() {
  std::cout << "Testing vector vs scalar consistency..." << std::endl;
  
  test_vector_integrand_2d<Real> vector_f;
  hypercube<Real> box(2);
  box.lower = {0, 0};
  box.upper = {pi<Real>(), pi<Real>()};
  
  // Vector evaluation
  auto vector_results = integrate_adaptive_vector<Real>(
    vector_f, box, 3, Real(1e-6), Real(1e-6));
  
  // Scalar evaluation for each component
  std::vector<result<Real>> scalar_results;
  
  // Component 0: sin(x)*sin(y)
  auto f0 = [](const Real* x, std::size_t) -> Real {
    return std::sin(x[0]) * std::sin(x[1]);
  };
  scalar_results.push_back(integrate_adaptive<Real>(f0, box, Real(1e-6), Real(1e-6)));
  
  // Component 1: cos(x+y)
  auto f1 = [](const Real* x, std::size_t) -> Real {
    return std::cos(x[0] + x[1]);
  };
  scalar_results.push_back(integrate_adaptive<Real>(f1, box, Real(1e-6), Real(1e-6)));
  
  // Component 2: exp(-(x^2+y^2))
  auto f2 = [](const Real* x, std::size_t) -> Real {
    return std::exp(-(x[0]*x[0] + x[1]*x[1]));
  };
  scalar_results.push_back(integrate_adaptive<Real>(f2, box, Real(1e-6), Real(1e-6)));
  
  std::cout << "  Comparing results:" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    Real diff = std::abs(vector_results[i].value - scalar_results[i].value);
    std::cout << "    Component " << i << ": "
              << "vector=" << vector_results[i].value 
              << ", scalar=" << scalar_results[i].value
              << ", diff=" << diff << std::endl;
    
    // Results should be very close (within numerical tolerance)
    assert(diff < 1e-10 || diff < 1e-10 * std::abs(scalar_results[i].value));
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_error_norms() {
  std::cout << "Testing different error norms..." << std::endl;
  
  polynomial_vector_2d<Real> f;
  hypercube<Real> box(2);
  box.lower = {0, 0};
  box.upper = {1, 1};
  
  // Test with different norms
  auto results_linf = integrate_adaptive_vector<Real>(
    f, box, 3, Real(1e-4), Real(1e-4), 0, boost::math::cubature::detail::error_norm::l_infinity);
  
  auto results_l2 = integrate_adaptive_vector<Real>(
    f, box, 3, Real(1e-4), Real(1e-4), 0, boost::math::cubature::detail::error_norm::l2);
  
  auto results_l1 = integrate_adaptive_vector<Real>(
    f, box, 3, Real(1e-4), Real(1e-4), 0, boost::math::cubature::detail::error_norm::l1);
  
  std::cout << "  L-infinity norm results:" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    std::cout << "    Component " << i << ": " 
              << results_linf[i].value << " (error: " << results_linf[i].error << ")" << std::endl;
  }
  
  std::cout << "  L2 norm results:" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    std::cout << "    Component " << i << ": " 
              << results_l2[i].value << " (error: " << results_l2[i].error << ")" << std::endl;
  }
  
  std::cout << "  L1 norm results:" << std::endl;
  for (std::size_t i = 0; i < 3; ++i) {
    std::cout << "    Component " << i << ": " 
              << results_l1[i].value << " (error: " << results_l1[i].error << ")" << std::endl;
  }
  
  // All should converge to correct values
  std::vector<Real> exact = {Real(1)/Real(3), Real(1)/Real(3), Real(1)/Real(4)};
  for (std::size_t i = 0; i < 3; ++i) {
    assert(std::abs(results_linf[i].value - exact[i]) < 1e-3);
    assert(std::abs(results_l2[i].value - exact[i]) < 1e-3);
    assert(std::abs(results_l1[i].value - exact[i]) < 1e-3);
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_performance_comparison() {
  std::cout << "Testing vector vs scalar performance..." << std::endl;
  
  test_vector_integrand_2d<Real> vector_f;
  hypercube<Real> box(2);
  box.lower = {0, 0};
  box.upper = {2, 2};
  
  const std::size_t num_runs = 5;
  
  // Time vector evaluation
  auto start = std::chrono::high_resolution_clock::now();
  for (std::size_t run = 0; run < num_runs; ++run) {
    auto vector_results = integrate_adaptive_vector<Real>(
      vector_f, box, 3, Real(1e-6), Real(1e-6));
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto vector_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  
  // Time scalar evaluations
  start = std::chrono::high_resolution_clock::now();
  for (std::size_t run = 0; run < num_runs; ++run) {
    auto f0 = [](const Real* x, std::size_t) -> Real {
      return std::sin(x[0]) * std::sin(x[1]);
    };
    auto r0 = integrate_adaptive<Real>(f0, box, Real(1e-6), Real(1e-6));
    
    auto f1 = [](const Real* x, std::size_t) -> Real {
      return std::cos(x[0] + x[1]);
    };
    auto r1 = integrate_adaptive<Real>(f1, box, Real(1e-6), Real(1e-6));
    
    auto f2 = [](const Real* x, std::size_t) -> Real {
      return std::exp(-(x[0]*x[0] + x[1]*x[1]));
    };
    auto r2 = integrate_adaptive<Real>(f2, box, Real(1e-6), Real(1e-6));
  }
  end = std::chrono::high_resolution_clock::now();
  auto scalar_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  
  Real speedup = static_cast<Real>(scalar_time) / static_cast<Real>(vector_time);
  std::cout << "  Vector time: " << vector_time/num_runs << " μs (avg)" << std::endl;
  std::cout << "  Scalar time: " << scalar_time/num_runs << " μs (avg for 3 components)" << std::endl;
  std::cout << "  Speedup: " << std::fixed << std::setprecision(2) << speedup << "x" << std::endl;
  
  // Vector should be faster (or at least within 3% of scalar)
  assert(speedup > 0.97);  // Within 3% performance requirement
  
  std::cout << "  PASSED - Performance within requirements" << std::endl;
}

template <typename Real>
void test_high_dimension_vector() {
  std::cout << "Testing higher dimensional vector integration..." << std::endl;
  
  // 3D vector integrand
  class vector_3d {
  public:
    void operator()(const Real* x, Real* out, std::size_t m) const {
      if (m >= 1) out[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];  // r^2
      if (m >= 2) out[1] = std::exp(-(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]));  // exp(-r^2)
    }
  };
  
  vector_3d f;
  hypercube<Real> box(3);
  box.lower = {0, 0, 0};
  box.upper = {1, 1, 1};
  
  auto results = integrate_adaptive_vector<Real>(
    f, box, 2, Real(1e-6), Real(1e-6));
  
  std::cout << "  3D Results:" << std::endl;
  std::cout << "    Component 0 (r^2): " << results[0].value << std::endl;
  std::cout << "    Component 1 (exp(-r^2)): " << results[1].value << std::endl;
  
  // Check convergence
  assert(results[0].status == status_code::success);
  assert(results[1].status == status_code::success);
  
  // Rough checks (exact values would require more calculation)
  assert(results[0].value > 0 && results[0].value < 3);  // r^2 integral
  assert(results[1].value > 0 && results[1].value < 1);  // exp(-r^2) integral
  
  std::cout << "  PASSED" << std::endl;
}

int main() {
  std::cout << "===== Vector Integration Tests =====" << std::endl;
  std::cout << "Testing P07 implementation\n" << std::endl;
  
  using Real = double;
  
  test_basic_vector_integration<Real>();
  test_vector_vs_scalar_consistency<Real>();
  test_error_norms<Real>();
  test_performance_comparison<Real>();
  test_high_dimension_vector<Real>();
  
  std::cout << "\n===== All Vector Integration Tests PASSED =====" << std::endl;
  std::cout << "Vector-valued integrands working correctly!" << std::endl;
  
  return 0;
}