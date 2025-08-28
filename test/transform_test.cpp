#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <boost/math/cubature/transforms.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/constants/constants.hpp>

using namespace boost::math::cubature;
using namespace boost::math::constants;

template <typename Real>
void test_rational_transform() {
  std::cout << "Testing rational transform..." << std::endl;
  
  rational_transform<Real> transform;
  
  // Test forward transform
  auto [x1, j1] = transform.forward(Real(0));
  assert(x1 == Real(0));
  assert(j1 == Real(1));
  
  auto [x2, j2] = transform.forward(Real(0.5));
  assert(std::abs(x2 - Real(1)) < 1e-10);
  assert(std::abs(j2 - Real(4)) < 1e-10);
  
  auto [x3, j3] = transform.forward(Real(0.75));
  assert(std::abs(x3 - Real(3)) < 1e-10);
  assert(std::abs(j3 - Real(16)) < 1e-10);
  
  // Test inverse transform
  assert(std::abs(transform.inverse(Real(0)) - Real(0)) < 1e-10);
  assert(std::abs(transform.inverse(Real(1)) - Real(0.5)) < 1e-10);
  assert(std::abs(transform.inverse(Real(3)) - Real(0.75)) < 1e-10);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_tangent_transform() {
  std::cout << "Testing tangent transform..." << std::endl;
  
  tangent_transform<Real> transform;
  
  // Test at u = 0.5 (should give x = 0)
  auto [x1, j1] = transform.forward(Real(0.5));
  assert(std::abs(x1) < 1e-10);
  assert(std::abs(j1 - pi<Real>()) < 1e-10);
  
  // Test inverse at x = 0
  assert(std::abs(transform.inverse(Real(0)) - Real(0.5)) < 1e-10);
  
  // Test symmetry
  auto [x2, j2] = transform.forward(Real(0.25));
  auto [x3, j3] = transform.forward(Real(0.75));
  assert(std::abs(x2 + x3) < 1e-10);  // Should be symmetric
  assert(std::abs(j2 - j3) < 1e-10);   // Jacobians should match
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_exponential_transform() {
  std::cout << "Testing exponential transform..." << std::endl;
  
  exponential_transform<Real> transform;
  
  // Test at u = 0
  auto [x1, j1] = transform.forward(Real(0));
  assert(x1 == Real(0));
  assert(j1 == Real(1));
  
  // Test at u = 1 - 1/e
  Real u = Real(1) - Real(1) / std::exp(Real(1));
  auto [x2, j2] = transform.forward(u);
  assert(std::abs(x2 - Real(1)) < 1e-10);
  assert(std::abs(j2 - std::exp(Real(1))) < 1e-10);
  
  // Test inverse
  assert(std::abs(transform.inverse(Real(0)) - Real(0)) < 1e-10);
  assert(std::abs(transform.inverse(Real(1)) - u) < 1e-10);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_duffy_transform() {
  std::cout << "Testing Duffy transform..." << std::endl;
  
  // Test 2D case
  {
    Real u[2] = {0.5, 0.5};
    Real x[2];
    Real jac = duffy_transform<Real>::apply(u, x, 2);
    
    // 2D: (u,v) -> (u(1-v), uv) with Jacobian = u
    assert(std::abs(x[0] - Real(0.25)) < 1e-10);   // u[0]*(1-u[1]) = 0.5*0.5
    assert(std::abs(x[1] - Real(0.25)) < 1e-10);   // u[0]*u[1] = 0.5*0.5
    assert(std::abs(jac - Real(0.5)) < 1e-10);     // u[0] = 0.5
  }
  
  // Test 3D case
  {
    Real u[3] = {0.5, 0.5, 0.5};
    Real x[3];
    Real jac = duffy_transform<Real>::apply(u, x, 3);
    
    // 3D: (u,v,w) -> (u(1-v), uv(1-w), uvw) with Jacobian = u²v
    assert(std::abs(x[0] - Real(0.25)) < 1e-10);    // u[0]*(1-u[1]) = 0.5*0.5
    assert(std::abs(x[1] - Real(0.125)) < 1e-10);   // u[0]*u[1]*(1-u[2]) = 0.5*0.5*0.5  
    assert(std::abs(x[2] - Real(0.125)) < 1e-10);   // u[0]*u[1]*u[2] = 0.5*0.5*0.5
    assert(std::abs(jac - Real(0.125)) < 1e-10);    // u[0]²*u[1] = 0.25*0.5
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_infinite_integral_gaussian() {
  std::cout << "Testing infinite integral (Gaussian)..." << std::endl;
  
  // Integral of exp(-x^2) from 0 to infinity = sqrt(pi)/2
  // Using rational transform: u/(1-u) maps [0,1] to [0,∞)
  rational_transform<Real> transform;
  
  // Create transformed integrand for 1D quadrature
  auto integrand = [&transform](Real u) -> Real {
    if (u >= Real(1)) return Real(0);  // Handle boundary
    auto [x, jac] = transform.forward(u);
    return std::exp(-x * x) * jac;
  };
  
  // Use tanh_sinh quadrature for [0,1] interval
  boost::math::quadrature::tanh_sinh<Real> quad;
  Real result = quad.integrate(integrand, Real(0), Real(1));
  
  Real exact = std::sqrt(pi<Real>()) / Real(2);
  Real error = std::abs(result - exact);
  
  std::cout << "  Computed: " << result << std::endl;
  std::cout << "  Exact: " << exact << std::endl;
  std::cout << "  Error: " << error << std::endl;
  assert(error < 1e-6);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_bi_infinite_integral() {
  std::cout << "Testing bi-infinite integral (sech)..." << std::endl;
  
  // Integral of sech(x) from -infinity to infinity = pi
  // Using tangent transform: tan(π(u-1/2)) maps [0,1] to (-∞,∞)
  tangent_transform<Real> transform;
  
  // Create transformed integrand for 1D quadrature
  auto integrand = [&transform](Real u) -> Real {
    auto [x, jac] = transform.forward(u);
    return (Real(1) / std::cosh(x)) * jac;
  };
  
  // Use tanh_sinh quadrature for [0,1] interval
  boost::math::quadrature::tanh_sinh<Real> quad;
  Real result = quad.integrate(integrand, Real(0), Real(1));
  
  Real exact = pi<Real>();
  Real error = std::abs(result - exact);
  
  std::cout << "  Computed: " << result << std::endl;
  std::cout << "  Exact: " << exact << std::endl;
  std::cout << "  Error: " << error << std::endl;
  assert(error < 1e-4);
  
  std::cout << "  PASSED" << std::endl;
}

int main() {
  std::cout << "===== Transform Tests =====" << std::endl;
  
  using Real = double;
  
  test_rational_transform<Real>();
  test_tangent_transform<Real>();
  test_exponential_transform<Real>();
  test_duffy_transform<Real>();
  test_infinite_integral_gaussian<Real>();
  test_bi_infinite_integral<Real>();
  
  std::cout << "\n===== All Tests PASSED =====" << std::endl;
  
  return 0;
}