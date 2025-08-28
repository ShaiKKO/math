#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <array>
#include <boost/math/cubature/simplex.hpp>
#include <boost/math/cubature/regions.hpp>

using namespace boost::math::cubature;

template <typename Real>
void test_triangle_constant() {
  std::cout << "Testing triangle integration (constant function)..." << std::endl;
  
  // Standard right triangle with vertices (0,0), (1,0), (0,1)
  std::array<std::array<Real, 2>, 3> vertices = {{
    {{0, 0}},
    {{1, 0}},
    {{0, 1}}
  }};
  
  // Constant function f(x,y) = 1
  auto f = [](const Real* x) -> Real {
    return Real(1);
  };
  
  auto result = integrate_triangle<Real>(f, vertices, 1e-8, 1e-8);
  
  // Area of triangle should be 0.5
  Real exact = Real(0.5);
  Real error = std::abs(result.value - exact);
  
  std::cout << "  Computed: " << result.value << std::endl;
  std::cout << "  Exact: " << exact << std::endl;
  std::cout << "  Error: " << error << std::endl;
  assert(error < 1e-10);
  assert(result.status == status_code::success || result.status == status_code::maxeval_reached);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_triangle_linear() {
  std::cout << "Testing triangle integration (linear function)..." << std::endl;
  
  // Standard right triangle
  std::array<std::array<Real, 2>, 3> vertices = {{
    {{0, 0}},
    {{1, 0}},
    {{0, 1}}
  }};
  
  // Linear function f(x,y) = x + y
  auto f = [](const Real* x) -> Real {
    return x[0] + x[1];
  };
  
  auto result = integrate_triangle<Real>(f, vertices, 1e-8, 1e-8);
  
  // Integral of (x + y) over standard triangle = 1/3
  Real exact = Real(1) / Real(3);
  Real error = std::abs(result.value - exact);
  
  std::cout << "  Computed: " << result.value << std::endl;
  std::cout << "  Exact: " << exact << std::endl;
  std::cout << "  Error: " << error << std::endl;
  assert(error < 1e-8);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_triangle_polynomial() {
  std::cout << "Testing triangle integration (polynomial)..." << std::endl;
  
  // Arbitrary triangle
  std::array<std::array<Real, 2>, 3> vertices = {{
    {{0, 0}},
    {{2, 0}},
    {{1, 1}}
  }};
  
  // Polynomial f(x,y) = x^2 + y^2
  auto f = [](const Real* x) -> Real {
    return x[0]*x[0] + x[1]*x[1];
  };
  
  auto result = integrate_triangle<Real>(f, vertices, 1e-8, 1e-8);
  
  // For this specific triangle with vertices (0,0), (2,0), (1,1)
  // Area = |det([[2-0, 1-0], [0-0, 1-0]])|/2 = |det([[2,1],[0,1]])|/2 = |2|/2 = 1
  // We can compute this analytically using parametric integration
  // The exact value is 4/3
  Real exact = Real(4) / Real(3);  // Verified analytically
  Real error = std::abs(result.value - exact);
  
  std::cout << "  Computed: " << result.value << std::endl;
  std::cout << "  Expected: ~" << exact << std::endl;
  std::cout << "  Error: " << error << std::endl;
  assert(error < 1e-6);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_tetrahedron_constant() {
  std::cout << "Testing tetrahedron integration (constant)..." << std::endl;
  
  // Standard tetrahedron with vertices at origin and unit vectors
  std::array<std::array<Real, 3>, 4> vertices = {{
    {{0, 0, 0}},
    {{1, 0, 0}},
    {{0, 1, 0}},
    {{0, 0, 1}}
  }};
  
  // Constant function f(x,y,z) = 1
  auto f = [](const Real* x) -> Real {
    return Real(1);
  };
  
  auto result = integrate_tetrahedron<Real>(f, vertices, 1e-8, 1e-8);
  
  // Volume of standard tetrahedron = 1/6
  Real exact = Real(1) / Real(6);
  Real error = std::abs(result.value - exact);
  
  std::cout << "  Computed: " << result.value << std::endl;
  std::cout << "  Exact: " << exact << std::endl;
  std::cout << "  Error: " << error << std::endl;
  assert(error < 1e-8);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_tetrahedron_linear() {
  std::cout << "Testing tetrahedron integration (linear)..." << std::endl;
  
  // Standard tetrahedron
  std::array<std::array<Real, 3>, 4> vertices = {{
    {{0, 0, 0}},
    {{1, 0, 0}},
    {{0, 1, 0}},
    {{0, 0, 1}}
  }};
  
  // Linear function f(x,y,z) = x + y + z
  auto f = [](const Real* x) -> Real {
    return x[0] + x[1] + x[2];
  };
  
  auto result = integrate_tetrahedron<Real>(f, vertices, 1e-8, 1e-8);
  
  // Integral of (x + y + z) over standard tetrahedron
  // The integral should be 3*(1/24) where each coordinate contributes 1/24
  // Actually checking: ∫∫∫(x+y+z) over tetrahedron = volume * avg(x+y+z) = 1/6 * (3/4) = 1/8
  Real exact = Real(1) / Real(8);
  Real error = std::abs(result.value - exact);
  
  std::cout << "  Computed: " << result.value << std::endl;
  std::cout << "  Exact: " << exact << std::endl;
  std::cout << "  Error: " << error << std::endl;
  assert(error < 1e-8);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_general_simplex() {
  std::cout << "Testing general simplex integration..." << std::endl;
  
  // 4D simplex (5 vertices in 4D)
  simplex<Real> simp;
  simp.vertices = {
    {0, 0, 0, 0},
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}
  };
  
  // Constant function
  auto f = [](const Real* x) -> Real {
    return Real(1);
  };
  
  auto result = integrate_simplex<Real>(f, simp, 1e-8, 1e-8);
  
  // Volume of 4D standard simplex = 1/4! = 1/24
  Real exact = Real(1) / Real(24);
  Real error = std::abs(result.value - exact);
  
  std::cout << "  Computed: " << result.value << std::endl;
  std::cout << "  Exact: " << exact << std::endl;
  std::cout << "  Error: " << error << std::endl;
  assert(error < 1e-8);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_simplex_with_different_signatures() {
  std::cout << "Testing simplex with different integrand signatures..." << std::endl;
  
  // Standard triangle
  simplex<Real> tri;
  tri.vertices = {
    {0, 0},
    {1, 0},
    {0, 1}
  };
  
  // Test with vector signature
  auto f_vec = [](const std::vector<Real>& x) -> Real {
    return Real(1);
  };
  auto r1 = integrate_simplex<Real>(f_vec, tri, 1e-8, 1e-8);
  assert(std::abs(r1.value - Real(0.5)) < 1e-8);
  
  // Test with pointer + size signature
  auto f_ptr = [](const Real* x, std::size_t /*n*/) -> Real {
    return Real(1);
  };
  auto r2 = integrate_simplex<Real>(f_ptr, tri, 1e-8, 1e-8);
  assert(std::abs(r2.value - Real(0.5)) < 1e-8);
  
  // Test with raw pointer signature
  auto f_raw = [](const Real* x) -> Real {
    return Real(1);
  };
  auto r3 = integrate_simplex<Real>(f_raw, tri, 1e-8, 1e-8);
  assert(std::abs(r3.value - Real(0.5)) < 1e-8);
  
  std::cout << "  All signatures work correctly" << std::endl;
  std::cout << "  PASSED" << std::endl;
}

int main() {
  std::cout << "===== Simplex Integration Tests =====" << std::endl;
  
  using Real = double;
  
  test_triangle_constant<Real>();
  test_triangle_linear<Real>();
  test_triangle_polynomial<Real>();
  test_tetrahedron_constant<Real>();
  test_tetrahedron_linear<Real>();
  test_general_simplex<Real>();
  test_simplex_with_different_signatures<Real>();
  
  std::cout << "\n===== All Tests PASSED =====" << std::endl;
  
  return 0;
}