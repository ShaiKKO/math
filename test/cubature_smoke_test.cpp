#include <cassert>
#include <vector>

#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/qmc.hpp>
#include <boost/math/cubature/simplex.hpp>
#include <boost/math/cubature/transforms.hpp>

using namespace boost::math::cubature;

int main() {
  using Real = double;
  hypercube<Real> box{{0.0, 0.0}, {1.0, 1.0}};

  // For C++17 portability in this smoke test, avoid std::span
  auto f_scalar = [](const std::vector<Real>& x) -> Real { return x[0] + x[1]; };
  auto r1 = integrate_adaptive<Real>(f_scalar, box, 1e-8, 1e-8);
  (void)r1;

  auto r2 = integrate_sparse_grid<Real>(f_scalar, box, 3);
  (void)r2;

  auto f_vec = [](const std::vector<Real>& x, Real* out, std::size_t m) {
    out[0] = x[0] + x[1];
    if (m > 1) out[1] = x[0] * x[1];
  };
  auto r3 = integrate_qmc<Real>(f_vec, box, 1 << 10, 2, true);
  (void)r3;

  simplex<Real> tri{ { {0.0,0.0}, {1.0,0.0}, {0.0,1.0} } };
  auto r4 = integrate_simplex<Real>(f_scalar, tri, 1e-8, 1e-8);
  (void)r4;

  return 0;
}

