// Compile-only tests for cubature API shapes
#include <vector>
#include <cstddef>

#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/qmc.hpp>
#include <boost/math/cubature/simplex.hpp>
#include <boost/math/cubature/workspace.hpp>

using namespace boost::math::cubature;

namespace {

// Several callable shapes

template <class Real>
struct scalar_vec {
  Real operator()(const std::vector<Real>& x) const { return x.empty()? Real(0): x[0]; }
};

template <class Real>
struct vector_vec {
  void operator()(const std::vector<Real>& x, Real* out, std::size_t m) const {
    if (m) out[0] = x.empty()? Real(0): x[0];
    if (m>1) out[1] = Real(1);
  }
};

template <class Real>
struct scalar_ptr {
  Real operator()(const Real* p, std::size_t n) const { return n? p[0]: Real(0); }
};

template <class Real>
struct vector_ptr {
  void operator()(const Real* p, std::size_t n, Real* out, std::size_t m) const {
    (void)n; if (m) out[0] = p? p[0]: Real(0);
  }
};

} // anonymous

int main() {
  using Real = double;
  hypercube<Real> box{{0.0, 0.0}, {1.0, 1.0}};
  simplex<Real> tri{ { {0.0,0.0}, {1.0,0.0}, {0.0,1.0} } };

  execution_options ex{}; workspace ws{};

  // Adaptive
  (void)integrate_adaptive<Real>(scalar_vec<Real>{}, box, 1e-8, 1e-8);
  (void)integrate_adaptive<Real>(vector_vec<Real>{}, box, 1e-8, 1e-8);
  (void)integrate_adaptive<Real>(scalar_ptr<Real>{}, box, 1e-8, 1e-8, 0);

  (void)integrate_adaptive<Real>(scalar_vec<Real>{}, box, 1e-8, 1e-8, 0);
  (void)integrate_adaptive<Real>(scalar_vec<Real>{}, box, 1e-8, 1e-8, ex, 0);

  // Sparse-grid
  (void)integrate_sparse_grid<Real>(scalar_vec<Real>{}, box, 3);
  (void)integrate_sparse_grid<Real>(scalar_vec<Real>{}, box, 3, ex, &ws);

  // QMC
  (void)integrate_qmc<Real>(vector_vec<Real>{}, box, 1<<10, 2, true);
  (void)integrate_qmc<Real>(vector_vec<Real>{}, box, 1<<10, 4, true, ex, &ws);

  // Simplex
  (void)integrate_simplex<Real>(scalar_vec<Real>{}, tri, 1e-8, 1e-8);
  (void)integrate_simplex<Real>(scalar_vec<Real>{}, tri, 1e-8, 1e-8, ex);

  return 0;
}

