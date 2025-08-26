// Example: scalar integrand with adaptive hypercube cubature (compile-only)
#include <cmath>
#include <vector>
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/regions.hpp>

using namespace boost::math::cubature;

int main() {
  using Real = double;
  hypercube<Real> box{{-5.0, -5.0}, {5.0, 5.0}};

  auto f = [](const std::vector<Real>& x)->Real { return std::exp(-(x[0]*x[0] + x[1]*x[1])); };
  auto res = integrate_adaptive<Real>(f, box, 1e-10, 1e-10);
  (void)res;
  return 0;
}

