#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <boost/math/cubature/detail/gm_rules.hpp>

using namespace boost::math::cubature::detail::gm;

static bool approx_eq(long double a, long double b, long double tol=1e-12L) {
  return std::abs(a-b) <= tol * std::max<long double>(1.0L, std::abs(b));
}

template <int Deg, int Dim>
void check_sum_weights() {
  if constexpr (!rule_fam<Deg, Dim, family_9_7>::available) return;
  auto ws = rule_fam<Deg, Dim, family_9_7>::template weights<double>();
  long double s = 0.0L;
  for (auto w : ws) s += static_cast<long double>(w);
  if (!approx_eq(s, 1.0L)) {
    std::cerr << "[gm97] Dim=" << Dim << " Deg=" << Deg << " weight sum = " << s << " (expected 1)\n";
    std::abort();
  }
}

int main() {
  // Degree-7 across d=2..6
  check_sum_weights<7,2>();
  check_sum_weights<7,3>();
  check_sum_weights<7,4>();
  check_sum_weights<7,5>();
  check_sum_weights<7,6>();
  // Degree-9 across d=2..6
  check_sum_weights<9,2>();
  check_sum_weights<9,3>();
  check_sum_weights<9,4>();
  check_sum_weights<9,5>();
  check_sum_weights<9,6>();
  return 0;
}

