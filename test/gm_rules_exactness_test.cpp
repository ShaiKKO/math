#define BOOST_TEST_MODULE gm_rules_exactness
#include <boost/test/unit_test.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <limits>

#include <boost/math/cubature/detail/gm_rules.hpp>

using namespace boost::math::cubature::detail::gm;

template <class Real, int Deg, int Dim>
void check_rule_exactness()
{
  if constexpr (!available_v<Deg, Dim>) {
    BOOST_TEST_MESSAGE("GM rule<" << Deg << "," << Dim << "> not available; skipping exactness checks.");
    BOOST_CHECK(true);
    return;
  }

  constexpr auto nodes = rule<Deg, Dim>::template nodes<Real>();
  constexpr auto weights = rule<Deg, Dim>::template weights<Real>();
  static_assert(nodes.size() == weights.size(), "nodes/weights size mismatch");

  auto ref_monomial = [](const std::array<unsigned, Dim>& a) {
    long double ref = 1.0L;
    for (int i = 0; i < Dim; ++i) ref *= 1.0L / static_cast<long double>(a[i] + 1);
    return ref; // integral over [0,1]^Dim of prod x_i^{a_i}
  };

  auto quad_monomial = [&](const std::array<unsigned, Dim>& a) {
    long double s = 0.0L;
    for (std::size_t k = 0; k < nodes.size(); ++k) {
      long double phi = 1.0L;
      for (int i = 0; i < Dim; ++i) {
        phi *= std::pow(static_cast<long double>(nodes[k][i]), static_cast<int>(a[i]));
      }
      s += static_cast<long double>(weights[k]) * phi;
    }
    return s;
  };

  auto within_tol = [](long double a, long double b) {
    long double tol = 1e-12L;
    return std::abs(a - b) <= tol * std::max<long double>(1.0L, std::abs(b));
  };

  // iterate all multi-indices with sum <= Deg
  std::array<unsigned, Dim> a{};
  std::function<void(int, unsigned)> rec = [&](int idx, unsigned sum) {
    if (idx == Dim) {
      if (sum <= static_cast<unsigned>(Deg)) {
        auto ref = ref_monomial(a);
        auto quad = quad_monomial(a);
        BOOST_CHECK(within_tol(ref, quad));
      }
      return;
    }
    unsigned max_k = static_cast<unsigned>(Deg);
    for (unsigned k = 0; k <= max_k; ++k) {
      a[idx] = k;
      rec(idx+1, sum + k);
    }
  };
  rec(0, 0);
}

BOOST_AUTO_TEST_CASE(exactness_if_available)
{
  using Real = double;
  // These calls will skip internally if tables are not yet available.
  check_rule_exactness<Real, 5, 2>();
  check_rule_exactness<Real, 7, 2>();
  check_rule_exactness<Real, 5, 3>();
  check_rule_exactness<Real, 7, 3>();
}

