#include <cassert>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <iostream>

#include <boost/math/cubature/detail/gm_rules.hpp>

using namespace boost::math::cubature::detail::gm;

template <class Real, int Deg, int Dim>
void check_rule_exactness()
{
#if __cplusplus >= 201703L
  if constexpr (!available_v<Deg, Dim>) {
    return; // skip silently
  }
#endif

  auto nodes = rule<Deg, Dim>::template nodes<Real>();
  auto weights = rule<Deg, Dim>::template weights<Real>();

  long double sumw = 0.0L;
  for (auto w : weights) sumw += static_cast<long double>(w);
  if (!(std::abs(sumw - 1.0L) < 1e-12L)) {
    std::cerr << "Weight sum mismatch for Deg=" << Deg << ", Dim=" << Dim
              << ": sumw=" << sumw << "\n";
    std::abort();
  }

  auto ref_monomial = [](const std::array<unsigned, Dim>& a) {
    long double ref = 1.0L;
    for (int i = 0; i < Dim; ++i) ref *= 1.0L / static_cast<long double>(a[i] + 1);
    return ref;
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

  std::array<unsigned, Dim> a{};
  struct Rec {
    decltype(quad_monomial)& quad; decltype(ref_monomial)& ref;
    std::array<unsigned,Dim>& a; int DimC; int DegC;
    void run(int idx, unsigned sum){
      if (idx == DimC) {
        if (sum <= static_cast<unsigned>(DegC)) {
          long double r = ref(a);
          long double q = quad(a);
          long double diff = std::abs(r - q);
          long double tol = 1e-12L * std::max<long double>(1.0L, std::abs(r));
          if (!(diff <= tol)) {
            std::cerr << "Exactness fail Deg=" << DegC << ", Dim=" << DimC << ", a=";
            for (int i=0;i<DimC;++i) std::cerr << (i?",":"[") << a[i];
            std::cerr << "] r=" << r << " q=" << q << " diff=" << diff << " tol=" << tol << "\n";
            std::abort();
          }
        }
        return;
      }
      for (unsigned k = 0; k <= static_cast<unsigned>(DegC); ++k) { a[idx]=k; run(idx+1, sum+k);}    }
  };
  Rec{quad_monomial, ref_monomial, a, Dim, Deg}.run(0,0);
}

int main()
{
  using Real = double;
  // Deg-5 (embedded) checks for 2D..5D
  check_rule_exactness<Real, 5, 2>();
  check_rule_exactness<Real, 5, 3>();
  check_rule_exactness<Real, 5, 4>();
  check_rule_exactness<Real, 5, 5>();
  // Deg-7 checks for 2D..5D
  check_rule_exactness<Real, 7, 2>();
  check_rule_exactness<Real, 7, 3>();
  check_rule_exactness<Real, 7, 4>();
  check_rule_exactness<Real, 7, 5>();
  return 0;
}

