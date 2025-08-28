#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <boost/math/cubature/detail/gm_rules.hpp>

using namespace boost::math::cubature::detail::gm;

static bool approx_eq(long double a, long double b, long double tol=1e-12L) {
  return std::abs(a-b) <= tol * std::max<long double>(1.0L, std::abs(b));
}

template <class Real, int Deg, int Dim>
void check_rule_exactness_fam()
{
  if constexpr (!rule_fam<Deg, Dim, family_9_7>::available) return;

  auto nodes = rule_fam<Deg, Dim, family_9_7>::template nodes<Real>();
  auto weights = rule_fam<Deg, Dim, family_9_7>::template weights<Real>();

  long double sumw = 0.0L;
  for (auto w : weights) sumw += static_cast<long double>(w);
  if (!approx_eq(sumw, 1.0L)) {
    std::cerr << "[gm97] weight sum mismatch: Dim=" << Dim << " Deg=" << Deg
              << " sumw=" << sumw << "\n";
    std::abort();
  }

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

  std::array<unsigned, Dim> a{};
  struct Rec {
    decltype(quad_monomial)& quad; decltype(ref_monomial)& ref;
    std::array<unsigned,Dim>& a; int DimC; int DegC;
    void run(int idx, unsigned sum){
      if (idx == DimC) {
        if (sum <= static_cast<unsigned>(DegC)) {
          long double r = ref(a);
          long double q = quad(a);
          if (!approx_eq(r, q)) {
            std::cerr << "[gm97] exactness fail Dim=" << DimC << " Deg=" << DegC << " a=";
            for (int i=0;i<DimC;++i) std::cerr << (i?",":"[") << a[i];
            std::cerr << "] r=" << r << " q=" << q << " diff=" << std::abs(r-q) << "\n";
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
  // Degree-7 exactness across dims 2..6
  check_rule_exactness_fam<Real, 7, 2>();
  check_rule_exactness_fam<Real, 7, 3>();
  check_rule_exactness_fam<Real, 7, 4>();
  check_rule_exactness_fam<Real, 7, 5>();
  check_rule_exactness_fam<Real, 7, 6>();
  // Degree-9 exactness across dims 2..6
  check_rule_exactness_fam<Real, 9, 2>();
  check_rule_exactness_fam<Real, 9, 3>();
  check_rule_exactness_fam<Real, 9, 4>();
  check_rule_exactness_fam<Real, 9, 5>();
  check_rule_exactness_fam<Real, 9, 6>();
  return 0;
}

