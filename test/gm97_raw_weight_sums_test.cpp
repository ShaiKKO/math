#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <boost/math/cubature/detail/gm_rules.hpp>

using namespace boost::math::cubature::detail::gm;

static bool approx_eq(long double a, long double b, long double tol=1e-12L) {
  return std::abs(a-b) <= tol * std::max<long double>(1.0L, std::abs(b));
}

int main(){
  // Only checking 2D and 3D raw (have raw specializations)
  if constexpr (raw_rule_fam<9,2,family_9_7>::available) {
    auto ws = raw_rule_fam<9,2,family_9_7>::template weights_with_zero<double>();
    long double s = 0.0L; for (auto w: ws) s += w; if (!approx_eq(s, 1.0L)) std::abort();
  }
  if constexpr (raw_rule_fam<9,3,family_9_7>::available) {
    auto ws = raw_rule_fam<9,3,family_9_7>::template weights_with_zero<double>();
    long double s = 0.0L; for (auto w: ws) s += w; if (!approx_eq(s, 1.0L)) std::abort();
  }
  return 0;
}

