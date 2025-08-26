// Compile-only checks for gm_rules API
#include <type_traits>
#include <boost/math/cubature/detail/gm_rules.hpp>

using namespace boost::math::cubature::detail::gm;

static_assert(!available_v<7, 3>, "No concrete tables yet");
static_assert(!pair_available_v<7,5, 3>, "No concrete tables yet");

int main() {
  // Ensure types instantiate
  using P = pair_7_5<3>;
  (void)P::available;
  return 0;
}

