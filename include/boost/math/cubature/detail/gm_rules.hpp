#ifndef BOOST_MATH_CUBATURE_DETAIL_GM_RULES_HPP
#define BOOST_MATH_CUBATURE_DETAIL_GM_RULES_HPP

#include <array>
#include <cstddef>
#include <type_traits>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \file gm_rules.hpp
/// \brief Fully symmetric Genz–Malik embedded rules (API + metadata).
/// \details Provides compile-time metadata and accessors for degree-7/5 and 9/7
///  embedded cubature rules on the unit hypercube [0,1]^d (affine scaled by caller).

namespace gm {

// Tags for rule degrees
struct deg5  { static constexpr int value = 5; };
struct deg7  { static constexpr int value = 7; };
struct deg9  { static constexpr int value = 9; };

// Rule descriptor: primary template (not available by default)
template <int Degree, int Dim>
struct rule {
  static constexpr bool available = false;
  static constexpr int degree     = Degree;
  static constexpr int dim        = Dim;
  static constexpr std::size_t size = 0; // number of nodes

  template <class Real>
  static constexpr std::array<std::array<Real, Dim>, size> nodes() { return {}; }

  template <class Real>
  static constexpr std::array<Real, size> weights() { return {}; }
};

// Embedded pair helper: exposes fine/coarse rules sharing nodes
template <int FineDeg, int CoarseDeg, int Dim>
struct embedded_pair {
  using fine   = rule<FineDeg, Dim>;
  using coarse = rule<CoarseDeg, Dim>;
  static constexpr bool available = fine::available && coarse::available;
};

// Convenience alias for common pairs
template <int Dim>
using pair_7_5 = embedded_pair<7,5,Dim>;

template <int Dim>
using pair_9_7 = embedded_pair<9,7,Dim>;

// Metadata helpers
template <int Degree, int Dim>
inline constexpr bool available_v = rule<Degree, Dim>::available;

template <int Fine, int Coarse, int Dim>
inline constexpr bool pair_available_v = embedded_pair<Fine, Coarse, Dim>::available;

// NOTE: Concrete specializations with constexpr nodes/weights for (7/5) and (9/7)
// up to Dim=15 will be added in a follow-up patch. The API above is stable.

} // namespace gm

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_GM_RULES_HPP

