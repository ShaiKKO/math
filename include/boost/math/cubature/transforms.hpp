#ifndef BOOST_MATH_CUBATURE_TRANSFORMS_HPP
#define BOOST_MATH_CUBATURE_TRANSFORMS_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>

namespace boost { namespace math { namespace cubature {

/// \file transforms.hpp
/// \brief Domain and variance-reduction transforms (API skeleton).
/// \details Provides finite↔infinite mappings, Duffy mapping helpers, and the
///  tent (Baker) transform used by QMC.

// Tent/Baker transform helper
template <class Real>
inline Real tent(Real u) {
  return Real(1) - Real(2) * (u > Real(0.5) ? (u - Real(0.5)) : (Real(0.5) - u));
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_TRANSFORMS_HPP

