#ifndef BOOST_MATH_CUBATURE_REGIONS_HPP
#define BOOST_MATH_CUBATURE_REGIONS_HPP

#include <vector>

namespace boost { namespace math { namespace cubature {

/// \file regions.hpp
/// \brief Domain types for cubature front-ends.

/// Hypercube domain [lower, upper] in R^d
template <class Real>
struct hypercube {
  std::vector<Real> lower; ///< lower[i]
  std::vector<Real> upper; ///< upper[i]
};

/// General simplex: d+1 vertices in R^d (each vertex is length d)
template <class Real>
struct simplex {
  std::vector<std::vector<Real>> vertices;
};

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_REGIONS_HPP

