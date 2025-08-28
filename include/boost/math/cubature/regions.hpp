// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


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
  
  // Default constructor
  hypercube() = default;
  
  // Constructor with dimension
  explicit hypercube(std::size_t dim) : lower(dim), upper(dim) {}
  
  // Get dimension
  std::size_t dimension() const { return lower.size(); }
  
  // Get volume
  Real volume() const {
    Real vol = Real(1);
    for (std::size_t i = 0; i < lower.size(); ++i) {
      vol *= (upper[i] - lower[i]);
    }
    return vol;
  }
};

/// General simplex: d+1 vertices in R^d (each vertex is length d)
template <class Real>
struct simplex {
  std::vector<std::vector<Real>> vertices;
};

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_REGIONS_HPP

