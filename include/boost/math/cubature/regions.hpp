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
/// \brief Domain types for cubature integration
/// \details Provides domain representations for hypercube and simplex regions.
///          These are the primary input types for integration functions.
/// 
/// \section regions_usage Usage
/// \code
/// // Create a 2D unit square [0,1] x [0,1]
/// hypercube<double> unit_square(2);
/// std::fill(unit_square.lower.begin(), unit_square.lower.end(), 0.0);
/// std::fill(unit_square.upper.begin(), unit_square.upper.end(), 1.0);
/// 
/// // Create a 2D triangle with vertices at (0,0), (1,0), (0,1)
/// simplex<double> triangle;
/// triangle.vertices = {{0, 0}, {1, 0}, {0, 1}};
/// 
/// // Integrate over regions
/// auto f = [](const double* x) { return x[0] * x[1]; };
/// auto result1 = integrate_adaptive(f, unit_square, 1e-6, 1e-6);
/// auto result2 = integrate_simplex(f, triangle, 1e-6, 1e-6);
/// \endcode

/// \brief Hypercube domain [lower, upper] in R^d
/// \tparam Real Floating point type
/// \details Represents an axis-aligned box in d-dimensional space.
///          Each dimension i has bounds [lower[i], upper[i]].
template <class Real>
struct hypercube {
  std::vector<Real> lower; ///< Lower bounds for each dimension
  std::vector<Real> upper; ///< Upper bounds for each dimension
  
  /// \brief Default constructor
  hypercube() = default;
  
  /// \brief Constructor with dimension
  /// \param dim Number of dimensions
  explicit hypercube(std::size_t dim) : lower(dim), upper(dim) {}
  
  /// \brief Get dimension of the hypercube
  /// \return Number of dimensions
  std::size_t dimension() const { return lower.size(); }
  
  /// \brief Compute volume of the hypercube
  /// \return Product of side lengths: ∏(upper[i] - lower[i])
  Real volume() const {
    Real vol = Real(1);
    for (std::size_t i = 0; i < lower.size(); ++i) {
      vol *= (upper[i] - lower[i]);
    }
    return vol;
  }
};

/// \brief General simplex: d+1 vertices in R^d
/// \tparam Real Floating point type
/// \details A d-dimensional simplex is the convex hull of d+1 vertices.
///          Examples: 2D triangle (3 vertices), 3D tetrahedron (4 vertices).
///          Each vertex is a d-dimensional point.
/// \note The vertices should be affinely independent (non-degenerate simplex)
template <class Real>
struct simplex {
  std::vector<std::vector<Real>> vertices;  ///< List of d+1 vertices, each of dimension d
  
  /// \brief Default constructor
  simplex() = default;
  
  /// \brief Get dimension of the simplex
  /// \return d where simplex is d-dimensional (number of vertices - 1)
  std::size_t dimension() const { 
    return vertices.empty() ? 0 : vertices.size() - 1; 
  }
  
  /// \brief Check if simplex is valid
  /// \return true if all vertices have consistent dimension
  bool is_valid() const {
    if (vertices.empty()) return true;
    std::size_t dim = vertices[0].size();
    for (const auto& v : vertices) {
      if (v.size() != dim) return false;
    }
    return vertices.size() == dim + 1;
  }
};

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_REGIONS_HPP

