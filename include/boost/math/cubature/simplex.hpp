// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_SIMPLEX_HPP
#define BOOST_MATH_CUBATURE_SIMPLEX_HPP

// STL headers first per Boost conventions
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// Project headers
#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/concepts.hpp>
#include <boost/math/cubature/workspace.hpp>
#include <boost/math/cubature/transforms.hpp>
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>

namespace boost { namespace math { namespace cubature {

/// \file simplex.hpp
/// \brief Simplex-domain integrator using Duffy transform
/// \details Maps simplex to unit hypercube via Duffy transform,
///          then calls adaptive or sparse grid backend

namespace detail {

/// \brief Compute affine transform from reference simplex to general simplex
template <typename Real>
class simplex_affine_map {
private:
  std::vector<std::vector<Real>> matrix_;  // Transform matrix
  std::vector<Real> offset_;              // Translation vector
  Real determinant_;                      // Jacobian determinant
  std::size_t dim_;
  
public:
  simplex_affine_map(const simplex<Real>& s) {
    // Get dimension (number of vertices - 1)
    dim_ = s.vertices.size() - 1;
    if (dim_ == 0 || s.vertices.empty()) {
      determinant_ = Real(1);
      return;
    }
    
    // Build transformation matrix from reference simplex to given simplex
    // Reference simplex has vertices at origin and unit vectors
    matrix_.resize(dim_, std::vector<Real>(dim_));
    offset_ = s.vertices[0];  // First vertex is the offset
    
    // Matrix columns are edge vectors from first vertex
    for (std::size_t i = 0; i < dim_; ++i) {
      for (std::size_t j = 0; j < dim_; ++j) {
        matrix_[j][i] = s.vertices[i + 1][j] - s.vertices[0][j];
      }
    }
    
    // Compute determinant (volume scaling factor)
    determinant_ = compute_determinant(matrix_);
  }
  
  /// \brief Map point from reference simplex to general simplex
  void map_to_simplex(const Real* ref_point, Real* point) const {
    // point = offset + matrix * ref_point
    for (std::size_t i = 0; i < dim_; ++i) {
      point[i] = offset_[i];
      for (std::size_t j = 0; j < dim_; ++j) {
        point[i] += matrix_[i][j] * ref_point[j];
      }
    }
  }
  
  Real jacobian() const { return std::abs(determinant_); }
  
private:
  /// \brief Compute determinant of matrix (LU decomposition)
  Real compute_determinant(const std::vector<std::vector<Real>>& A) const {
    if (A.empty()) return Real(1);
    if (A.size() == 1) return A[0][0];
    if (A.size() == 2) {
      return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    }
    
    // LU decomposition for larger matrices
    std::vector<std::vector<Real>> LU = A;
    Real det = Real(1);
    
    for (std::size_t k = 0; k < A.size(); ++k) {
      // Find pivot
      std::size_t pivot = k;
      for (std::size_t i = k + 1; i < A.size(); ++i) {
        if (std::abs(LU[i][k]) > std::abs(LU[pivot][k])) {
          pivot = i;
        }
      }
      
      // Swap rows if needed
      if (pivot != k) {
        std::swap(LU[k], LU[pivot]);
        det = -det;  // Swap changes sign
      }
      
      // Check for singularity
      if (std::abs(LU[k][k]) < std::numeric_limits<Real>::epsilon()) {
        return Real(0);
      }
      
      det *= LU[k][k];
      
      // Eliminate column
      for (std::size_t i = k + 1; i < A.size(); ++i) {
        LU[i][k] /= LU[k][k];
        for (std::size_t j = k + 1; j < A.size(); ++j) {
          LU[i][j] -= LU[i][k] * LU[k][j];
        }
      }
    }
    
    return det;
  }
};

/// \brief Simplex integrand wrapper with Duffy transform
template <typename Real, typename F>
class duffy_simplex_integrand {
private:
  const F& f_;
  const simplex_affine_map<Real>& affine_;
  std::size_t dim_;
  mutable std::vector<Real> ref_point_;
  mutable std::vector<Real> simplex_point_;
  
public:
  duffy_simplex_integrand(const F& f, const simplex_affine_map<Real>& map, 
                          std::size_t dim)
    : f_(f), affine_(map), dim_(dim), 
      ref_point_(dim), simplex_point_(dim) {}
  
  Real operator()(const Real* u) const {
    // Apply Duffy transform: [0,1]^d → reference simplex
    Real duffy_jac = duffy_transform<Real>::apply(u, ref_point_.data(), dim_);
    
    // Apply affine map: reference simplex → general simplex
    affine_.map_to_simplex(ref_point_.data(), simplex_point_.data());
    Real affine_jac = affine_.jacobian();
    
    // Evaluate integrand with combined Jacobian
    Real value = evaluate_integrand(simplex_point_.data());
    return value * duffy_jac * affine_jac;
  }
  
  Real operator()(const Real* u, std::size_t n) const {
    return (*this)(u);
  }
  
  Real operator()(const std::vector<Real>& u) const {
    return (*this)(u.data());
  }
  
private:
  Real evaluate_integrand(const Real* x) const {
    // Support multiple integrand signatures
    if constexpr (std::is_invocable_v<F, const Real*, std::size_t>) {
      return f_(x, dim_);
    } else if constexpr (std::is_invocable_v<F, const Real*>) {
      return f_(x);
    } else {
      std::vector<Real> vec(x, x + dim_);
      return f_(vec);
    }
  }
};

} // namespace detail

/// \brief Integrate over a simplex domain
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_simplex(const F& f, const simplex<Real>& simp,
                                      Real abs_tol, Real rel_tol,
                                      Policy const& pol = Policy{})
{
  // Validate simplex
  if (simp.vertices.empty()) {
    result<Real> res;
    res.status = status_code::invalid_input;
    res.error = std::numeric_limits<Real>::max();
    return res;
  }
  
  std::size_t dim = simp.vertices.size() - 1;
  if (dim == 0) {
    // 0-dimensional simplex (point)
    result<Real> res;
    res.value = Real(1);
    res.error = Real(0);
    res.evaluations = 1;
    res.status = status_code::success;
    return res;
  }
  
  // Check that all vertices have correct dimension
  for (const auto& vertex : simp.vertices) {
    if (vertex.size() != dim) {
      result<Real> res;
      res.status = status_code::invalid_input;
      res.error = std::numeric_limits<Real>::max();
      return res;
    }
  }
  
  // Create affine map from reference simplex to given simplex
  detail::simplex_affine_map<Real> affine_map(simp);
  
  // Create transformed integrand with Duffy transform
  detail::duffy_simplex_integrand<Real, F> transformed_f(f, affine_map, dim);
  
  // Create unit hypercube for integration
  hypercube<Real> unit_cube(dim);
  std::fill(unit_cube.lower.begin(), unit_cube.lower.end(), Real(0));
  std::fill(unit_cube.upper.begin(), unit_cube.upper.end(), Real(1));
  
  // Call adaptive integration backend
  // Could also use sparse grid for smooth integrands
  return integrate_adaptive<Real>(transformed_f, unit_cube, abs_tol, rel_tol, 
                                  100000, pol);
}

/// \brief Integrate over a simplex with execution options
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_simplex(const F& f, const simplex<Real>& simp,
                                      Real abs_tol, Real rel_tol,
                                      execution_options const& /*ex*/,
                                      Policy const& pol = Policy{})
{
  // For now, ignore execution options and forward to main implementation
  return integrate_simplex<Real, F, Policy>(f, simp, abs_tol, rel_tol, pol);
}

/// \brief Convenience function for 2D triangle integration
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_triangle(const F& f,
                                       const std::array<std::array<Real, 2>, 3>& vertices,
                                       Real abs_tol, Real rel_tol,
                                       Policy const& pol = Policy{})
{
  simplex<Real> tri;
  tri.vertices.resize(3);
  for (std::size_t i = 0; i < 3; ++i) {
    tri.vertices[i] = std::vector<Real>(vertices[i].begin(), vertices[i].end());
  }
  return integrate_simplex<Real, F, Policy>(f, tri, abs_tol, rel_tol, pol);
}

/// \brief Convenience function for 3D tetrahedron integration
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_tetrahedron(const F& f,
                                          const std::array<std::array<Real, 3>, 4>& vertices,
                                          Real abs_tol, Real rel_tol,
                                          Policy const& pol = Policy{})
{
  simplex<Real> tet;
  tet.vertices.resize(4);
  for (std::size_t i = 0; i < 4; ++i) {
    tet.vertices[i] = std::vector<Real>(vertices[i].begin(), vertices[i].end());
  }
  return integrate_simplex<Real, F, Policy>(f, tet, abs_tol, rel_tol, pol);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_SIMPLEX_HPP

