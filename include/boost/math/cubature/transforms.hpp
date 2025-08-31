// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_TRANSFORMS_HPP
#define BOOST_MATH_CUBATURE_TRANSFORMS_HPP

// STL headers first per Boost conventions
#include <cmath>
#include <vector>
#include <array>
#include <limits>
#include <algorithm>
#include <numeric>

// Boost headers
#include <boost/math/constants/constants.hpp>
#include <boost/math/cubature/constants.hpp>
#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>

namespace boost { namespace math { namespace cubature {

/// \file transforms.hpp
/// \brief Domain transforms for infinite bounds and simplex integration
/// \details Provides bijective mappings for infinite/semi-infinite domains,
///          Duffy transform for simplex domains, and the tent transform for QMC

/// \brief Transform types for infinite domain integration
enum class infinite_transform_type {
  rational,     // u/(1-u) for [0,∞), best conditioning
  tangent,      // tan(π(u-1/2)) for (-∞,∞), symmetric
  exponential,  // -log(1-u) for [0,∞), slow growth
  algebraic     // u^α/(1-u)^β for endpoint singularities
};

/// \brief Tent/Baker transform for variance reduction in QMC
template <class Real>
inline Real tent(Real u) {
  return Real(1) - Real(2) * std::abs(u - Real(0.5));
}

/// \brief Base class for domain transforms
template <typename Real>
class domain_transform {
public:
  virtual ~domain_transform() = default;
  
  /// \brief Apply forward transform from [0,1] to target domain
  /// \param u Unit interval value
  /// \return Transformed value and Jacobian
  virtual std::pair<Real, Real> forward(Real u) const = 0;
  
  /// \brief Apply inverse transform from target domain to [0,1]
  virtual Real inverse(Real x) const = 0;
};

/// \brief Rational transform: [0,1] → [0,∞)
/// \details x = u/(1-u), jacobian = 1/(1-u)²
template <typename Real>
class rational_transform : public domain_transform<Real> {
public:
  std::pair<Real, Real> forward(Real u) const override {
    if (u >= Real(1)) {
      return {std::numeric_limits<Real>::infinity(), 
              std::numeric_limits<Real>::infinity()};
    }
    Real one_minus_u = Real(1) - u;
    Real x = u / one_minus_u;
    Real jacobian = Real(1) / (one_minus_u * one_minus_u);
    return {x, jacobian};
  }
  
  Real inverse(Real x) const override {
    return x / (Real(1) + x);
  }
};

/// \brief Tangent transform: [0,1] → (-∞,∞)
/// \details x = tan(π(u-1/2)), jacobian = π/cos²(π(u-1/2))
template <typename Real>
class tangent_transform : public domain_transform<Real> {
public:
  std::pair<Real, Real> forward(Real u) const override {
    using std::tan;
    using std::cos;
    
    Real pi = boost::math::constants::pi<Real>();
    Real arg = pi * (u - Real(0.5));
    Real x = tan(arg);
    Real cos_arg = cos(arg);
    Real jacobian = pi / (cos_arg * cos_arg);
    return {x, jacobian};
  }
  
  Real inverse(Real x) const override {
    using std::atan;
    Real pi = boost::math::constants::pi<Real>();
    return atan(x) / pi + Real(0.5);
  }
};

/// \brief Exponential transform: [0,1] → [0,∞)
/// \details x = -log(1-u), jacobian = 1/(1-u)
template <typename Real>
class exponential_transform : public domain_transform<Real> {
public:
  std::pair<Real, Real> forward(Real u) const override {
    using std::log;
    
    if (u >= Real(1)) {
      return {std::numeric_limits<Real>::infinity(),
              std::numeric_limits<Real>::infinity()};
    }
    Real one_minus_u = Real(1) - u;
    Real x = -log(one_minus_u);
    Real jacobian = Real(1) / one_minus_u;
    return {x, jacobian};
  }
  
  Real inverse(Real x) const override {
    using std::exp;
    return Real(1) - exp(-x);
  }
};

/// \brief Algebraic transform with endpoint control
/// \details x = u^α/(1-u)^β, for controlling endpoint singularities
template <typename Real>
class algebraic_transform : public domain_transform<Real> {
private:
  Real alpha_;
  Real beta_;
  
public:
  algebraic_transform(Real alpha = Real(1), Real beta = Real(1))
    : alpha_(alpha), beta_(beta) {}
  
  std::pair<Real, Real> forward(Real u) const override {
    using std::pow;
    
    if (u <= Real(0)) {
      return {Real(0), Real(0)};
    }
    if (u >= Real(1)) {
      return {std::numeric_limits<Real>::infinity(),
              std::numeric_limits<Real>::infinity()};
    }
    
    Real one_minus_u = Real(1) - u;
    Real x = pow(u, alpha_) / pow(one_minus_u, beta_);
    
    // Jacobian = u^(α-1) * (1-u)^(-β-1) * [α(1-u) + βu]
    Real jacobian = pow(u, alpha_ - Real(1)) * 
                   pow(one_minus_u, -beta_ - Real(1)) *
                   (alpha_ * one_minus_u + beta_ * u);
    
    return {x, jacobian};
  }
  
  Real inverse(Real x) const override {
    // Solve x = u^α/(1-u)^β for u
    // For α=β=1, closed form: u = x/(1+x)
    if (alpha_ == Real(1) && beta_ == Real(1)) {
      return x / (Real(1) + x);
    }
    
    // General case: Newton-Raphson iteration
    // We need to solve f(u) = u^α/(1-u)^β - x = 0
    // f'(u) = α*u^(α-1)*(1-u)^(-β-1)*[(1-u) + β*u/(α)]
    //       = u^(α-1)*(1-u)^(-β-1)*[α*(1-u) + β*u]
    
    // Initial guess using linear approximation
    Real u = x / (Real(1) + x);  // Works well for moderate x
    if (x > Real(10)) {
      // For large x, u is close to 1
      u = Real(1) - Real(1) / std::pow(x, Real(1) / beta_);
    } else if (x < Real(0.1)) {
      // For small x, u is close to 0
      u = std::pow(x, Real(1) / alpha_);
    }
    
    // Ensure initial guess is in valid range
    const Real eps = std::numeric_limits<Real>::epsilon();
    u = std::max(eps, std::min(u, Real(1) - eps));
    
    // Newton iteration
    const int max_iter = 50;
    const Real tol = Real(boost::math::cubature::constants::machine_epsilon_safety_factor) * eps;
    
    for (int iter = 0; iter < max_iter; ++iter) {
      Real one_minus_u = Real(1) - u;
      
      // Compute f(u) = u^α/(1-u)^β - x
      Real f_val = std::pow(u, alpha_) / std::pow(one_minus_u, beta_) - x;
      
      // Compute f'(u)
      Real f_deriv = std::pow(u, alpha_ - Real(1)) * 
                     std::pow(one_minus_u, -beta_ - Real(1)) *
                     (alpha_ * one_minus_u + beta_ * u);
      
      // Check for convergence
      if (std::abs(f_val) < tol * (Real(1) + std::abs(x))) {
        break;
      }
      
      // Newton step with safeguarding
      Real delta = f_val / f_deriv;
      Real u_new = u - delta;
      
      // Safeguard to keep u in (0,1)
      if (u_new <= Real(0) || u_new >= Real(1)) {
        // Bisection step as fallback
        if (f_val > Real(0)) {
          // u is too large
          u = u * Real(0.5);
        } else {
          // u is too small
          u = Real(0.5) * (u + Real(1));
        }
      } else {
        u = u_new;
      }
      
      // Additional bounds check
      u = std::max(eps, std::min(u, Real(1) - eps));
    }
    
    return u;
  }
};

/// \brief Duffy transform for simplex integration
/// \details Maps unit hypercube [0,1]^d to standard simplex
template <typename Real>
class duffy_transform {
public:
  /// \brief Apply Duffy transform for d-dimensional simplex
  /// \param u Input point in [0,1]^d
  /// \param x Output point in simplex
  /// \param dim Dimension
  /// \return Jacobian of transformation
  static Real apply(const Real* u, Real* x, std::size_t dim) {
    if (dim == 0) return Real(1);
    
    // Standard Duffy transform for reference simplex
    // Maps [0,1]^d to d-simplex with vertices at origin and unit vectors
    // 2D: (u,v) -> (u(1-v), uv) with Jacobian = u
    // 3D: (u,v,w) -> (u(1-v)(1-w), uv(1-w), uvw) with Jacobian = u²v
    
    if (dim == 1) {
      x[0] = u[0];
      return Real(1);
    }
    
    if (dim == 2) {
      x[0] = u[0] * (Real(1) - u[1]);
      x[1] = u[0] * u[1];
      return u[0];
    }
    
    if (dim == 3) {
      // Standard mapping for 3-simplex
      // (u,v,w) -> (u(1-v), uv(1-w), uvw)
      x[0] = u[0] * (Real(1) - u[1]);
      x[1] = u[0] * u[1] * (Real(1) - u[2]);
      x[2] = u[0] * u[1] * u[2];
      return u[0] * u[0] * u[1];  // Jacobian: u²v
    }
    
    // General case
    std::vector<Real> prod(dim);
    prod[0] = u[0];
    for (std::size_t i = 1; i < dim; ++i) {
      prod[i] = prod[i-1] * u[i];
    }
    
    // Build coordinates
    for (std::size_t i = 0; i < dim; ++i) {
      x[i] = (i == 0 ? u[0] : prod[i-1] * u[i]);
      for (std::size_t j = i + 1; j < dim; ++j) {
        x[i] *= (Real(1) - u[j]);
      }
    }
    
    // Compute Jacobian
    Real jacobian = Real(1);
    for (std::size_t i = 0; i < dim - 1; ++i) {
      Real power = static_cast<Real>(dim - 1 - i);
      jacobian *= std::pow(u[i], power);
    }
    
    return jacobian;
  }
  
  /// \brief Apply inverse Duffy transform
  /// \details Computes u from x where the forward transform is:
  ///   x[0] = u[0] * (1 - u[1]) * (1 - u[2]) * ... * (1 - u[dim-1])
  ///   x[1] = u[0] * u[1] * (1 - u[2]) * ... * (1 - u[dim-1])
  ///   x[2] = u[0] * u[1] * u[2] * (1 - u[3]) * ... * (1 - u[dim-1])
  ///   etc.
  static void inverse(const Real* x, Real* u, std::size_t dim) {
    if (dim == 0) return;
    
    if (dim == 1) {
      u[0] = x[0];
      return;
    }
    
    // Compute sum of all x[i] to get partial products
    Real sum = Real(0);
    for (std::size_t i = 0; i < dim; ++i) {
      sum += x[i];
    }
    
    // u[0] is the sum of all x components
    u[0] = sum;
    
    // For higher dimensions, compute u[i] recursively
    if (std::abs(sum) > std::numeric_limits<Real>::epsilon()) {
      // Compute partial sums for reconstruction
      std::vector<Real> partial_sum(dim);
      partial_sum[dim-1] = x[dim-1];
      for (int i = dim-2; i >= 0; --i) {
        partial_sum[i] = partial_sum[i+1] + x[i];
      }
      
      // Compute u[i] values
      for (std::size_t i = 1; i < dim; ++i) {
        if (std::abs(partial_sum[i-1]) > std::numeric_limits<Real>::epsilon()) {
          u[i] = partial_sum[i] / partial_sum[i-1];
        } else {
          u[i] = Real(0);
        }
      }
    } else {
      // Degenerate case
      for (std::size_t i = 1; i < dim; ++i) {
        u[i] = Real(0);
      }
    }
  }
};

/// \brief Transform wrapper for infinite domain integrands
template <typename Real, typename F, typename Transform>
class transformed_integrand {
private:
  const F& integrand_function_;
  const Transform& coordinate_transform_;
  std::size_t dimension_;
  std::vector<std::size_t> transformed_dimensions_;
  
public:
  transformed_integrand(const F& f, const Transform& t, 
                       std::size_t dim,
                       const std::vector<std::size_t>& inf_dims = {})
    : integrand_function_(f), coordinate_transform_(t), dimension_(dim), transformed_dimensions_(inf_dims) {
    // If no dimensions specified, transform all
    if (transformed_dimensions_.empty()) {
      transformed_dimensions_.resize(dim);
      std::iota(transformed_dimensions_.begin(), transformed_dimensions_.end(), 0);
    }
  }
  
  Real operator()(const Real* u) const {
    std::vector<Real> x(dimension_);
    Real jacobian = Real(1);
    
    // Copy input
    std::copy(u, u + dimension_, x.begin());
    
    // Apply transform to specified dimensions
    for (std::size_t i : transformed_dimensions_) {
      auto result = coordinate_transform_.forward(u[i]);
      x[i] = result.first;
      jacobian *= result.second;
    }
    
    // Evaluate original integrand with Jacobian
    return integrand_function_(x.data(), dimension_) * jacobian;
  }
  
  // Support other integrand signatures
  Real operator()(const Real* u, std::size_t /* dimension */) const {
    return (*this)(u);
  }
  
  Real operator()(const std::vector<Real>& u) const {
    return (*this)(u.data());
  }
};

/// \brief Helper to create transformed integrand
template <typename Real, typename F, typename Transform>
auto make_transformed_integrand(const F& f, const Transform& t, 
                               std::size_t dim,
                               const std::vector<std::size_t>& inf_dims = {}) {
  return transformed_integrand<Real, F, Transform>(f, t, dim, inf_dims);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_TRANSFORMS_HPP

