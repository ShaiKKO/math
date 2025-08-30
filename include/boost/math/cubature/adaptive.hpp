// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_ADAPTIVE_HPP
#define BOOST_MATH_CUBATURE_ADAPTIVE_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/concepts.hpp>
#include <boost/math/cubature/workspace.hpp>
#include <boost/math/cubature/detail/adaptive_impl.hpp>
#include <boost/math/cubature/detail/adaptive_integrator_vector.hpp>
#include <boost/math/cubature/detail/vector_adapter.hpp>
#include <vector>

namespace boost { namespace math { namespace cubature {

/// \file adaptive.hpp
/// \brief Adaptive integration with embedded error estimation
/// \details Implements the DCUHRE (Double Cone Uniform Hypercube Rule Enhancement)
///          algorithm for adaptive multidimensional integration. Uses embedded
///          Genz-Malik rules for error estimation and axis-aligned subdivision.
///
/// Key features:
///  - Global adaptive refinement (always subdivides region with largest error)
///  - Embedded rule pairs (7/5 or 9/7) for reliable error estimation
///  - Fourth-difference based axis selection for optimal subdivision
///  - Safeguarded error estimates to prevent spurious convergence
///  - Support for vector-valued integrands
///
/// References:
///  - Berntsen, Espelid, Genz (1991): "An adaptive algorithm for the approximate
///    calculation of multiple integrals"
///  - Genz, Malik (1980): "An adaptive algorithm for numerical integration over
///    an N-dimensional rectangular region"
///
/// \author Colin MacRitchie
/// \date 2025

/// \brief Integrate using adaptive subdivision with embedded error estimation
/// \tparam Real Floating-point type (float, double, long double, or multiprecision)
/// \tparam F Integrand type callable as f(const std::vector<Real>&) -> Real
/// \tparam Policy Boost.Math policy type for error handling
/// \param f Integrand function accepting std::vector<Real> and returning Real
/// \param box Integration domain (hypercube)
/// \param abs_tol Absolute error tolerance
/// \param rel_tol Relative error tolerance
/// \param max_eval Maximum number of function evaluations (0 = automatic)
/// \param pol Boost.Math policy object
/// \return result<Real> containing integral estimate, error bound, and diagnostics
/// \note Dimensions 2-15 are supported with optimized Genz-Malik rules
template <class Real, class F, class Policy = default_policy
#if !__cpp_concepts
  , typename = typename std::enable_if<is_scalar_integrand<F, Real>::value>::type
#endif
>
#if __cpp_concepts
  requires ScalarIntegrand<F, Real>
#endif
inline result<Real> integrate_adaptive(const F& f, const hypercube<Real>& box,
                                       Real abs_tol, Real rel_tol,
                                       std::size_t max_eval = 0,
                                       Policy const& pol = Policy{})
{
  // Create integrator and run
  detail::adaptive_integrator<Real, F, Policy> integrator(
    f, box, abs_tol, rel_tol, max_eval, pol);
  
  return integrator.integrate();
}

/// \brief Integrate with explicit execution options
/// \param ex Execution options controlling parallelism and memory usage
/// \note Other parameters as above
template <class Real, class F, class Policy = default_policy
#if !__cpp_concepts
  , typename = typename std::enable_if<is_scalar_integrand<F, Real>::value>::type
#endif
>
#if __cpp_concepts
  requires ScalarIntegrand<F, Real>
#endif
inline result<Real> integrate_adaptive(const F& f, const hypercube<Real>& box,
                                       Real abs_tol, Real rel_tol,
                                       execution_options const& ex,
                                       std::size_t max_eval = 0,
                                       Policy const& pol = Policy{})
{
  // For now, ignore execution options (parallelism will be added in P09)
  // Just forward to main implementation
  return integrate_adaptive<Real, F, Policy>(f, box, abs_tol, rel_tol, 
                                             max_eval > 0 ? max_eval : ex.max_eval, pol);
}

// Overload with workspace (for amortizing allocations)
template <class Real, class F, class Policy = default_policy
#if !__cpp_concepts
  , typename = typename std::enable_if<is_scalar_integrand<F, Real>::value>::type
#endif
>
#if __cpp_concepts
  requires ScalarIntegrand<F, Real>
#endif
inline result<Real> integrate_adaptive(const F& f, const hypercube<Real>& box,
                                       Real abs_tol, Real rel_tol,
                                       workspace<Real>& /*ws*/,
                                       std::size_t max_eval = 0,
                                       Policy const& pol = Policy{})
{
  // Workspace support will be added in P09
  // For now, just forward to main implementation
  return integrate_adaptive<Real, F, Policy>(f, box, abs_tol, rel_tol, max_eval, pol);
}

// Vector-valued integrand support
template <class Real, class F, class Policy = default_policy
#if !__cpp_concepts
  , typename = typename std::enable_if<is_vector_integrand<F, Real>::value>::type
#endif
>
#if __cpp_concepts
  requires VectorIntegrand<F, Real>
#endif
inline std::vector<result<Real>> integrate_adaptive_vector(
    const F& f, const hypercube<Real>& box,
    std::size_t num_components,
    Real abs_tol, Real rel_tol,
    std::size_t max_eval = 0,
    detail::error_norm norm = detail::error_norm::l_infinity,
    Policy const& pol = Policy{})
{
  // Create vector integrator and run
  detail::adaptive_integrator_vector<Real, F, Policy> integrator(
    f, box, abs_tol, rel_tol, num_components, max_eval, norm, pol);
  
  return integrator.integrate();
}

// Convenience overload for vector integrands with automatic dispatch
template <class Real, class F, class Policy = default_policy
#if !__cpp_concepts
  , typename = typename std::enable_if<is_vector_integrand<F, Real>::value>::type
#endif
>
#if __cpp_concepts
  requires VectorIntegrand<F, Real>
#endif
inline std::vector<result<Real>> integrate_adaptive(
    const F& f, const hypercube<Real>& box,
    std::size_t num_components,
    Real abs_tol, Real rel_tol,
    std::size_t max_eval = 0,
    Policy const& pol = Policy{})
{
  return integrate_adaptive_vector<Real, F, Policy>(
    f, box, num_components, abs_tol, rel_tol, max_eval, 
    detail::error_norm::l_infinity, pol);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_ADAPTIVE_HPP