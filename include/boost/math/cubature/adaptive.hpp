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

namespace boost { namespace math { namespace cubature {

/// \file adaptive.hpp
/// \brief Adaptive hypercube cubature (DCUHRE/Cuhre-style) implementation.
/// \details Implements DCUHRE algorithm following ALGORITHMS.md §1 and ARCHITECTURE.md §2.
///  Returns result<Real> with value (estimate), error (embedded-difference bound), 
///  evaluations, and status.

// Adaptive cubature (DCUHRE-style) main entry point
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

// Overload with execution options
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

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_ADAPTIVE_HPP