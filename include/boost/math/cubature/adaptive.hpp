#ifndef BOOST_MATH_CUBATURE_ADAPTIVE_HPP
#define BOOST_MATH_CUBATURE_ADAPTIVE_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/concepts.hpp>
#include <boost/math/cubature/workspace.hpp>

namespace boost { namespace math { namespace cubature {

/// \file adaptive.hpp
/// \brief Adaptive hypercube cubature (DCUHRE/Cuhre-style) public API skeleton.
/// \details See ALGORITHMS.md §1 and ARCHITECTURE.md §2. Returns result<Real> with
///  value (estimate), error (embedded-difference bound), evaluations, and status.

// Adaptive cubature (DCUHRE-style): API skeleton only

template <class Real, class F, class Policy = default_policy
#if !__cpp_concepts
  , typename = typename std::enable_if<is_scalar_integrand<F, Real>::value || is_vector_integrand<F, Real>::value>::type
#endif
>
#if __cpp_concepts
  requires (ScalarIntegrand<F, Real> || VectorIntegrand<F, Real>)
#endif
inline result<Real> integrate_adaptive(const F& /*f*/, const hypercube<Real>& /*box*/,
                                       Real /*abs_tol*/, Real /*rel_tol*/,
                                       std::size_t /*max_eval*/ = 0,
                                       Policy const& /*pol*/ = Policy{})
{
  return result<Real>{};
}

// Overload with execution options

template <class Real, class F, class Policy = default_policy
#if !__cpp_concepts
  , typename = typename std::enable_if<is_scalar_integrand<F, Real>::value || is_vector_integrand<F, Real>::value>::type
#endif
>
#if __cpp_concepts
  requires (ScalarIntegrand<F, Real> || VectorIntegrand<F, Real>)
#endif
inline result<Real> integrate_adaptive(const F& f, const hypercube<Real>& box,
                                       Real abs_tol, Real rel_tol,
                                       execution_options const& /*ex*/,
                                       std::size_t max_eval = 0,
                                       Policy const& pol = Policy{})
{
  return integrate_adaptive<Real, F, Policy>(f, box, abs_tol, rel_tol, max_eval, pol);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_ADAPTIVE_HPP

