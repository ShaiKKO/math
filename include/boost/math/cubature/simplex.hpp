#ifndef BOOST_MATH_CUBATURE_SIMPLEX_HPP
#define BOOST_MATH_CUBATURE_SIMPLEX_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/concepts.hpp>
#include <boost/math/cubature/workspace.hpp>

namespace boost { namespace math { namespace cubature {

/// \file simplex.hpp
/// \brief Simplex-domain integrator front-end using Duffy mapping + backend.
/// \details Accepts a general simplex (d+1 vertices of length d). The initial
///  implementation will map to a unit hypercube and call the adaptive/QMC backend.

template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_simplex(const F& /*f*/, const simplex<Real>& /*simp*/,
                                      Real /*abs_tol*/, Real /*rel_tol*/,
                                      Policy const& /*pol*/ = Policy{})
{
  return result<Real>{};
}

// Overload with execution options

template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_simplex(const F& f, const simplex<Real>& simp,
                                      Real abs_tol, Real rel_tol,
                                      execution_options const& /*ex*/,
                                      Policy const& pol = Policy{})
{
  return integrate_simplex<Real, F, Policy>(f, simp, abs_tol, rel_tol, pol);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_SIMPLEX_HPP

