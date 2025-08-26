#ifndef BOOST_MATH_CUBATURE_QMC_HPP
#define BOOST_MATH_CUBATURE_QMC_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/concepts.hpp>
#include <boost/math/cubature/workspace.hpp>

#ifndef BOOST_MATH_STANDALONE
#include <boost/random/sobol.hpp>
#endif

namespace boost { namespace math { namespace cubature {

/// \file qmc.hpp
/// \brief Quasi-Monte-Carlo (Sobol') / RQMC front-end public API skeleton.
/// \details Uses Boost.Random sobol_engine when available; result<Real>::error
///  reports replicate-based standard error across k replicates.

#ifdef BOOST_MATH_STANDALONE
template <class Real, class F, class Policy = default_policy, class Sobol = void>
#else
template <class Real, class F, class Policy = default_policy, class Sobol = boost::random::sobol_engine<>>
#endif
inline result<Real> integrate_qmc(const F& /*f*/, const hypercube<Real>& /*box*/,
                                  std::size_t /*n_points_per_rep*/, unsigned /*replicates*/ = 1,
                                  bool /*tent_transform*/ = true,
                                  Policy const& /*pol*/ = Policy{})
{
  return result<Real>{};
}

// Overload with execution options and workspace
#ifdef BOOST_MATH_STANDALONE
template <class Real, class F, class Policy = default_policy, class Sobol = void>
#else
template <class Real, class F, class Policy = default_policy, class Sobol = boost::random::sobol_engine<>>
#endif
inline result<Real> integrate_qmc(const F& f, const hypercube<Real>& box,
                                  std::size_t n_points_per_rep, unsigned replicates,
                                  bool tent_transform,
                                  execution_options const& /*ex*/,
                                  workspace* /*ws*/ = nullptr,
                                  Policy const& pol = Policy{})
{
  return integrate_qmc<Real, F, Policy, Sobol>(f, box, n_points_per_rep, replicates, tent_transform, pol);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_QMC_HPP

