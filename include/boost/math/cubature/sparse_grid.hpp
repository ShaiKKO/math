#ifndef BOOST_MATH_CUBATURE_SPARSE_GRID_HPP
#define BOOST_MATH_CUBATURE_SPARSE_GRID_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/concepts.hpp>
#include <boost/math/cubature/workspace.hpp>

namespace boost { namespace math { namespace cubature {

/// \file sparse_grid.hpp
/// \brief Smolyak sparse-grid cubature public API skeleton.
/// \details Implements nested Clenshaw–Curtis (default) or Gauss–Patterson 1-D
///  rules with node de-duplication; returns result<Real> with simple error proxy.

// Sparse-grid Smolyak: API skeleton only

template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_sparse_grid(const F& /*f*/, const hypercube<Real>& /*box*/,
                                          unsigned /*level*/, Policy const& /*pol*/ = Policy{})
{
  return result<Real>{};
}

// Overload with execution options and workspace

template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_sparse_grid(const F& f, const hypercube<Real>& box,
                                          unsigned level,
                                          execution_options const& /*ex*/,
                                          workspace* /*ws*/ = nullptr,
                                          Policy const& pol = Policy{})
{
  return integrate_sparse_grid<Real, F, Policy>(f, box, level, pol);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_SPARSE_GRID_HPP

