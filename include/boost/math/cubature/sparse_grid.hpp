// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_SPARSE_GRID_HPP
#define BOOST_MATH_CUBATURE_SPARSE_GRID_HPP

#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/concepts.hpp>
#include <boost/math/cubature/workspace.hpp>
#include <boost/math/cubature/detail/sparse_grid_impl.hpp>

namespace boost { namespace math { namespace cubature {

/// \file sparse_grid.hpp
/// \brief Smolyak sparse-grid cubature public API skeleton.
/// \details Implements nested Clenshaw–Curtis (default) or Gauss–Patterson 1-D
///  rules with node de-duplication; returns result<Real> with simple error proxy.
///
/// \warning Sparse grids can struggle with highly localized functions (sharp peaks).
///          The Smolyak construction uses alternating signs (inclusion-exclusion principle)
///          which can lead to negative weights. For functions with sharp peaks away from
///          the Clenshaw-Curtis nodes, this can produce incorrect results including
///          negative values for positive integrands. Consider using adaptive or QMC
///          integration for such functions.
///
/// \note Performance characteristics:
///       - Excellent for smooth functions with bounded mixed derivatives
///       - Polynomial exactness of degree 2k-1 for level k
///       - Curse of dimensionality is mitigated but not eliminated
///       - Best suited for moderate dimensions (d < 10) and smooth integrands

/// \brief Integrate using Smolyak sparse grid
/// \param f Integrand function
/// \param box Integration domain  
/// \param level Sparse grid level (controls accuracy)
/// \param pol Boost.Math policy
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_sparse_grid(const F& f, const hypercube<Real>& box,
                                          unsigned level, Policy const& /*pol*/ = Policy{})
{
  return detail::integrate_sparse_grid_impl<Real>(f, box, level);
}

/// \brief Integrate with execution options and optional workspace
/// \param ex Execution options for parallelism and memory control
/// \param ws Optional workspace to reuse allocations across calls
/// \note Other parameters as above

template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_sparse_grid(const F& f, const hypercube<Real>& box,
                                          unsigned level,
                                          execution_options const& /*ex*/,
                                          workspace<Real>* /*ws*/ = nullptr,
                                          Policy const& pol = Policy{})
{
  return integrate_sparse_grid<Real, F, Policy>(f, box, level, pol);
}

/// \brief Integrate vector-valued function using Smolyak sparse grid
/// \param f Vector integrand function (const Real*, Real*, std::size_t)
/// \param box Integration domain
/// \param num_components Number of components in the vector function
/// \param level Sparse grid level (controls accuracy)
/// \param pol Boost.Math policy
template <class Real, class F, class Policy = default_policy>
inline std::vector<result<Real>> integrate_sparse_grid_vector(
    const F& f, const hypercube<Real>& box,
    std::size_t num_components,
    unsigned level, 
    Policy const& /*pol*/ = Policy{})
{
  return detail::integrate_sparse_grid_vector_impl<Real>(f, box, num_components, level);
}

/// \brief Vector integration with execution options and workspace
/// \note Parameters as above for vector-valued functions

template <class Real, class F, class Policy = default_policy>
inline std::vector<result<Real>> integrate_sparse_grid_vector(
    const F& f, const hypercube<Real>& box,
    std::size_t num_components,
    unsigned level,
    execution_options const& /*ex*/,
    workspace<Real>* /*ws*/ = nullptr,
    Policy const& pol = Policy{})
{
  return integrate_sparse_grid_vector<Real, F, Policy>(f, box, num_components, level, pol);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_SPARSE_GRID_HPP

