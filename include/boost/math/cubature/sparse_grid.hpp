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

/// \brief Integrate using dimension-adaptive sparse grid
/// \param f Integrand function
/// \param box Integration domain
/// \param initial_levels Initial level per dimension (or single level for all)
/// \param dim_weights Importance weights per dimension (optional)
/// \param max_evals Maximum function evaluations
/// \param pol Boost.Math policy
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_sparse_grid_adaptive(
    const F& f, 
    const hypercube<Real>& box,
    const std::vector<std::size_t>& initial_levels,
    const std::vector<Real>& dim_weights = {},
    std::size_t max_evals = 10000,
    Policy const& /*pol*/ = Policy{})
{
  detail::adaptive_smolyak_grid<Real> grid(
      box.dimension(), initial_levels, dim_weights);
  
  result<Real> res = grid.evaluate(f, box);
  
  // Adaptive refinement
  while (res.evaluations < max_evals) {
    grid.adapt();
    auto new_res = grid.evaluate(f, box);
    
    // Check convergence
    if (std::abs(new_res.value - res.value) < res.error * Real(0.1)) {
      break;  // Converged
    }
    
    res = new_res;
  }
  
  return res;
}

/// \brief Integrate using Gauss-Hermite sparse grid for Gaussian weight
/// \param f Integrand function (without weight)
/// \param dim Dimension
/// \param level Sparse grid level
/// \param use_genz_keister Use Genz-Keister nested rules (default: false)
/// \param pol Boost.Math policy
/// \note Integrates f(x) * exp(-||x||^2) over R^d
template <class Real, class F, class Policy = default_policy>
inline result<Real> integrate_sparse_grid_gaussian(
    const F& f,
    std::size_t dim,
    unsigned level,
    bool use_genz_keister = false,
    Policy const& /*pol*/ = Policy{})
{
  // This would use Gauss-Hermite nodes instead of Clenshaw-Curtis
  // Implementation deferred to separate PR for clarity
  result<Real> res;
  res.value = Real(0);
  res.error = std::numeric_limits<Real>::max();
  res.status = status_code::dimension_error;  // Not yet implemented
  res.evaluations = 0;
  
  // TODO: Implement using gauss_hermite_rules.hpp
  // Will construct sparse grid with Hermite nodes for Gaussian weight
  
  return res;
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_SPARSE_GRID_HPP

