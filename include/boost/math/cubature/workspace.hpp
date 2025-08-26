#ifndef BOOST_MATH_CUBATURE_WORKSPACE_HPP
#define BOOST_MATH_CUBATURE_WORKSPACE_HPP

#include <boost/math/cubature/policies.hpp>
#include <cstddef>

namespace boost { namespace math { namespace cubature {

struct execution_options {
  std::size_t max_threads = 0; // 0 -> use hardware_concurrency at call-site
  std::uint64_t max_eval  = 0;
  bool deterministic      = true;
};

// Placeholder workspace that will later contain scratch buffers and state
struct workspace {
  // TODO: add buffers for abscissas, weights, PQ, QMC state, etc.
};

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_WORKSPACE_HPP

