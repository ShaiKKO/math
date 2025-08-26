#ifndef BOOST_MATH_CUBATURE_POLICIES_HPP
#define BOOST_MATH_CUBATURE_POLICIES_HPP

// Experimental toggle for this preview module
#ifndef BOOST_MATH_CUBATURE_EXPERIMENTAL
#define BOOST_MATH_CUBATURE_EXPERIMENTAL 1
#endif

#include <boost/math/policies/policy.hpp>
#include <cstdint>

namespace boost { namespace math { namespace cubature {

/// \file policies.hpp
/// \brief Policy aliasing and result/status types for cubature front-ends.

/// Status codes for integration routines
enum class status_code { success, maxeval_reached, cancelled };

/// Result type returned by all algorithms
/// - value: integral estimate
/// - error: absolute error bound or standard error (QMC)
/// - evaluations: number of integrand evaluations performed
/// - status: status_code indicating outcome
template <class Real>
struct result {
  Real value{};
  Real error{};
  std::uint64_t evaluations{0};
  status_code status{status_code::success};
};

/// Policy alias forwarded from Boost.Math policies
using default_policy = boost::math::policies::policy<>;

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_POLICIES_HPP

