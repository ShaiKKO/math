// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_QMC_HPP
#define BOOST_MATH_CUBATURE_QMC_HPP

// STL headers first per Boost conventions
#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm>
#include <numeric>

// Boost headers - conditionally include if available
#ifdef BOOST_MATH_HAS_BOOST_RANDOM
#include <boost/random/sobol.hpp>
#endif

// Project headers
#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/detail/sobol_owen.hpp>

namespace boost { namespace math { namespace cubature {

// Forward declarations
namespace detail {
  template <typename Real, typename F, typename Policy>
  result<Real> integrate_qmc_impl(
      const F& f,
      const hypercube<Real>& box,
      std::size_t n_points,
      bool use_scrambling,
      std::size_t n_replicates,
      const Policy& pol);
}

/// \file qmc.hpp
/// \brief Quasi-Monte Carlo integration with optional randomization
/// \details Implements QMC using Sobol sequences with tent transform and
///          optional Owen scrambling for variance estimation
///          Following ALGORITHMS.md section 3 and TECHNICAL_BLUEPRINT.md

/// \brief Integrate using Quasi-Monte Carlo with Sobol sequences
/// \param f Integrand function
/// \param box Integration domain
/// \param n_points Number of QMC points
/// \param pol Boost.Math policy
template <class Real, class F, class Policy = default_policy>
result<Real> integrate_qmc(
    const F& f, 
    const hypercube<Real>& box,
    std::size_t n_points,
    const Policy& pol = Policy{})
{
    return detail::integrate_qmc_impl<Real>(f, box, n_points, false, 1, pol);
}

/// \brief Integrate using Randomized QMC with Owen scrambling
/// \param f Integrand function
/// \param box Integration domain  
/// \param n_points Number of QMC points per replicate
/// \param n_replicates Number of independent replicates for variance estimation
/// \param pol Boost.Math policy
template <class Real, class F, class Policy = default_policy>
result<Real> integrate_rqmc(
    const F& f,
    const hypercube<Real>& box,
    std::size_t n_points,
    std::size_t n_replicates = 10,
    const Policy& pol = Policy{})
{
    return detail::integrate_qmc_impl<Real>(f, box, n_points, true, n_replicates, pol);
}

namespace detail {

/// \brief Tent (Baker) transform for variance reduction
/// \details Maps u to 1 - 2|u - 1/2| to enforce periodicity
template <typename Real>
inline Real tent_transform(Real u) {
    return Real(1) - Real(2) * std::abs(u - Real(0.5));
}

/// \brief QMC implementation with optional randomization
template <typename Real, typename F, typename Policy>
result<Real> integrate_qmc_impl(
    const F& f,
    const hypercube<Real>& box,
    std::size_t n_points,
    bool use_scrambling,
    std::size_t n_replicates,
    const Policy& /*pol*/)
{
#ifndef BOOST_MATH_HAS_BOOST_RANDOM
    // If Boost.Random is not available, return an error
    (void)f; (void)box; (void)n_points; (void)use_scrambling; (void)n_replicates;
    result<Real> res;
    res.status = status_code::not_implemented;
    res.error = std::numeric_limits<Real>::max();
    res.value = 0;
    res.evaluations = 0;
    return res;
#else
    const std::size_t dim = box.dimension();
    
    // Validate inputs
    if (dim == 0 || dim > 40) {  // Sobol supports up to 40 dimensions efficiently
        result<Real> res;
        res.status = status_code::dimension_error;
        res.error = std::numeric_limits<Real>::max();
        return res;
    }
    
    if (n_points == 0) {
        result<Real> res;
        res.status = status_code::invalid_input;
        res.error = std::numeric_limits<Real>::max();
        return res;
    }
    
    // Storage for replicate values (for RQMC variance estimation)
    std::vector<Real> replicate_values;
    if (use_scrambling) {
        replicate_values.reserve(n_replicates);
    }
    
    // Main QMC integration loop
    Real total_sum = 0;
    std::size_t total_evals = 0;
    
    for (std::size_t rep = 0; rep < n_replicates; ++rep) {
        // Create Sobol engine
        boost::random::sobol_engine<std::uint64_t, dim> sobol(dim);
        
        // Apply Owen scrambling if requested
        std::uint32_t scramble_seed = 0;
        if (use_scrambling) {
            // Deterministic seed based on replicate number
            // Use large primes to ensure good seed distribution
            scramble_seed = static_cast<std::uint32_t>(rep * 65521 + 32749);
        }
        
        // Deterministic index partitioning for parallel execution
        // Each thread gets a contiguous block of indices
        std::size_t thread_id = 0;  // Would be omp_get_thread_num() in parallel context
        std::size_t num_threads = 1;  // Would be omp_get_num_threads() in parallel context
        std::size_t points_per_thread = n_points / num_threads;
        std::size_t remainder = n_points % num_threads;
        
        // Calculate this thread's index range
        std::size_t start_idx = thread_id * points_per_thread + std::min(thread_id, remainder);
        std::size_t end_idx = start_idx + points_per_thread + (thread_id < remainder ? 1 : 0);
        
        // Skip Sobol sequence to starting index
        for (std::size_t skip = 0; skip < start_idx; ++skip) {
            sobol();  // Advance without using the value
        }
        
        // Accumulate QMC sum for this replicate
        Real replicate_sum = 0;
        std::vector<Real> point(dim);
        
        for (std::size_t i = start_idx; i < end_idx; ++i) {
            // Generate Sobol point
            auto sobol_point = sobol();
            
            // Convert to [0,1)^d and apply scrambling if needed
            for (std::size_t d = 0; d < dim; ++d) {
                Real u = Real(sobol_point[d]) / Real(std::numeric_limits<std::uint64_t>::max());
                
                if (use_scrambling) {
                    u = detail::owen_scramble(u, scramble_seed + d);
                }
                
                // Apply tent transform for variance reduction
                u = tent_transform(u);
                
                // Map from [0,1] to [a,b]
                point[d] = box.lower[d] + u * (box.upper[d] - box.lower[d]);
            }
            
            // Evaluate integrand - support multiple signatures
            Real value;
            if constexpr (std::is_invocable_v<F, decltype(point)>) {
                // f(std::vector<Real>)
                value = f(point);
            } else if constexpr (std::is_invocable_v<F, const Real*, std::size_t>) {
                // f(const Real*, std::size_t)
                value = f(point.data(), dim);
            } else if constexpr (std::is_invocable_v<F, const Real*>) {
                // f(const Real*) - legacy support
                value = f(point.data());
            } else {
                static_assert(sizeof(F) == 0, "Integrand must be callable with vector<Real>, (const Real*, size_t), or const Real*");
            }
            replicate_sum += value;
            ++total_evals;
        }
        
        // Scale by volume / n_points
        Real volume = 1.0;
        for (std::size_t d = 0; d < dim; ++d) {
            volume *= (box.upper[d] - box.lower[d]);
        }
        replicate_sum = (replicate_sum * volume) / Real(n_points);
        
        if (use_scrambling) {
            replicate_values.push_back(replicate_sum);
        }
        total_sum += replicate_sum;
    }
    
    // Prepare result
    result<Real> res;
    res.evaluations = total_evals;
    res.status = status_code::success;
    
    if (use_scrambling && n_replicates > 1) {
        // RQMC: Use replicate mean and standard error
        res.value = total_sum / Real(n_replicates);
        
        // Compute standard error
        Real mean = res.value;
        Real variance = 0;
        for (const auto& val : replicate_values) {
            Real diff = val - mean;
            variance += diff * diff;
        }
        variance /= Real(n_replicates - 1);
        
        // Standard error of the mean
        Real std_error = std::sqrt(variance / Real(n_replicates));
        
        // Use t-distribution critical value for 95% CI (approximate as 2 for simplicity)
        res.error = Real(2) * std_error;
    } else {
        // Plain QMC: No error estimate available
        res.value = total_sum / Real(n_replicates);
        res.error = 0;  // Cannot estimate error without randomization
    }
    
    return res;
#endif  // BOOST_MATH_HAS_BOOST_RANDOM
}

} // namespace detail

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_QMC_HPP