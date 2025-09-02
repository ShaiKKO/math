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
#include <type_traits>

// Check for Boost.Random availability
#ifndef BOOST_MATH_HAS_BOOST_RANDOM
#ifdef __has_include
#if __has_include(<boost/random/sobol.hpp>)
#define BOOST_MATH_HAS_BOOST_RANDOM 1
#endif
#else
// Assume Boost.Random is available if we can't check
#define BOOST_MATH_HAS_BOOST_RANDOM 1
#endif
#endif

// Boost headers - conditionally include if available
#ifdef BOOST_MATH_HAS_BOOST_RANDOM
#include <boost/random/sobol.hpp>
#endif

// Project headers
#include <boost/math/cubature/policies.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/detail/sobol_owen.hpp>
#include <boost/math/cubature/detail/qmc_sequences.hpp>

namespace boost { namespace math { namespace cubature {

/// \brief Options for QMC integration
/// \details Controls sequence type and other parameters for QMC integration
template <typename Real = double>
struct qmc_options {
    /// \brief Type of QMC sequence to use
    detail::qmc_sequence_type sequence_type = detail::qmc_sequence_type::sobol;
    
    /// \brief Number of points for lattice rules (ignored for other sequences)
    std::size_t lattice_points = 0;
    
    /// \brief Use scrambling for randomized QMC
    bool use_scrambling = false;
    
    /// \brief Number of independent replicates for RQMC
    std::size_t n_replicates = 10;
    
    /// \brief Random seed for scrambling
    std::uint32_t scramble_seed = 0;
    
    /// \brief Apply tent transform for variance reduction
    bool use_tent_transform = true;
    
    /// \brief Leap value for Halton sequence (1 = standard, >1 = leaped)
    std::size_t halton_leap = 1;
    
    /// \brief Default constructor
    qmc_options() = default;
    
    /// \brief Constructor for backward compatibility (Sobol-only)
    qmc_options(bool scrambling, std::size_t replicates = 10)
        : use_scrambling(scrambling), n_replicates(replicates) {}
};

/// \file qmc.hpp
/// \brief Quasi-Monte Carlo integration with optional randomization
/// \details Implements (Randomized) Quasi-Monte Carlo integration using:
///          - Sobol sequences for low-discrepancy point generation
///          - Owen scrambling for randomization (RQMC)
///          - Tent (Baker) transformation for variance reduction
///          - Statistical error estimation via independent replicates
///
/// The QMC method is particularly effective for:
///  - High-dimensional integrals (d > 10)
///  - Smooth integrands with low effective dimension
///  - Problems where Monte Carlo convergence is too slow
///
/// Convergence rate: O(n^{-1}(log n)^d) for QMC vs O(n^{-1/2}) for Monte Carlo
///
/// References:
///  - Owen, A.B. (2003). "Variance with alternative scramblings of digital nets"
///  - Dick, J. & Pillichshammer, F. (2010). "Digital Nets and Sequences"
///  - Sobol, I.M. (1967). "Distribution of points in a cube and approximate 
///    evaluation of integrals"

// Forward declarations for implementation details
namespace detail {
  template <typename Real, typename F, typename Policy>
  result<Real> integrate_qmc_impl(
      const F& f,
      const hypercube<Real>& box,
      std::size_t n_points,
      bool use_scrambling,
      std::size_t n_replicates,
      const Policy& pol);
      
  template <typename Real, typename F, typename Policy>
  result<Real> integrate_qmc_impl_ex(
      const F& f,
      const hypercube<Real>& box,
      std::size_t n_points,
      const qmc_options<Real>& options,
      const Policy& pol);
      
  template <typename Real, typename F, typename Policy>
  std::vector<result<Real>> integrate_qmc_vector_impl(
      const F& f,
      const hypercube<Real>& box,
      std::size_t num_components,
      std::size_t n_points,
      bool use_scrambling,
      std::size_t n_replicates,
      const Policy& pol);
}

/// \brief Integrate using Quasi-Monte Carlo with Sobol sequences
/// \details Uses low-discrepancy Sobol sequences to achieve faster convergence
///          than Monte Carlo for smooth integrands. The tent transform is
///          applied to improve performance for integrands with boundary peaks.
///
/// \tparam Real Floating-point type
/// \tparam F Integrand type callable as f(const Real*, std::size_t) -> Real
/// \tparam Policy Boost.Math policy type
///
/// \param f Integrand function
/// \param box Integration domain (hypercube)
/// \param n_points Number of QMC points to use
/// \param pol Boost.Math policy object
///
/// \return result<Real> with integral estimate (no error estimate for plain QMC)
///
/// \note Error estimate is zero for plain QMC. Use integrate_rqmc for error bars.
template <class Real, class F, class Policy = default_policy>
result<Real> integrate_qmc(
    const F& f, 
    const hypercube<Real>& box,
    std::size_t n_points,
    const Policy& pol = Policy{})
{
    // Use the new unified implementation via qmc_options
    qmc_options<Real> opts;
    opts.sequence_type = detail::qmc_sequence_type::sobol;
    opts.use_scrambling = false;
    opts.n_replicates = 1;
    return detail::integrate_qmc_impl_ex<Real>(f, box, n_points, opts, pol);
}

/// \brief Integrate using Randomized Quasi-Monte Carlo (RQMC)
/// \details Applies Owen scrambling to Sobol sequences to enable statistical
///          error estimation while preserving the superior convergence of QMC.
///          Each replicate uses an independently scrambled sequence.
///
/// \tparam Real Floating-point type
/// \tparam F Integrand type callable as f(const Real*, std::size_t) -> Real
/// \tparam Policy Boost.Math policy type
///
/// \param f Integrand function
/// \param box Integration domain (hypercube)
/// \param n_points Number of QMC points per replicate
/// \param n_replicates Number of independent replicates (typically 10-30)
/// \param pol Boost.Math policy object
///
/// \return result<Real> with integral estimate and statistical error (2σ)
///
/// \note The error estimate is computed as 2 * standard error of the replicate means,
///       providing approximately 95% confidence interval for normally distributed errors.
/// \param pol Boost.Math policy
template <class Real, class F, class Policy = default_policy>
result<Real> integrate_rqmc(
    const F& f,
    const hypercube<Real>& box,
    std::size_t n_points,
    std::size_t n_replicates = 10,
    const Policy& pol = Policy{})
{
    // Use the new unified implementation via qmc_options
    qmc_options<Real> opts;
    opts.sequence_type = detail::qmc_sequence_type::sobol;
    opts.use_scrambling = true;
    opts.n_replicates = n_replicates;
    return detail::integrate_qmc_impl_ex<Real>(f, box, n_points, opts, pol);
}

/// \brief Integrate using QMC with configurable sequence type
/// \details Uses specified QMC sequence type for integration with full control
///          over all parameters including sequence type, scrambling, and replicates.
///
/// \tparam Real Floating-point type
/// \tparam F Integrand type callable as f(const Real*, std::size_t) -> Real
/// \tparam Policy Boost.Math policy type
///
/// \param f Integrand function
/// \param box Integration domain (hypercube)
/// \param n_points Number of QMC points (per replicate if using RQMC)
/// \param options QMC options controlling sequence type and parameters
/// \param pol Boost.Math policy object
///
/// \return result<Real> with integral estimate and error (if using RQMC)
template <class Real, class F, class Policy = default_policy>
result<Real> integrate_qmc_ex(
    const F& f,
    const hypercube<Real>& box,
    std::size_t n_points,
    const qmc_options<Real>& options,
    const Policy& pol = Policy{})
{
    return detail::integrate_qmc_impl_ex<Real>(f, box, n_points, options, pol);
}

/// \brief Integrate vector-valued function using Quasi-Monte Carlo
/// \param f Vector integrand function (const Real*, Real*, std::size_t)
/// \param box Integration domain
/// \param num_components Number of components in the vector function
/// \param n_points Number of QMC points
/// \param pol Boost.Math policy
template <class Real, class F, class Policy = default_policy>
std::vector<result<Real>> integrate_qmc_vector(
    const F& f,
    const hypercube<Real>& box,
    std::size_t num_components,
    std::size_t n_points,
    const Policy& pol = Policy{})
{
    return detail::integrate_qmc_vector_impl<Real>(f, box, num_components, n_points, false, 1, pol);
}

/// \brief Integrate vector-valued function using Randomized QMC
/// \param f Vector integrand function (const Real*, Real*, std::size_t)
/// \param box Integration domain
/// \param num_components Number of components in the vector function
/// \param n_points Number of QMC points per replicate
/// \param n_replicates Number of independent replicates for variance estimation
/// \param pol Boost.Math policy
template <class Real, class F, class Policy = default_policy>
std::vector<result<Real>> integrate_rqmc_vector(
    const F& f,
    const hypercube<Real>& box,
    std::size_t num_components,
    std::size_t n_points,
    std::size_t n_replicates = 10,
    const Policy& pol = Policy{})
{
    return detail::integrate_qmc_vector_impl<Real>(f, box, num_components, n_points, true, n_replicates, pol);
}

namespace detail {

/// \brief Tent (Baker) transform for variance reduction
/// \details Maps u to 1 - 2|u - 1/2| to enforce periodicity
template <typename Real>
inline Real tent_transform(Real u) {
    return Real(1) - Real(2) * std::abs(u - Real(0.5));
}

// SFINAE helper for detecting integrand signature
template <typename F, typename Real, typename = void>
struct has_pointer_size_signature : std::false_type {};

template <typename F, typename Real>
struct has_pointer_size_signature<F, Real,
    std::void_t<decltype(std::declval<F>()(std::declval<const Real*>(), std::declval<std::size_t>()))>>
    : std::true_type {};

template <typename F, typename Real, typename = void>
struct has_pointer_signature : std::false_type {};

template <typename F, typename Real>
struct has_pointer_signature<F, Real,
    std::void_t<decltype(std::declval<F>()(std::declval<const Real*>()))>>
    : std::true_type {};

template <typename F, typename Point, typename = void>
struct has_vector_signature : std::false_type {};

template <typename F, typename Point>
struct has_vector_signature<F, Point,
    std::void_t<decltype(std::declval<F>()(std::declval<Point>()))>>
    : std::true_type {};

// Unified integrand evaluation
template <typename F, typename Point>
inline auto evaluate_integrand_impl(const F& f, const Point& point) 
    -> typename Point::value_type 
{
    using Real = typename Point::value_type;
    
    // Try pointer + size signature first (most common)
    if constexpr (has_pointer_size_signature<F, Real>::value) {
        return f(point.data(), point.size());
    }
    // Try pointer-only signature
    else if constexpr (has_pointer_signature<F, Real>::value) {
        return f(point.data());
    }
    // Try vector/array signature
    else if constexpr (has_vector_signature<F, Point>::value) {
        return f(point);
    }
    else {
        static_assert(has_pointer_size_signature<F, Real>::value ||
                     has_pointer_signature<F, Real>::value ||
                     has_vector_signature<F, Point>::value,
                     "Integrand must be callable with (const Real*, size_t), (const Real*), or vector<Real>");
        return Real{};
    }
}

// Helper to dispatch to correct Sobol dimension at compile time
template <typename Real, typename F, std::size_t Dim>
result<Real> integrate_qmc_impl_fixed_dim(
    const F& f,
    const hypercube<Real>& box,
    std::size_t n_points,
    bool use_scrambling,
    std::size_t n_replicates)
{
#ifdef BOOST_MATH_HAS_BOOST_RANDOM
    // Storage for replicate values (for RQMC variance estimation)
    std::vector<Real> replicate_values;
    if (use_scrambling) {
        replicate_values.reserve(n_replicates);
    }
    
    // Main QMC integration loop
    Real total_sum = 0;
    std::size_t total_evals = 0;
    
    for (std::size_t rep = 0; rep < n_replicates; ++rep) {
        // Create Sobol engine with compile-time dimension
        boost::random::sobol_engine<std::uint64_t, Dim> sobol(0);
        
        // Apply Owen scrambling if requested
        std::uint32_t scramble_seed = 0;
        if (use_scrambling) {
            scramble_seed = static_cast<std::uint32_t>(rep * 65521 + 32749);
        }
        
        // Accumulate QMC sum for this replicate
        Real replicate_sum = 0;
        std::vector<Real> point(Dim);
        
        for (std::size_t i = 0; i < n_points; ++i) {
            // Generate Sobol point
            auto sobol_point = sobol();
            
            // Convert to [0,1)^d and apply scrambling if needed
            for (std::size_t d = 0; d < Dim; ++d) {
                Real u = Real(sobol_point[d]) / Real(std::numeric_limits<std::uint64_t>::max());
                
                if (use_scrambling) {
                    u = detail::owen_scramble(u, scramble_seed + static_cast<std::uint32_t>(d));
                }
                
                // Apply tent transform for variance reduction
                u = tent_transform(u);
                
                // Map from [0,1] to [a,b]
                point[d] = box.lower[d] + u * (box.upper[d] - box.lower[d]);
            }
            
            // Evaluate integrand
            Real value = evaluate_integrand_impl(f, point);
            replicate_sum += value;
            ++total_evals;
        }
        
        // Scale by volume / n_points
        Real volume = 1.0;
        for (std::size_t d = 0; d < Dim; ++d) {
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
        res.value = total_sum / Real(n_replicates);
        
        // Compute standard error
        Real mean = res.value;
        Real variance = 0;
        for (const auto& val : replicate_values) {
            Real diff = val - mean;
            variance += diff * diff;
        }
        variance /= Real(n_replicates - 1);
        res.error = Real(2) * std::sqrt(variance / Real(n_replicates));
    } else {
        res.value = total_sum / Real(n_replicates);
        res.error = 0;
    }
    
    return res;
#else
    (void)f; (void)box; (void)n_points; (void)use_scrambling; (void)n_replicates;
    result<Real> res;
    res.status = status_code::not_implemented;
    res.error = std::numeric_limits<Real>::max();
    res.value = 0;
    res.evaluations = 0;
    return res;
#endif
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
    const std::size_t dim = box.dimension();
    
    // Validate inputs
    if (dim == 0 || dim > 20) {  // Limit to 20 dimensions for practical dispatch
        result<Real> res;
        res.status = status_code::dimension_error;
        res.error = std::numeric_limits<Real>::max();
        res.value = 0;
        res.evaluations = 0;
        return res;
    }
    
    if (n_points == 0) {
        result<Real> res;
        res.status = status_code::invalid_input;
        res.error = std::numeric_limits<Real>::max();
        res.value = 0;
        res.evaluations = 0;
        return res;
    }
    
    // Dispatch to compile-time dimension
    switch(dim) {
        case 1: return integrate_qmc_impl_fixed_dim<Real, F, 1>(f, box, n_points, use_scrambling, n_replicates);
        case 2: return integrate_qmc_impl_fixed_dim<Real, F, 2>(f, box, n_points, use_scrambling, n_replicates);
        case 3: return integrate_qmc_impl_fixed_dim<Real, F, 3>(f, box, n_points, use_scrambling, n_replicates);
        case 4: return integrate_qmc_impl_fixed_dim<Real, F, 4>(f, box, n_points, use_scrambling, n_replicates);
        case 5: return integrate_qmc_impl_fixed_dim<Real, F, 5>(f, box, n_points, use_scrambling, n_replicates);
        case 6: return integrate_qmc_impl_fixed_dim<Real, F, 6>(f, box, n_points, use_scrambling, n_replicates);
        case 7: return integrate_qmc_impl_fixed_dim<Real, F, 7>(f, box, n_points, use_scrambling, n_replicates);
        case 8: return integrate_qmc_impl_fixed_dim<Real, F, 8>(f, box, n_points, use_scrambling, n_replicates);
        case 9: return integrate_qmc_impl_fixed_dim<Real, F, 9>(f, box, n_points, use_scrambling, n_replicates);
        case 10: return integrate_qmc_impl_fixed_dim<Real, F, 10>(f, box, n_points, use_scrambling, n_replicates);
        case 11: return integrate_qmc_impl_fixed_dim<Real, F, 11>(f, box, n_points, use_scrambling, n_replicates);
        case 12: return integrate_qmc_impl_fixed_dim<Real, F, 12>(f, box, n_points, use_scrambling, n_replicates);
        case 13: return integrate_qmc_impl_fixed_dim<Real, F, 13>(f, box, n_points, use_scrambling, n_replicates);
        case 14: return integrate_qmc_impl_fixed_dim<Real, F, 14>(f, box, n_points, use_scrambling, n_replicates);
        case 15: return integrate_qmc_impl_fixed_dim<Real, F, 15>(f, box, n_points, use_scrambling, n_replicates);
        case 16: return integrate_qmc_impl_fixed_dim<Real, F, 16>(f, box, n_points, use_scrambling, n_replicates);
        case 17: return integrate_qmc_impl_fixed_dim<Real, F, 17>(f, box, n_points, use_scrambling, n_replicates);
        case 18: return integrate_qmc_impl_fixed_dim<Real, F, 18>(f, box, n_points, use_scrambling, n_replicates);
        case 19: return integrate_qmc_impl_fixed_dim<Real, F, 19>(f, box, n_points, use_scrambling, n_replicates);
        case 20: return integrate_qmc_impl_fixed_dim<Real, F, 20>(f, box, n_points, use_scrambling, n_replicates);
        default: {
            result<Real> res;
            res.status = status_code::dimension_error;
            res.error = std::numeric_limits<Real>::max();
            res.value = 0;
            res.evaluations = 0;
            return res;
        }
    }
}

// Helper for vector QMC with fixed dimension
template <typename Real, typename F, std::size_t Dim>
std::vector<result<Real>> integrate_qmc_vector_impl_fixed_dim(
    const F& f,
    const hypercube<Real>& box,
    std::size_t num_components,
    std::size_t n_points,
    bool use_scrambling,
    std::size_t n_replicates)
{
    std::vector<result<Real>> results(num_components);
    
#ifdef BOOST_MATH_HAS_BOOST_RANDOM
    // Storage for replicate values (for RQMC variance estimation)
    std::vector<std::vector<Real>> replicate_values;
    if (use_scrambling) {
        replicate_values.resize(num_components);
        for (auto& rep_vals : replicate_values) {
            rep_vals.reserve(n_replicates);
        }
    }
    
    // Component sums
    std::vector<Real> total_sums(num_components, 0);
    std::size_t total_evals = 0;
    
    // Temporary buffer for function values
    std::vector<Real> values(num_components);
    
    for (std::size_t rep = 0; rep < n_replicates; ++rep) {
        // Create Sobol engine with compile-time dimension
        boost::random::sobol_engine<std::uint64_t, Dim> sobol(0);
        
        // Apply Owen scrambling if requested
        std::uint32_t scramble_seed = 0;
        if (use_scrambling) {
            scramble_seed = static_cast<std::uint32_t>(rep * 65521 + 32749);
        }
        
        // Accumulate QMC sum for this replicate
        std::vector<Real> replicate_sums(num_components, 0);
        std::vector<Real> point(Dim);
        
        for (std::size_t i = 0; i < n_points; ++i) {
            // Generate Sobol point
            auto sobol_point = sobol();
            
            // Convert to [0,1)^d and apply scrambling if needed
            for (std::size_t d = 0; d < Dim; ++d) {
                Real u = Real(sobol_point[d]) / Real(std::numeric_limits<std::uint64_t>::max());
                
                if (use_scrambling) {
                    u = detail::owen_scramble(u, scramble_seed + static_cast<std::uint32_t>(d));
                }
                
                // Apply tent transform for variance reduction
                u = tent_transform(u);
                
                // Map from [0,1] to [a,b]
                point[d] = box.lower[d] + u * (box.upper[d] - box.lower[d]);
            }
            
            // Evaluate vector integrand (single call for all components)
            f(point.data(), values.data(), num_components);
            
            // Accumulate for each component
            for (std::size_t c = 0; c < num_components; ++c) {
                replicate_sums[c] += values[c];
            }
            ++total_evals;
        }
        
        // Scale by volume / n_points
        Real volume = 1.0;
        for (std::size_t d = 0; d < Dim; ++d) {
            volume *= (box.upper[d] - box.lower[d]);
        }
        
        for (std::size_t c = 0; c < num_components; ++c) {
            replicate_sums[c] = (replicate_sums[c] * volume) / Real(n_points);
            
            if (use_scrambling) {
                replicate_values[c].push_back(replicate_sums[c]);
            }
            total_sums[c] += replicate_sums[c];
        }
    }
    
    // Prepare results for each component
    for (std::size_t c = 0; c < num_components; ++c) {
        results[c].evaluations = total_evals;
        results[c].status = status_code::success;
        
        if (use_scrambling && n_replicates > 1) {
            results[c].value = total_sums[c] / Real(n_replicates);
            
            // Compute standard error
            Real mean = results[c].value;
            Real variance = 0;
            for (const auto& val : replicate_values[c]) {
                Real diff = val - mean;
                variance += diff * diff;
            }
            variance /= Real(n_replicates - 1);
            results[c].error = Real(2) * std::sqrt(variance / Real(n_replicates));
        } else {
            results[c].value = total_sums[c] / Real(n_replicates);
            results[c].error = 0;
        }
    }
    
    return results;
#else
    for (auto& res : results) {
        res.status = status_code::not_implemented;
        res.error = std::numeric_limits<Real>::max();
        res.value = 0;
        res.evaluations = 0;
    }
    return results;
#endif
}

/// \brief Extended QMC implementation with configurable sequence types
template <typename Real, typename F, typename Policy>
result<Real> integrate_qmc_impl_ex(
    const F& f,
    const hypercube<Real>& box,
    std::size_t n_points,
    const qmc_options<Real>& options,
    const Policy& /*pol*/)
{
    const std::size_t dim = box.dimension();
    
    // Validate inputs
    if (dim == 0) {
        result<Real> res;
        res.status = status_code::dimension_error;
        res.error = std::numeric_limits<Real>::max();
        res.value = 0;
        res.evaluations = 0;
        return res;
    }
    
    if (n_points == 0) {
        result<Real> res;
        res.status = status_code::invalid_input;
        res.error = std::numeric_limits<Real>::max();
        res.value = 0;
        res.evaluations = 0;
        return res;
    }
    
    // Storage for replicate values (for RQMC variance estimation)
    std::vector<Real> replicate_values;
    if (options.use_scrambling && options.n_replicates > 1) {
        replicate_values.reserve(options.n_replicates);
    }
    
    // Main QMC integration loop
    Real total_sum = 0;
    std::size_t total_evals = 0;
    
    // Compute hypercube volume
    Real volume = 1.0;
    for (std::size_t d = 0; d < dim; ++d) {
        volume *= (box.upper[d] - box.lower[d]);
    }
    
    for (std::size_t rep = 0; rep < options.n_replicates; ++rep) {
        // Create QMC sequence based on options
        std::unique_ptr<qmc_sequence_base<Real>> sequence;
        
        if (options.sequence_type == qmc_sequence_type::lattice_rank1) {
            // Lattice rules need the total number of points
            sequence = create_qmc_sequence<Real>(
                options.sequence_type, dim, options.lattice_points > 0 ? options.lattice_points : n_points);
        } else {
            sequence = create_qmc_sequence<Real>(
                options.sequence_type, dim);
        }
        
        // Special handling for Halton leap
        if (options.sequence_type == qmc_sequence_type::halton && options.halton_leap > 1) {
            auto* halton = dynamic_cast<halton_sequence<Real>*>(sequence.get());
            if (halton) {
                halton->set_leap(options.halton_leap);
            }
        }
        
        // Apply scrambling if requested
        if (options.use_scrambling) {
            std::uint32_t scramble_seed = options.scramble_seed;
            if (scramble_seed == 0) {
                // Generate different seed for each replicate
                scramble_seed = static_cast<std::uint32_t>(rep * 65521 + 32749);
            }
            
            // For lattice rules, apply random shift
            if (options.sequence_type == qmc_sequence_type::lattice_rank1) {
                auto* lattice = dynamic_cast<lattice_rank1_sequence<Real>*>(sequence.get());
                if (lattice) {
                    std::mt19937 rng(scramble_seed);
                    lattice->apply_random_shift(rng);
                }
            } else {
                // For digital sequences, apply Owen scrambling
                sequence = std::make_unique<scrambled_sequence<Real>>(
                    std::move(sequence), scramble_seed);
            }
        }
        
        // Accumulate QMC sum for this replicate
        Real replicate_sum = 0;
        std::vector<Real> point(dim);
        
        // Determine actual number of points for this sequence type
        std::size_t actual_n_points = n_points;
        if (options.sequence_type == qmc_sequence_type::lattice_rank1) {
            actual_n_points = options.lattice_points > 0 ? options.lattice_points : n_points;
        }
        
        for (std::size_t i = 0; i < actual_n_points; ++i) {
            // Check if sequence has more points (relevant for lattice rules)
            if (!sequence->has_next()) {
                break;
            }
            
            // Generate next QMC point
            auto qmc_point = sequence->next();
            
            // Apply tent transform if requested
            for (std::size_t d = 0; d < dim; ++d) {
                Real u = qmc_point[d];
                
                if (options.use_tent_transform) {
                    u = tent_transform(u);
                }
                
                // Map from [0,1] to [a,b]
                point[d] = box.lower[d] + u * (box.upper[d] - box.lower[d]);
            }
            
            // Evaluate integrand
            Real value = evaluate_integrand_impl(f, point);
            replicate_sum += value;
            ++total_evals;
        }
        
        // Scale by volume / actual_n_points
        replicate_sum = (replicate_sum * volume) / Real(actual_n_points);
        
        if (options.use_scrambling && options.n_replicates > 1) {
            replicate_values.push_back(replicate_sum);
        }
        total_sum += replicate_sum;
    }
    
    // Prepare result
    result<Real> res;
    res.evaluations = total_evals;
    res.status = status_code::success;
    
    if (options.use_scrambling && options.n_replicates > 1) {
        res.value = total_sum / Real(options.n_replicates);
        
        // Compute standard error
        Real mean = res.value;
        Real variance = 0;
        for (const auto& val : replicate_values) {
            Real diff = val - mean;
            variance += diff * diff;
        }
        variance /= Real(options.n_replicates - 1);
        res.error = Real(2) * std::sqrt(variance / Real(options.n_replicates));
    } else {
        res.value = total_sum / Real(options.n_replicates);
        res.error = 0;
    }
    
    return res;
}

/// \brief QMC vector implementation with optional randomization
template <typename Real, typename F, typename Policy>
std::vector<result<Real>> integrate_qmc_vector_impl(
    const F& f,
    const hypercube<Real>& box,
    std::size_t num_components,
    std::size_t n_points,
    bool use_scrambling,
    std::size_t n_replicates,
    const Policy& /*pol*/)
{
    std::vector<result<Real>> results(num_components);
    const std::size_t dim = box.dimension();
    
    // Validate inputs
    if (dim == 0 || dim > 20 || num_components == 0) {
        for (auto& res : results) {
            res.status = status_code::dimension_error;
            res.error = std::numeric_limits<Real>::max();
            res.value = 0;
            res.evaluations = 0;
        }
        return results;
    }
    
    if (n_points == 0) {
        for (auto& res : results) {
            res.status = status_code::invalid_input;
            res.error = std::numeric_limits<Real>::max();
            res.value = 0;
            res.evaluations = 0;
        }
        return results;
    }
    
    // Dispatch to compile-time dimension
    switch(dim) {
        case 1: return integrate_qmc_vector_impl_fixed_dim<Real, F, 1>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 2: return integrate_qmc_vector_impl_fixed_dim<Real, F, 2>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 3: return integrate_qmc_vector_impl_fixed_dim<Real, F, 3>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 4: return integrate_qmc_vector_impl_fixed_dim<Real, F, 4>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 5: return integrate_qmc_vector_impl_fixed_dim<Real, F, 5>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 6: return integrate_qmc_vector_impl_fixed_dim<Real, F, 6>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 7: return integrate_qmc_vector_impl_fixed_dim<Real, F, 7>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 8: return integrate_qmc_vector_impl_fixed_dim<Real, F, 8>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 9: return integrate_qmc_vector_impl_fixed_dim<Real, F, 9>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 10: return integrate_qmc_vector_impl_fixed_dim<Real, F, 10>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 11: return integrate_qmc_vector_impl_fixed_dim<Real, F, 11>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 12: return integrate_qmc_vector_impl_fixed_dim<Real, F, 12>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 13: return integrate_qmc_vector_impl_fixed_dim<Real, F, 13>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 14: return integrate_qmc_vector_impl_fixed_dim<Real, F, 14>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 15: return integrate_qmc_vector_impl_fixed_dim<Real, F, 15>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 16: return integrate_qmc_vector_impl_fixed_dim<Real, F, 16>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 17: return integrate_qmc_vector_impl_fixed_dim<Real, F, 17>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 18: return integrate_qmc_vector_impl_fixed_dim<Real, F, 18>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 19: return integrate_qmc_vector_impl_fixed_dim<Real, F, 19>(f, box, num_components, n_points, use_scrambling, n_replicates);
        case 20: return integrate_qmc_vector_impl_fixed_dim<Real, F, 20>(f, box, num_components, n_points, use_scrambling, n_replicates);
        default: {
            for (auto& res : results) {
                res.status = status_code::dimension_error;
                res.error = std::numeric_limits<Real>::max();
                res.value = 0;
                res.evaluations = 0;
            }
            return results;
        }
    }
}

} // namespace detail

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_QMC_HPP
