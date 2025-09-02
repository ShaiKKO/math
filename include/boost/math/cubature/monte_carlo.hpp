// Copyright 2025 Boost.Math Contributors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_MONTE_CARLO_HPP
#define BOOST_MATH_CUBATURE_MONTE_CARLO_HPP

/// \file monte_carlo.hpp
/// \brief Monte Carlo integration with variance reduction techniques
///
/// This header provides Monte Carlo integration methods for high-dimensional
/// integration problems. It includes:
/// - Standard Monte Carlo with uniform sampling
/// - Stratified sampling for variance reduction
/// - Antithetic variates
/// - Latin hypercube sampling
/// - Importance sampling
/// - Adaptive sampling with automatic convergence checking
///
/// Example usage:
/// \code
/// #include <boost/math/cubature/monte_carlo.hpp>
/// 
/// // Integrate exp(-x^2-y^2) over [-1,1]^2
/// auto f = [](const double* x, std::size_t) { 
///     return std::exp(-x[0]*x[0] - x[1]*x[1]); 
/// };
/// 
/// std::vector<double> lower = {-1, -1};
/// std::vector<double> upper = {1, 1};
/// 
/// auto result = boost::math::cubature::integrate_monte_carlo(
///     f, lower, upper, 1e-3, 1e-3);
/// 
/// std::cout << "Integral: " << result.value 
///           << " +/- " << result.error << std::endl;
/// \endcode

#include <boost/math/cubature/detail/monte_carlo_base.hpp>
#include <boost/math/cubature/detail/variance_reduction.hpp>
#include <boost/math/cubature/detail/adaptive_monte_carlo.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <type_traits>
#include <limits>

namespace boost { namespace math { namespace cubature {

/// \brief Monte Carlo integration options
template <typename Real = double>
struct monte_carlo_options {
    /// Sampling strategy to use
    enum class strategy {
        automatic,         ///< Automatically select based on dimension
        uniform,          ///< Standard Monte Carlo
        stratified,       ///< Stratified sampling
        antithetic,       ///< Antithetic variates
        latin_hypercube,  ///< Latin hypercube sampling
        importance,       ///< Importance sampling
        randomized_qmc    ///< Randomized quasi-Monte Carlo
    };
    
    strategy method = strategy::automatic;
    std::size_t min_samples = 0;        ///< Minimum samples (0 for auto)
    std::size_t max_samples = 0;        ///< Maximum samples (0 for auto)
    std::size_t batch_size = 0;         ///< Batch size (0 for auto)
    Real confidence_level = Real(0.95); ///< Confidence level for error estimates
    std::uint64_t seed = 0;             ///< Random seed (0 for random)
    std::size_t max_time_ms = 0;        ///< Time limit in milliseconds (0 for none)
    bool adaptive = true;               ///< Use adaptive sampling
    
    /// \brief Get default options for given dimension
    static monte_carlo_options defaults(std::size_t dim) {
        monte_carlo_options opts;
        
        // Auto-select strategy based on dimension
        if (dim <= 3) {
            opts.method = strategy::stratified;
            opts.min_samples = 1000;
            opts.max_samples = 1000000;
            opts.batch_size = 1000;
        } else if (dim <= 10) {
            opts.method = strategy::latin_hypercube;
            opts.min_samples = 5000;
            opts.max_samples = 5000000;
            opts.batch_size = 5000;
        } else if (dim <= 20) {
            opts.method = strategy::antithetic;
            opts.min_samples = 10000;
            opts.max_samples = 10000000;
            opts.batch_size = 10000;
        } else {
            // High dimensions: standard MC often best
            opts.method = strategy::uniform;
            opts.min_samples = 50000;
            opts.max_samples = 50000000;
            opts.batch_size = 50000;
        }
        
        return opts;
    }
};

/// \brief Result type for Monte Carlo integration
template <typename Real>
using mc_result = detail::mc_result<Real>;

/// \brief Main Monte Carlo integration function
/// \param f Function to integrate (callable with signature Real(const Real*, std::size_t))
/// \param lower Lower bounds of integration domain
/// \param upper Upper bounds of integration domain
/// \param abs_tol Absolute tolerance for convergence
/// \param rel_tol Relative tolerance for convergence
/// \param options Monte Carlo options (optional)
/// \param pol Boost.Math policy (optional)
/// \return Monte Carlo result with value, error estimate, and statistics
template <typename F, typename Real = double, 
          typename Policy = policies::policy<>>
mc_result<Real> integrate_monte_carlo(
    F&& f,
    const std::vector<Real>& lower,
    const std::vector<Real>& upper,
    Real abs_tol,
    Real rel_tol,
    const monte_carlo_options<Real>& options = monte_carlo_options<Real>(),
    const Policy& pol = policies::policy<>()) {
    
    // Input validation
    if (lower.size() != upper.size()) {
        policies::raise_domain_error(
            "integrate_monte_carlo",
            "Lower and upper bounds must have same dimension, got %1%",
            static_cast<Real>(lower.size()), pol);
        return mc_result<Real>{Real(0), std::numeric_limits<Real>::infinity(), 
                               Real(0), 0, Real(0)};
    }
    
    const std::size_t dim = lower.size();
    if (dim == 0) {
        policies::raise_domain_error(
            "integrate_monte_carlo",
            "Dimension must be positive, got %1%",
            static_cast<Real>(0), pol);
        return mc_result<Real>{Real(0), std::numeric_limits<Real>::infinity(), 
                               Real(0), 0, Real(0)};
    }
    
    // Check bounds validity
    for (std::size_t i = 0; i < dim; ++i) {
        if (lower[i] >= upper[i]) {
            policies::raise_domain_error(
                "integrate_monte_carlo",
                "Invalid bounds: lower[%1%] >= upper[%1%]",
                static_cast<Real>(i), pol);
            return mc_result<Real>{Real(0), std::numeric_limits<Real>::infinity(), 
                                   Real(0), 0, Real(0)};
        }
    }
    
    // Get options with defaults if needed
    monte_carlo_options<Real> opts = options;
    if (opts.method == monte_carlo_options<Real>::strategy::automatic) {
        opts = monte_carlo_options<Real>::defaults(dim);
        // Preserve user-specified values if set
        if (options.min_samples > 0) opts.min_samples = options.min_samples;
        if (options.max_samples > 0) opts.max_samples = options.max_samples;
        if (options.batch_size > 0) opts.batch_size = options.batch_size;
        opts.confidence_level = options.confidence_level;
        opts.seed = options.seed;
        opts.max_time_ms = options.max_time_ms;
        opts.adaptive = options.adaptive;
    }
    
    // Set defaults for unspecified parameters
    if (opts.min_samples == 0) {
        opts.min_samples = std::max(std::size_t(100), std::size_t(10 * dim));
    }
    if (opts.max_samples == 0) {
        opts.max_samples = std::max(std::size_t(1000000), std::size_t(100000 * dim));
    }
    if (opts.batch_size == 0) {
        opts.batch_size = std::max(std::size_t(1000), opts.min_samples / 10);
    }
    
    // Choose integration method
    if (opts.adaptive) {
        // Use adaptive Monte Carlo
        using strategy_t = typename detail::adaptive_monte_carlo<Real>::strategy;
        
        strategy_t strat;
        switch (opts.method) {
            case monte_carlo_options<Real>::strategy::uniform:
                strat = strategy_t::uniform;
                break;
            case monte_carlo_options<Real>::strategy::stratified:
                strat = strategy_t::stratified;
                break;
            case monte_carlo_options<Real>::strategy::antithetic:
                strat = strategy_t::antithetic;
                break;
            case monte_carlo_options<Real>::strategy::latin_hypercube:
                strat = strategy_t::latin_hypercube;
                break;
            case monte_carlo_options<Real>::strategy::importance:
                strat = strategy_t::importance;
                break;
            case monte_carlo_options<Real>::strategy::randomized_qmc:
                strat = strategy_t::randomized_qmc;
                break;
            default:
                strat = strategy_t::uniform;
        }
        
        detail::adaptive_monte_carlo<Real> mc(dim, strat, opts.seed);
        mc.set_parameters(opts.min_samples, opts.max_samples, opts.batch_size);
        mc.set_confidence_level(opts.confidence_level);
        
        return mc.integrate(std::forward<F>(f), lower, upper, 
                           abs_tol, rel_tol, opts.max_time_ms);
    } else {
        // Use fixed-sample Monte Carlo
        std::unique_ptr<detail::monte_carlo_sampler<Real>> sampler;
        
        switch (opts.method) {
            case monte_carlo_options<Real>::strategy::stratified:
                sampler = std::make_unique<detail::stratified_sampler<Real>>(dim, 0, opts.seed);
                break;
            case monte_carlo_options<Real>::strategy::antithetic:
                sampler = std::make_unique<detail::antithetic_sampler<Real>>(dim, opts.seed);
                break;
            case monte_carlo_options<Real>::strategy::latin_hypercube:
                sampler = std::make_unique<detail::latin_hypercube_sampler<Real>>(dim, opts.seed);
                break;
            case monte_carlo_options<Real>::strategy::importance:
                sampler = std::make_unique<detail::importance_sampler<Real>>(dim, opts.seed);
                break;
            case monte_carlo_options<Real>::strategy::randomized_qmc:
                sampler = std::make_unique<detail::randomized_qmc_sampler<Real>>(dim, 10, opts.seed);
                break;
            default:
                sampler = std::make_unique<detail::uniform_sampler<Real>>(dim, opts.seed);
        }
        
        return sampler->sample(std::forward<F>(f), lower, upper, opts.min_samples);
    }
}

/// \brief Simplified Monte Carlo integration with automatic settings
/// \param f Function to integrate
/// \param lower Lower bounds
/// \param upper Upper bounds
/// \param n_samples Number of samples to use
/// \param seed Random seed (0 for random)
template <typename F, typename Real = double>
mc_result<Real> monte_carlo(F&& f,
                           const std::vector<Real>& lower,
                           const std::vector<Real>& upper,
                           std::size_t n_samples,
                           std::uint64_t seed = 0) {
    monte_carlo_options<Real> opts;
    opts.min_samples = n_samples;
    opts.max_samples = n_samples;
    opts.adaptive = false;
    opts.seed = seed;
    
    // Use high tolerances since we're doing fixed sampling
    Real abs_tol = std::numeric_limits<Real>::infinity();
    Real rel_tol = std::numeric_limits<Real>::infinity();
    
    return integrate_monte_carlo(std::forward<F>(f), lower, upper, 
                                abs_tol, rel_tol, opts);
}

/// \brief Monte Carlo with control variates
/// \param f Target function
/// \param g Control function with known integral
/// \param g_integral Known integral of g
/// \param lower Lower bounds
/// \param upper Upper bounds
/// \param n_samples Number of samples
/// \param seed Random seed
template <typename F, typename G, typename Real = double>
mc_result<Real> monte_carlo_control_variate(
    F&& f, G&& g, Real g_integral,
    const std::vector<Real>& lower,
    const std::vector<Real>& upper,
    std::size_t n_samples,
    std::uint64_t seed = 0) {
    
    const std::size_t dim = lower.size();
    detail::control_variate_estimator<Real> cv(dim, seed);
    
    return cv.estimate(std::forward<F>(f), std::forward<G>(g), g_integral,
                      lower, upper, n_samples, true);
}

/// \brief Estimate sample size needed for target accuracy
/// \param dim Dimension of problem
/// \param variance_estimate Estimated variance (or 0 for unknown)
/// \param target_error Target error
/// \param confidence_level Confidence level (default 0.95)
/// \return Estimated number of samples needed
template <typename Real = double>
std::size_t estimate_sample_size(std::size_t dim,
                                Real variance_estimate,
                                Real target_error,
                                Real confidence_level = Real(0.95)) {
    // If variance unknown, use pessimistic estimate
    if (variance_estimate <= Real(0)) {
        // Assume coefficient of variation of 1
        variance_estimate = Real(1);
    }
    
    // Z-score for confidence level
    Real z_score = Real(1.96);  // 95% confidence
    if (confidence_level <= Real(0.9)) z_score = Real(1.645);
    else if (confidence_level >= Real(0.99)) z_score = Real(2.576);
    
    // n = (z * sigma / error)^2
    Real n_real = std::pow(z_score * std::sqrt(variance_estimate) / target_error, 2);
    
    // Add safety factor for high dimensions
    if (dim > 10) {
        n_real *= Real(1) + Real(dim) / Real(100);
    }
    
    return static_cast<std::size_t>(std::ceil(n_real));
}

/// \brief Helper to create importance sampler with Gaussian proposal
/// \param f Function to integrate
/// \param lower Lower bounds
/// \param upper Upper bounds
/// \param center Center of Gaussian proposal (in [0,1]^d)
/// \param std_dev Standard deviation
/// \param n_samples Number of samples
/// \param seed Random seed
template <typename F, typename Real = double>
mc_result<Real> monte_carlo_importance_gaussian(
    F&& f,
    const std::vector<Real>& lower,
    const std::vector<Real>& upper,
    const std::vector<Real>& center,
    Real std_dev,
    std::size_t n_samples,
    std::uint64_t seed = 0) {
    
    const std::size_t dim = lower.size();
    if (center.size() != dim) {
        throw std::invalid_argument("Center dimension mismatch");
    }
    
    detail::importance_sampler<Real> sampler(dim, seed);
    
    // Transform center to unit cube
    std::vector<Real> unit_center(dim);
    for (std::size_t i = 0; i < dim; ++i) {
        unit_center[i] = (center[i] - lower[i]) / (upper[i] - lower[i]);
    }
    
    sampler.set_gaussian_proposal(unit_center, std_dev);
    
    return sampler.sample(std::forward<F>(f), lower, upper, n_samples);
}

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_MONTE_CARLO_HPP