// Copyright 2025 Shaiko
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_MONTE_CARLO_BASE_HPP
#define BOOST_MATH_CUBATURE_DETAIL_MONTE_CARLO_BASE_HPP

#include <vector>
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>
#include <functional>
#include <array>
#include <type_traits>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Result structure for Monte Carlo integration
/// \tparam Real Floating point type
template <typename Real>
struct mc_result {
    Real value;              ///< Estimated integral value
    Real error;              ///< Standard error estimate
    Real variance;           ///< Sample variance
    std::size_t n_samples;   ///< Number of samples used
    Real volume;             ///< Domain volume
    
    /// \brief Check if result meets tolerance
    bool converged(Real abs_tol, Real rel_tol) const {
        return error < abs_tol || error < std::abs(value) * rel_tol;
    }
    
    /// \brief Confidence interval (default 95%)
    std::pair<Real, Real> confidence_interval(Real z_score = Real(1.96)) const {
        Real margin = z_score * error;
        return {value - margin, value + margin};
    }
};

/// \brief Base class for Monte Carlo sampling strategies
/// \tparam Real Floating point type
template <typename Real>
class monte_carlo_sampler {
public:
    using result_type = mc_result<Real>;
    
    /// \brief Constructor
    /// \param dim Dimension of the integration domain
    /// \param seed Random seed (0 for random seed)
    monte_carlo_sampler(std::size_t dim, std::uint64_t seed = 0)
        : dim_(dim), gen_(seed ? seed : std::random_device{}()),
          uniform_dist_(Real(0), Real(1)) {}
    
    virtual ~monte_carlo_sampler() = default;
    
    /// \brief Sample the function and estimate integral
    /// \param f Function to integrate
    /// \param lower Lower bounds of integration domain
    /// \param upper Upper bounds of integration domain
    /// \param n_samples Number of samples to use
    /// \return Monte Carlo result with value and error estimate
    template <typename F>
    result_type sample(F&& f, const std::vector<Real>& lower,
                      const std::vector<Real>& upper, std::size_t n_samples) {
        if (lower.size() != dim_ || upper.size() != dim_) {
            throw std::invalid_argument("Dimension mismatch in bounds");
        }
        
        // Compute domain volume
        Real volume = Real(1);
        for (std::size_t i = 0; i < dim_; ++i) {
            volume *= (upper[i] - lower[i]);
        }
        
        // Generate samples and evaluate function
        std::vector<Real> samples;
        samples.reserve(n_samples);
        std::vector<Real> point(dim_);
        
        for (std::size_t i = 0; i < n_samples; ++i) {
            generate_point(point, lower, upper, i, n_samples);
            samples.push_back(f(point.data(), dim_));
        }
        
        // Compute statistics
        return compute_statistics(samples, volume);
    }
    
    /// \brief Get the dimension
    std::size_t dimension() const { return dim_; }
    
    /// \brief Reset random number generator
    void reset_generator(std::uint64_t seed = 0) {
        gen_.seed(seed ? seed : std::random_device{}());
    }
    
    /// \brief Generate a sample point (to be overridden by derived classes)
    /// \param point Output point coordinates
    /// \param lower Lower bounds
    /// \param upper Upper bounds
    /// \param sample_idx Current sample index
    /// \param total_samples Total number of samples
    virtual void generate_point(std::vector<Real>& point,
                               const std::vector<Real>& lower,
                               const std::vector<Real>& upper,
                               std::size_t sample_idx,
                               std::size_t total_samples) = 0;
    
protected:
    /// \brief Compute mean, variance, and error from samples
    result_type compute_statistics(const std::vector<Real>& samples, Real volume) const {
        const std::size_t n = samples.size();
        if (n == 0) {
            return {Real(0), std::numeric_limits<Real>::infinity(), Real(0), 0, volume};
        }
        
        // Compute mean using Kahan summation for numerical stability
        Real sum = Real(0);
        Real c = Real(0);  // Compensation for lost low-order bits
        for (const auto& val : samples) {
            Real y = val - c;
            Real t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        Real mean = sum / Real(n);
        
        // Compute variance using two-pass algorithm for stability
        Real var_sum = Real(0);
        Real var_c = Real(0);
        for (const auto& val : samples) {
            Real diff = val - mean;
            Real y = diff * diff - var_c;
            Real t = var_sum + y;
            var_c = (t - var_sum) - y;
            var_sum = t;
        }
        Real variance = var_sum / Real(n - 1);  // Sample variance
        
        // Standard error of the mean
        Real error = std::sqrt(variance / Real(n));
        
        // Scale by volume
        return {mean * volume, error * volume, variance * volume * volume, n, volume};
    }
    
    /// \brief Generate uniform random number in [0,1)
    Real uniform_random() {
        return uniform_dist_(gen_);
    }
    
    std::size_t dim_;
    std::mt19937_64 gen_;
    std::uniform_real_distribution<Real> uniform_dist_;
};

/// \brief Standard Monte Carlo sampler with uniform random sampling
template <typename Real>
class uniform_sampler : public monte_carlo_sampler<Real> {
public:
    using base_type = monte_carlo_sampler<Real>;
    using base_type::monte_carlo_sampler;
    
protected:
    void generate_point(std::vector<Real>& point,
                       const std::vector<Real>& lower,
                       const std::vector<Real>& upper,
                       std::size_t /*sample_idx*/,
                       std::size_t /*total_samples*/) override {
        for (std::size_t i = 0; i < this->dim_; ++i) {
            point[i] = lower[i] + (upper[i] - lower[i]) * this->uniform_random();
        }
    }
};

/// \brief Stratified sampling for variance reduction
template <typename Real>
class stratified_sampler : public monte_carlo_sampler<Real> {
public:
    using base_type = monte_carlo_sampler<Real>;
    
    /// \brief Constructor
    /// \param dim Dimension
    /// \param strata_per_dim Number of strata per dimension
    /// \param seed Random seed
    stratified_sampler(std::size_t dim, std::size_t strata_per_dim = 0,
                      std::uint64_t seed = 0)
        : base_type(dim, seed), strata_per_dim_(strata_per_dim) {
        if (strata_per_dim_ == 0) {
            // Auto-select based on dimension to avoid too many strata
            strata_per_dim_ = (dim <= 3) ? 10 : (dim <= 5) ? 5 : 2;
        }
    }
    
protected:
    void generate_point(std::vector<Real>& point,
                       const std::vector<Real>& lower,
                       const std::vector<Real>& upper,
                       std::size_t sample_idx,
                       std::size_t total_samples) override {
        // Calculate total number of strata
        std::size_t total_strata = 1;
        for (std::size_t i = 0; i < this->dim_; ++i) {
            total_strata *= strata_per_dim_;
        }
        
        // Determine samples per stratum
        std::size_t samples_per_stratum = std::max(std::size_t(1), 
                                                   total_samples / total_strata);
        
        // Find which stratum this sample belongs to
        std::size_t stratum_idx = sample_idx / samples_per_stratum;
        if (stratum_idx >= total_strata) {
            stratum_idx = total_strata - 1;
        }
        
        // Convert stratum index to multi-dimensional indices
        std::vector<std::size_t> stratum_coords(this->dim_);
        std::size_t temp = stratum_idx;
        for (std::size_t i = 0; i < this->dim_; ++i) {
            stratum_coords[i] = temp % strata_per_dim_;
            temp /= strata_per_dim_;
        }
        
        // Generate point within the stratum with jittering
        for (std::size_t i = 0; i < this->dim_; ++i) {
            Real width = (upper[i] - lower[i]) / Real(strata_per_dim_);
            Real stratum_lower = lower[i] + stratum_coords[i] * width;
            Real stratum_upper = stratum_lower + width;
            
            // Random point within stratum
            point[i] = stratum_lower + (stratum_upper - stratum_lower) * this->uniform_random();
        }
    }
    
private:
    std::size_t strata_per_dim_;
};

/// \brief Latin Hypercube Sampling (LHS) for better space-filling
template <typename Real>
class latin_hypercube_sampler : public monte_carlo_sampler<Real> {
public:
    using base_type = monte_carlo_sampler<Real>;
    
    latin_hypercube_sampler(std::size_t dim, std::uint64_t seed = 0)
        : base_type(dim, seed), current_batch_size_(0), batch_index_(0) {}
    
protected:
    void generate_point(std::vector<Real>& point,
                       const std::vector<Real>& lower,
                       const std::vector<Real>& upper,
                       std::size_t sample_idx,
                       std::size_t total_samples) override {
        // Generate new batch if needed
        if (sample_idx == 0 || sample_idx >= current_batch_size_) {
            generate_lhs_batch(total_samples, lower, upper);
            batch_index_ = 0;
        }
        
        // Return point from batch
        for (std::size_t i = 0; i < this->dim_; ++i) {
            point[i] = lhs_points_[batch_index_][i];
        }
        ++batch_index_;
    }
    
private:
    void generate_lhs_batch(std::size_t n_samples,
                           const std::vector<Real>& lower,
                           const std::vector<Real>& upper) {
        current_batch_size_ = n_samples;
        lhs_points_.clear();
        lhs_points_.reserve(n_samples);
        
        // Create permutations for each dimension
        std::vector<std::vector<std::size_t>> permutations(this->dim_);
        for (std::size_t d = 0; d < this->dim_; ++d) {
            permutations[d].resize(n_samples);
            std::iota(permutations[d].begin(), permutations[d].end(), 0);
            std::shuffle(permutations[d].begin(), permutations[d].end(), this->gen_);
        }
        
        // Generate LHS points
        for (std::size_t i = 0; i < n_samples; ++i) {
            std::vector<Real> point(this->dim_);
            for (std::size_t d = 0; d < this->dim_; ++d) {
                Real cell_width = (upper[d] - lower[d]) / Real(n_samples);
                Real jitter = this->uniform_random();
                point[d] = lower[d] + cell_width * (Real(permutations[d][i]) + jitter);
            }
            lhs_points_.push_back(std::move(point));
        }
    }
    
    std::vector<std::vector<Real>> lhs_points_;
    std::size_t current_batch_size_;
    std::size_t batch_index_;
};

/// \brief Utility functions for error estimation
namespace mc_utils {

/// \brief Compute effective sample size for weighted samples
template <typename Real>
Real effective_sample_size(const std::vector<Real>& weights) {
    Real sum_w = std::accumulate(weights.begin(), weights.end(), Real(0));
    Real sum_w2 = std::inner_product(weights.begin(), weights.end(),
                                     weights.begin(), Real(0));
    return (sum_w * sum_w) / sum_w2;
}

/// \brief Batch means method for error estimation
template <typename Real>
Real batch_means_error(const std::vector<Real>& samples, std::size_t n_batches = 10) {
    if (samples.size() < n_batches * 2) {
        return std::numeric_limits<Real>::infinity();
    }
    
    std::size_t batch_size = samples.size() / n_batches;
    std::vector<Real> batch_means;
    batch_means.reserve(n_batches);
    
    for (std::size_t b = 0; b < n_batches; ++b) {
        Real sum = Real(0);
        for (std::size_t i = b * batch_size; i < (b + 1) * batch_size; ++i) {
            sum += samples[i];
        }
        batch_means.push_back(sum / Real(batch_size));
    }
    
    // Compute variance of batch means
    Real mean = std::accumulate(batch_means.begin(), batch_means.end(), Real(0)) 
                / Real(n_batches);
    Real var = Real(0);
    for (const auto& bm : batch_means) {
        Real diff = bm - mean;
        var += diff * diff;
    }
    var /= Real(n_batches - 1);
    
    return std::sqrt(var / Real(n_batches));
}

/// \brief Jackknife error estimation
template <typename Real>
Real jackknife_error(const std::vector<Real>& samples) {
    const std::size_t n = samples.size();
    if (n < 2) return std::numeric_limits<Real>::infinity();
    
    Real total_sum = std::accumulate(samples.begin(), samples.end(), Real(0));
    Real full_mean = total_sum / Real(n);
    
    std::vector<Real> jackknife_estimates;
    jackknife_estimates.reserve(n);
    
    for (std::size_t i = 0; i < n; ++i) {
        Real jack_sum = total_sum - samples[i];
        jackknife_estimates.push_back(jack_sum / Real(n - 1));
    }
    
    // Compute variance of jackknife estimates
    Real var = Real(0);
    for (const auto& est : jackknife_estimates) {
        Real diff = est - full_mean;
        var += diff * diff;
    }
    var *= Real(n - 1) / Real(n);
    
    return std::sqrt(var);
}

} // namespace mc_utils

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_MONTE_CARLO_BASE_HPP