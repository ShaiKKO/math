// Copyright 2025 Boost.Math Contributors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_VARIANCE_REDUCTION_HPP
#define BOOST_MATH_CUBATURE_DETAIL_VARIANCE_REDUCTION_HPP

#include <boost/math/cubature/detail/monte_carlo_base.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Antithetic variates sampler for variance reduction
/// \details Uses negatively correlated sample pairs to reduce variance
/// \tparam Real Floating point type
template <typename Real>
class antithetic_sampler : public monte_carlo_sampler<Real> {
public:
    using base_type = monte_carlo_sampler<Real>;
    using result_type = typename base_type::result_type;
    
    /// \brief Constructor
    /// \param dim Dimension of integration domain
    /// \param seed Random seed (0 for random seed)
    antithetic_sampler(std::size_t dim, std::uint64_t seed = 0)
        : base_type(dim, seed), use_antithetic_(false), last_u_(dim) {}
    
    /// \brief Sample with antithetic variates
    template <typename F>
    result_type sample(F&& f, const std::vector<Real>& lower,
                      const std::vector<Real>& upper, std::size_t n_samples) {
        if (lower.size() != this->dim_ || upper.size() != this->dim_) {
            throw std::invalid_argument("Dimension mismatch in bounds");
        }
        
        // Ensure even number of samples for pairing
        n_samples = (n_samples + 1) & ~std::size_t(1);
        
        // Compute domain volume
        Real volume = Real(1);
        std::vector<Real> widths(this->dim_);
        for (std::size_t i = 0; i < this->dim_; ++i) {
            widths[i] = upper[i] - lower[i];
            volume *= widths[i];
        }
        
        // Generate antithetic pairs and evaluate
        std::vector<Real> samples;
        samples.reserve(n_samples);
        std::vector<Real> point(this->dim_);
        std::vector<Real> anti_point(this->dim_);
        
        for (std::size_t i = 0; i < n_samples; i += 2) {
            // Generate primary point
            for (std::size_t d = 0; d < this->dim_; ++d) {
                Real u = this->uniform_random();
                point[d] = lower[d] + widths[d] * u;
                // Antithetic point uses 1-u
                anti_point[d] = lower[d] + widths[d] * (Real(1) - u);
            }
            
            // Evaluate both points
            Real f1 = f(point.data(), this->dim_);
            Real f2 = f(anti_point.data(), this->dim_);
            
            // Average of antithetic pair reduces variance
            samples.push_back((f1 + f2) / Real(2));
        }
        
        // Compute statistics on paired samples
        return compute_antithetic_statistics(samples, volume, 2);
    }
    
protected:
    void generate_point(std::vector<Real>& point,
                       const std::vector<Real>& lower,
                       const std::vector<Real>& upper,
                       std::size_t sample_idx,
                       std::size_t /*total_samples*/) override {
        // Toggle between regular and antithetic points
        if (sample_idx % 2 == 0) {
            // Generate new random point
            for (std::size_t i = 0; i < this->dim_; ++i) {
                last_u_[i] = this->uniform_random();
                point[i] = lower[i] + (upper[i] - lower[i]) * last_u_[i];
            }
        } else {
            // Generate antithetic point
            for (std::size_t i = 0; i < this->dim_; ++i) {
                point[i] = lower[i] + (upper[i] - lower[i]) * (Real(1) - last_u_[i]);
            }
        }
    }
    
private:
    /// \brief Compute statistics accounting for antithetic pairing
    result_type compute_antithetic_statistics(const std::vector<Real>& paired_samples, 
                                             Real volume, std::size_t pair_size) const {
        const std::size_t n = paired_samples.size();
        if (n == 0) {
            return {Real(0), std::numeric_limits<Real>::infinity(), Real(0), 0, volume};
        }
        
        // Compute mean with Kahan summation
        Real sum = Real(0);
        Real c = Real(0);
        for (const auto& val : paired_samples) {
            Real y = val - c;
            Real t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        Real mean = sum / Real(n);
        
        // Compute variance of paired samples
        Real var_sum = Real(0);
        Real var_c = Real(0);
        for (const auto& val : paired_samples) {
            Real diff = val - mean;
            Real y = diff * diff - var_c;
            Real t = var_sum + y;
            var_c = (t - var_sum) - y;
            var_sum = t;
        }
        
        // Adjust variance for paired sampling
        Real variance = var_sum / Real(n - 1);
        Real error = std::sqrt(variance / Real(n));
        
        // Scale by volume and account for actual sample count
        return {mean * volume, error * volume, variance * volume * volume, 
                n * pair_size, volume};
    }
    
    bool use_antithetic_;
    std::vector<Real> last_u_;  // Store last uniform samples for antithetic generation
};

/// \brief Control variate estimator for variance reduction
/// \details Uses a control function with known integral to reduce variance
/// \tparam Real Floating point type
template <typename Real>
class control_variate_estimator {
public:
    using result_type = mc_result<Real>;
    
    /// \brief Constructor
    /// \param dim Dimension
    /// \param seed Random seed
    control_variate_estimator(std::size_t dim, std::uint64_t seed = 0)
        : dim_(dim), gen_(seed ? seed : std::random_device{}()),
          uniform_dist_(Real(0), Real(1)), beta_(Real(1)) {}
    
    /// \brief Estimate integral using control variates
    /// \param f Target function
    /// \param g Control function with known integral
    /// \param g_integral Known integral of g
    /// \param lower Lower bounds
    /// \param upper Upper bounds
    /// \param n_samples Number of samples
    /// \param estimate_beta Whether to estimate optimal beta
    template <typename F, typename G>
    result_type estimate(F&& f, G&& g, Real g_integral,
                        const std::vector<Real>& lower,
                        const std::vector<Real>& upper,
                        std::size_t n_samples,
                        bool estimate_beta = true) {
        if (lower.size() != dim_ || upper.size() != dim_) {
            throw std::invalid_argument("Dimension mismatch in bounds");
        }
        
        // Compute volume
        Real volume = Real(1);
        for (std::size_t i = 0; i < dim_; ++i) {
            volume *= (upper[i] - lower[i]);
        }
        
        // Generate samples
        std::vector<Real> f_samples, g_samples;
        f_samples.reserve(n_samples);
        g_samples.reserve(n_samples);
        
        std::vector<Real> point(dim_);
        for (std::size_t i = 0; i < n_samples; ++i) {
            // Generate uniform random point
            for (std::size_t d = 0; d < dim_; ++d) {
                point[d] = lower[d] + (upper[d] - lower[d]) * uniform_dist_(gen_);
            }
            
            f_samples.push_back(f(point.data(), dim_));
            g_samples.push_back(g(point.data(), dim_));
        }
        
        // Estimate optimal beta if requested
        if (estimate_beta && n_samples > 10) {
            beta_ = estimate_optimal_beta(f_samples, g_samples);
        }
        
        // Apply control variate correction
        std::vector<Real> corrected_samples(n_samples);
        
        for (std::size_t i = 0; i < n_samples; ++i) {
            corrected_samples[i] = f_samples[i] - beta_ * (g_samples[i] - g_integral/volume);
        }
        
        // Compute statistics on corrected samples
        return compute_statistics(corrected_samples, volume);
    }
    
    /// \brief Set the control coefficient manually
    void set_beta(Real beta) { beta_ = beta; }
    
    /// \brief Get current control coefficient
    Real get_beta() const { return beta_; }
    
private:
    /// \brief Estimate optimal beta coefficient
    Real estimate_optimal_beta(const std::vector<Real>& f_samples,
                               const std::vector<Real>& g_samples) {
        const std::size_t n = f_samples.size();
        
        // Compute means
        Real f_mean = std::accumulate(f_samples.begin(), f_samples.end(), Real(0)) / Real(n);
        Real g_mean = std::accumulate(g_samples.begin(), g_samples.end(), Real(0)) / Real(n);
        
        // Compute covariance and variance
        Real cov = Real(0);
        Real var_g = Real(0);
        
        for (std::size_t i = 0; i < n; ++i) {
            Real f_diff = f_samples[i] - f_mean;
            Real g_diff = g_samples[i] - g_mean;
            cov += f_diff * g_diff;
            var_g += g_diff * g_diff;
        }
        
        // Avoid division by zero
        if (std::abs(var_g) < std::numeric_limits<Real>::epsilon() * Real(n)) {
            return Real(0);
        }
        
        return cov / var_g;
    }
    
    /// \brief Compute statistics from samples
    result_type compute_statistics(const std::vector<Real>& samples, Real volume) const {
        const std::size_t n = samples.size();
        if (n == 0) {
            return {Real(0), std::numeric_limits<Real>::infinity(), Real(0), 0, volume};
        }
        
        // Compute mean using Kahan summation
        Real sum = Real(0);
        Real c = Real(0);
        for (const auto& val : samples) {
            Real y = val - c;
            Real t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        Real mean = sum / Real(n);
        
        // Compute variance
        Real var_sum = Real(0);
        Real var_c = Real(0);
        for (const auto& val : samples) {
            Real diff = val - mean;
            Real y = diff * diff - var_c;
            Real t = var_sum + y;
            var_c = (t - var_sum) - y;
            var_sum = t;
        }
        Real variance = var_sum / Real(n - 1);
        Real error = std::sqrt(variance / Real(n));
        
        return {mean * volume, error * volume, variance * volume * volume, n, volume};
    }
    
    std::size_t dim_;
    std::mt19937_64 gen_;
    std::uniform_real_distribution<Real> uniform_dist_;
    Real beta_;  // Control coefficient
};

/// \brief Importance sampling for variance reduction
/// \details Sample from a proposal distribution and reweight
/// \tparam Real Floating point type
template <typename Real>
class importance_sampler : public monte_carlo_sampler<Real> {
public:
    using base_type = monte_carlo_sampler<Real>;
    using result_type = typename base_type::result_type;
    using proposal_fn = std::function<Real(const Real*, std::size_t)>;
    using sample_fn = std::function<void(Real*, std::size_t)>;
    
    /// \brief Constructor with custom proposal
    /// \param dim Dimension
    /// \param proposal_pdf Proposal probability density function
    /// \param proposal_sampler Function to generate samples from proposal
    /// \param seed Random seed
    importance_sampler(std::size_t dim, 
                      proposal_fn proposal_pdf,
                      sample_fn proposal_sampler,
                      std::uint64_t seed = 0)
        : base_type(dim, seed), 
          proposal_pdf_(proposal_pdf),
          proposal_sampler_(proposal_sampler),
          use_default_proposal_(false) {}
    
    /// \brief Constructor with default proposal (uniform or Gaussian)
    /// \param dim Dimension
    /// \param seed Random seed
    importance_sampler(std::size_t dim, std::uint64_t seed = 0)
        : base_type(dim, seed), use_default_proposal_(true) {
        // Default to uniform proposal
        proposal_pdf_ = [](const Real*, std::size_t) { return Real(1); };
        proposal_sampler_ = [this](Real* point, std::size_t d) {
            for (std::size_t i = 0; i < d; ++i) {
                point[i] = this->uniform_random();
            }
        };
    }
    
    /// \brief Sample with importance weighting
    template <typename F>
    result_type sample(F&& f, const std::vector<Real>& lower,
                      const std::vector<Real>& upper, std::size_t n_samples) {
        if (lower.size() != this->dim_ || upper.size() != this->dim_) {
            throw std::invalid_argument("Dimension mismatch in bounds");
        }
        
        // Compute volume
        Real volume = Real(1);
        for (std::size_t i = 0; i < this->dim_; ++i) {
            volume *= (upper[i] - lower[i]);
        }
        
        std::vector<Real> weighted_samples;
        std::vector<Real> weights;
        weighted_samples.reserve(n_samples);
        weights.reserve(n_samples);
        
        std::vector<Real> point(this->dim_);
        std::vector<Real> unit_point(this->dim_);
        
        for (std::size_t i = 0; i < n_samples; ++i) {
            // Sample from proposal distribution
            if (use_default_proposal_) {
                // Default uniform sampling
                for (std::size_t d = 0; d < this->dim_; ++d) {
                    unit_point[d] = this->uniform_random();
                    point[d] = lower[d] + (upper[d] - lower[d]) * unit_point[d];
                }
            } else {
                // Custom proposal sampler
                proposal_sampler_(unit_point.data(), this->dim_);
                // Transform to domain
                for (std::size_t d = 0; d < this->dim_; ++d) {
                    point[d] = lower[d] + (upper[d] - lower[d]) * unit_point[d];
                }
            }
            
            // Compute importance weight
            Real proposal_density = proposal_pdf_(unit_point.data(), this->dim_);
            if (proposal_density <= Real(0)) {
                continue;  // Skip zero-probability samples
            }
            
            Real weight = Real(1) / proposal_density;  // Uniform target / proposal
            Real f_val = f(point.data(), this->dim_);
            
            weighted_samples.push_back(f_val * weight);
            weights.push_back(weight);
        }
        
        // Self-normalized importance sampling
        return compute_weighted_statistics(weighted_samples, weights, volume);
    }
    
    /// \brief Set Gaussian proposal centered at a point
    /// \param center Center of Gaussian proposal
    /// \param std_dev Standard deviation for each dimension
    void set_gaussian_proposal(const std::vector<Real>& center, Real std_dev) {
        if (center.size() != this->dim_) {
            throw std::invalid_argument("Center dimension mismatch");
        }
        
        center_ = center;
        std_dev_ = std_dev;
        use_default_proposal_ = false;
        
        // Normal distribution for sampling
        normal_dist_ = std::normal_distribution<Real>(Real(0), std_dev);
        
        // Gaussian proposal PDF
        proposal_pdf_ = [this](const Real* x, std::size_t d) {
            Real log_prob = Real(0);
            const Real norm_const = std::pow(Real(2) * M_PI * std_dev_ * std_dev_, 
                                           -Real(d) / Real(2));
            
            for (std::size_t i = 0; i < d; ++i) {
                Real diff = x[i] - center_[i];
                log_prob -= diff * diff / (Real(2) * std_dev_ * std_dev_);
            }
            
            return norm_const * std::exp(log_prob);
        };
        
        // Gaussian proposal sampler
        proposal_sampler_ = [this](Real* point, std::size_t d) {
            for (std::size_t i = 0; i < d; ++i) {
                point[i] = center_[i] + normal_dist_(this->gen_);
                // Clip to [0,1]
                point[i] = std::max(Real(0), std::min(Real(1), point[i]));
            }
        };
    }
    
protected:
    void generate_point(std::vector<Real>& point,
                       const std::vector<Real>& lower,
                       const std::vector<Real>& upper,
                       std::size_t /*sample_idx*/,
                       std::size_t /*total_samples*/) override {
        std::vector<Real> unit_point(this->dim_);
        
        if (use_default_proposal_) {
            for (std::size_t i = 0; i < this->dim_; ++i) {
                unit_point[i] = this->uniform_random();
                point[i] = lower[i] + (upper[i] - lower[i]) * unit_point[i];
            }
        } else {
            proposal_sampler_(unit_point.data(), this->dim_);
            for (std::size_t i = 0; i < this->dim_; ++i) {
                point[i] = lower[i] + (upper[i] - lower[i]) * unit_point[i];
            }
        }
    }
    
private:
    /// \brief Compute statistics for weighted samples
    result_type compute_weighted_statistics(const std::vector<Real>& weighted_samples,
                                           const std::vector<Real>& weights,
                                           Real volume) const {
        const std::size_t n = weighted_samples.size();
        if (n == 0) {
            return {Real(0), std::numeric_limits<Real>::infinity(), Real(0), 0, volume};
        }
        
        // Compute weighted mean
        Real sum_wx = Real(0);
        Real sum_w = Real(0);
        
        for (std::size_t i = 0; i < n; ++i) {
            sum_wx += weighted_samples[i];
            sum_w += weights[i];
        }
        
        if (sum_w <= Real(0)) {
            return {Real(0), std::numeric_limits<Real>::infinity(), Real(0), n, volume};
        }
        
        Real mean = sum_wx / sum_w;
        
        // Compute effective sample size
        Real sum_w2 = std::inner_product(weights.begin(), weights.end(),
                                        weights.begin(), Real(0));
        Real n_eff = (sum_w * sum_w) / sum_w2;
        
        // Compute weighted variance
        Real var_sum = Real(0);
        for (std::size_t i = 0; i < n; ++i) {
            Real diff = weighted_samples[i] / weights[i] - mean;
            var_sum += weights[i] * diff * diff;
        }
        Real variance = var_sum / sum_w;
        
        // Standard error accounting for effective sample size
        Real error = std::sqrt(variance / n_eff);
        
        // Scale by volume
        return {mean * volume, error * volume, variance * volume * volume, n, volume};
    }
    
    proposal_fn proposal_pdf_;
    sample_fn proposal_sampler_;
    bool use_default_proposal_;
    
    // For Gaussian proposal
    std::vector<Real> center_;
    Real std_dev_;
    mutable std::normal_distribution<Real> normal_dist_;
};

/// \brief Quasi-Monte Carlo with randomization for error estimation
/// \details Combines low-discrepancy sequences with random shifts
template <typename Real>
class randomized_qmc_sampler : public monte_carlo_sampler<Real> {
public:
    using base_type = monte_carlo_sampler<Real>;
    using result_type = typename base_type::result_type;
    
    /// \brief Constructor
    /// \param dim Dimension
    /// \param n_randomizations Number of random shifts
    /// \param seed Random seed
    randomized_qmc_sampler(std::size_t dim, std::size_t n_randomizations = 10,
                          std::uint64_t seed = 0)
        : base_type(dim, seed), n_randomizations_(n_randomizations) {
        shifts_.resize(n_randomizations);
        for (auto& shift : shifts_) {
            shift.resize(dim);
            for (std::size_t d = 0; d < dim; ++d) {
                shift[d] = this->uniform_random();
            }
        }
    }
    
    /// \brief Sample using randomized QMC
    template <typename F>
    result_type sample(F&& f, const std::vector<Real>& lower,
                      const std::vector<Real>& upper, std::size_t n_samples) {
        if (lower.size() != this->dim_ || upper.size() != this->dim_) {
            throw std::invalid_argument("Dimension mismatch in bounds");
        }
        
        // Compute volume
        Real volume = Real(1);
        for (std::size_t i = 0; i < this->dim_; ++i) {
            volume *= (upper[i] - lower[i]);
        }
        
        // Run multiple randomized QMC estimates
        std::vector<Real> estimates;
        estimates.reserve(n_randomizations_);
        
        for (std::size_t r = 0; r < n_randomizations_; ++r) {
            Real sum = Real(0);
            std::vector<Real> point(this->dim_);
            
            for (std::size_t i = 0; i < n_samples; ++i) {
                // Generate QMC point with random shift
                generate_shifted_qmc_point(point, lower, upper, i, r);
                sum += f(point.data(), this->dim_);
            }
            
            estimates.push_back(sum / Real(n_samples));
        }
        
        // Compute mean and error from randomizations
        Real mean = std::accumulate(estimates.begin(), estimates.end(), Real(0)) 
                    / Real(n_randomizations_);
        
        Real var = Real(0);
        for (const auto& est : estimates) {
            Real diff = est - mean;
            var += diff * diff;
        }
        var /= Real(n_randomizations_ - 1);
        
        Real error = std::sqrt(var / Real(n_randomizations_));
        
        return {mean * volume, error * volume, var * volume * volume, 
                n_samples * n_randomizations_, volume};
    }
    
protected:
    void generate_point(std::vector<Real>& point,
                       const std::vector<Real>& lower,
                       const std::vector<Real>& upper,
                       std::size_t sample_idx,
                       std::size_t /*total_samples*/) override {
        generate_shifted_qmc_point(point, lower, upper, sample_idx, 0);
    }
    
private:
    /// \brief Generate QMC point with random shift (simplified Sobol-like)
    void generate_shifted_qmc_point(std::vector<Real>& point,
                                   const std::vector<Real>& lower,
                                   const std::vector<Real>& upper,
                                   std::size_t index,
                                   std::size_t shift_idx) {
        // Simplified radical inverse sequence (van der Corput)
        for (std::size_t d = 0; d < this->dim_; ++d) {
            Real qmc_val = radical_inverse(index + 1, get_prime(d));
            
            // Apply random shift with wraparound
            qmc_val = std::fmod(qmc_val + shifts_[shift_idx][d], Real(1));
            
            // Map to domain
            point[d] = lower[d] + (upper[d] - lower[d]) * qmc_val;
        }
    }
    
    /// \brief Compute radical inverse (van der Corput sequence)
    Real radical_inverse(std::size_t n, std::size_t base) const {
        Real result = Real(0);
        Real f = Real(1) / Real(base);
        
        while (n > 0) {
            result += (n % base) * f;
            n /= base;
            f /= Real(base);
        }
        
        return result;
    }
    
    /// \brief Get the d-th prime number (simplified)
    std::size_t get_prime(std::size_t d) const {
        static const std::size_t primes[] = {
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
            53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113
        };
        const std::size_t n_primes = sizeof(primes) / sizeof(primes[0]);
        return primes[d % n_primes];
    }
    
    std::size_t n_randomizations_;
    std::vector<std::vector<Real>> shifts_;
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_VARIANCE_REDUCTION_HPP