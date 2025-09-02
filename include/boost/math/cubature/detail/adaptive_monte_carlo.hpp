// Copyright 2025 Boost.Math Contributors
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_ADAPTIVE_MONTE_CARLO_HPP
#define BOOST_MATH_CUBATURE_DETAIL_ADAPTIVE_MONTE_CARLO_HPP

#include <boost/math/cubature/detail/monte_carlo_base.hpp>
#include <boost/math/cubature/detail/variance_reduction.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <chrono>
#include <deque>
#include <memory>
#include <thread>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Running statistics for adaptive Monte Carlo
/// \details Maintains online statistics with numerical stability
template <typename Real>
class running_statistics {
public:
    running_statistics() : n_(0), mean_(0), M2_(0), min_(std::numeric_limits<Real>::max()),
                          max_(std::numeric_limits<Real>::lowest()) {}
    
    /// \brief Add a new sample
    void add_sample(Real x) {
        ++n_;
        Real delta = x - mean_;
        mean_ += delta / Real(n_);
        Real delta2 = x - mean_;
        M2_ += delta * delta2;
        
        min_ = std::min(min_, x);
        max_ = std::max(max_, x);
    }
    
    /// \brief Add multiple samples
    void add_samples(const std::vector<Real>& samples) {
        for (const auto& x : samples) {
            add_sample(x);
        }
    }
    
    /// \brief Get current mean
    Real mean() const { return mean_; }
    
    /// \brief Get sample variance
    Real variance() const {
        return n_ > 1 ? M2_ / Real(n_ - 1) : Real(0);
    }
    
    /// \brief Get standard deviation
    Real std_dev() const { return std::sqrt(variance()); }
    
    /// \brief Get standard error of the mean
    Real std_error() const {
        return n_ > 0 ? std_dev() / std::sqrt(Real(n_)) : std::numeric_limits<Real>::infinity();
    }
    
    /// \brief Get sample count
    std::size_t count() const { return n_; }
    
    /// \brief Get minimum value seen
    Real min() const { return min_; }
    
    /// \brief Get maximum value seen
    Real max() const { return max_; }
    
    /// \brief Reset statistics
    void reset() {
        n_ = 0;
        mean_ = Real(0);
        M2_ = Real(0);
        min_ = std::numeric_limits<Real>::max();
        max_ = std::numeric_limits<Real>::lowest();
    }
    
    /// \brief Combine with another statistics object
    void merge(const running_statistics& other) {
        if (other.n_ == 0) return;
        if (n_ == 0) {
            *this = other;
            return;
        }
        
        Real combined_n = Real(n_ + other.n_);
        Real delta = other.mean_ - mean_;
        
        // Combined mean
        Real new_mean = (Real(n_) * mean_ + Real(other.n_) * other.mean_) / combined_n;
        
        // Combined variance (parallel algorithm)
        Real new_M2 = M2_ + other.M2_ + delta * delta * Real(n_) * Real(other.n_) / combined_n;
        
        mean_ = new_mean;
        M2_ = new_M2;
        n_ += other.n_;
        min_ = std::min(min_, other.min_);
        max_ = std::max(max_, other.max_);
    }
    
private:
    std::size_t n_;
    Real mean_;
    Real M2_;  // Sum of squared differences from mean
    Real min_;
    Real max_;
};

/// \brief Convergence monitor for adaptive sampling
template <typename Real>
class convergence_monitor {
public:
    /// \brief Constructor
    /// \param window_size Size of moving window for convergence check
    /// \param check_interval How often to check convergence
    convergence_monitor(std::size_t window_size = 100, std::size_t check_interval = 50)
        : window_size_(window_size), check_interval_(check_interval),
          samples_since_check_(0), is_converged_(false) {}
    
    /// \brief Add new estimate and check convergence
    /// \param value Current integral estimate
    /// \param error Current error estimate
    /// \param abs_tol Absolute tolerance
    /// \param rel_tol Relative tolerance
    /// \return true if converged
    bool update(Real value, Real error, Real abs_tol, Real rel_tol) {
        estimates_.push_back(value);
        errors_.push_back(error);
        
        // Maintain window size
        if (estimates_.size() > window_size_) {
            estimates_.pop_front();
            errors_.pop_front();
        }
        
        ++samples_since_check_;
        
        // Check convergence periodically
        if (samples_since_check_ >= check_interval_) {
            samples_since_check_ = 0;
            is_converged_ = check_convergence(abs_tol, rel_tol);
        }
        
        return is_converged_;
    }
    
    /// \brief Check if estimates have stabilized
    bool has_stabilized() const {
        if (estimates_.size() < window_size_ / 2) return false;
        
        // Check coefficient of variation of recent estimates
        running_statistics<Real> stats;
        for (const auto& est : estimates_) {
            stats.add_sample(est);
        }
        
        Real cv = stats.std_dev() / std::abs(stats.mean());
        return cv < Real(0.01);  // Less than 1% variation
    }
    
    /// \brief Get convergence rate estimate
    Real convergence_rate() const {
        if (errors_.size() < 10) return Real(0);
        
        // Estimate convergence rate from error decay
        std::size_t n = errors_.size();
        std::size_t half = n / 2;
        
        Real early_error = Real(0);
        Real late_error = Real(0);
        
        for (std::size_t i = 0; i < half; ++i) {
            early_error += errors_[i];
            late_error += errors_[i + half];
        }
        
        early_error /= Real(half);
        late_error /= Real(half);
        
        if (early_error <= Real(0) || late_error <= Real(0)) return Real(0);
        
        return std::log(early_error / late_error) / std::log(Real(2));
    }
    
    /// \brief Reset convergence monitor
    void reset() {
        estimates_.clear();
        errors_.clear();
        samples_since_check_ = 0;
        is_converged_ = false;
    }
    
private:
    bool check_convergence(Real abs_tol, Real rel_tol) const {
        if (estimates_.empty() || errors_.empty()) return false;
        
        Real current_value = estimates_.back();
        Real current_error = errors_.back();
        
        // Basic tolerance check
        bool tol_met = current_error < abs_tol || 
                      current_error < std::abs(current_value) * rel_tol;
        
        if (!tol_met) return false;
        
        // Additional stability check
        return has_stabilized();
    }
    
    std::size_t window_size_;
    std::size_t check_interval_;
    std::size_t samples_since_check_;
    bool is_converged_;
    std::deque<Real> estimates_;
    std::deque<Real> errors_;
};

/// \brief Adaptive Monte Carlo integrator with automatic convergence
template <typename Real>
class adaptive_monte_carlo {
public:
    using result_type = mc_result<Real>;
    
    /// \brief Sampling strategy enum
    enum class strategy {
        uniform,           ///< Standard Monte Carlo
        stratified,        ///< Stratified sampling
        antithetic,        ///< Antithetic variates
        latin_hypercube,   ///< Latin hypercube sampling
        importance,        ///< Importance sampling
        randomized_qmc     ///< Randomized quasi-Monte Carlo
    };
    
    /// \brief Constructor
    /// \param dim Dimension of integration domain
    /// \param strat Sampling strategy to use
    /// \param seed Random seed
    adaptive_monte_carlo(std::size_t dim, strategy strat = strategy::uniform,
                        std::uint64_t seed = 0)
        : dim_(dim), strategy_(strat), seed_(seed),
          min_samples_(100), max_samples_(1000000),
          batch_size_(1000), confidence_level_(Real(0.95)) {
        create_sampler();
    }
    
    /// \brief Integrate with adaptive sampling
    /// \param f Function to integrate
    /// \param lower Lower bounds
    /// \param upper Upper bounds
    /// \param abs_tol Absolute tolerance
    /// \param rel_tol Relative tolerance
    /// \param max_time_ms Maximum time in milliseconds (0 for no limit)
    template <typename F>
    result_type integrate(F&& f, const std::vector<Real>& lower,
                         const std::vector<Real>& upper,
                         Real abs_tol, Real rel_tol,
                         std::size_t max_time_ms = 0) {
        if (lower.size() != dim_ || upper.size() != dim_) {
            throw std::invalid_argument("Dimension mismatch in bounds");
        }
        
        auto start_time = std::chrono::steady_clock::now();
        
        // Initialize statistics and convergence monitor
        running_statistics<Real> stats;
        convergence_monitor<Real> monitor(100, batch_size_ / 10);
        
        // Compute volume
        Real volume = Real(1);
        for (std::size_t i = 0; i < dim_; ++i) {
            volume *= (upper[i] - lower[i]);
        }
        
        std::size_t total_samples = 0;
        std::vector<Real> point(dim_);
        
        // Initial batch
        std::vector<Real> batch_values;
        batch_values.reserve(batch_size_);
        
        while (total_samples < max_samples_) {
            // Check time limit
            if (max_time_ms > 0) {
                auto elapsed = std::chrono::steady_clock::now() - start_time;
                auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
                if (elapsed_ms > static_cast<long>(max_time_ms)) {
                    break;
                }
            }
            
            // Adaptive batch size
            std::size_t current_batch = batch_size_;
            if (total_samples == 0) {
                current_batch = min_samples_;  // Start with minimum
            } else if (monitor.has_stabilized()) {
                // Increase batch size when stable
                current_batch = std::min(batch_size_ * 2, 
                                        max_samples_ - total_samples);
            }
            
            // Generate batch of samples
            batch_values.clear();
            for (std::size_t i = 0; i < current_batch; ++i) {
                sampler_->generate_point(point, lower, upper, 
                                        total_samples + i, total_samples + current_batch);
                batch_values.push_back(f(point.data(), dim_));
            }
            
            // Update statistics
            stats.add_samples(batch_values);
            total_samples += current_batch;
            
            // Check convergence
            Real current_mean = stats.mean() * volume;
            Real current_error = stats.std_error() * volume;
            
            if (monitor.update(current_mean, current_error, abs_tol, rel_tol)) {
                break;  // Converged
            }
            
            // Early termination if error is very small
            if (current_error < abs_tol * Real(0.1) && total_samples >= min_samples_) {
                break;
            }
            
            // Adaptive strategy switching (optional)
            if (total_samples > 10000 && !monitor.has_stabilized()) {
                consider_strategy_switch(stats);
            }
        }
        
        // Final result
        result_type result;
        result.value = stats.mean() * volume;
        result.error = stats.std_error() * volume;
        result.variance = stats.variance() * volume * volume;
        result.n_samples = total_samples;
        result.volume = volume;
        
        return result;
    }
    
    /// \brief Set sampling parameters
    void set_parameters(std::size_t min_samples, std::size_t max_samples,
                       std::size_t batch_size) {
        min_samples_ = min_samples;
        max_samples_ = max_samples;
        batch_size_ = batch_size;
    }
    
    /// \brief Set confidence level for error estimates
    void set_confidence_level(Real level) {
        if (level <= Real(0) || level >= Real(1)) {
            throw std::invalid_argument("Confidence level must be in (0,1)");
        }
        confidence_level_ = level;
    }
    
    /// \brief Get current sampling strategy
    strategy get_strategy() const { return strategy_; }
    
    /// \brief Switch sampling strategy
    void set_strategy(strategy strat) {
        if (strat != strategy_) {
            strategy_ = strat;
            create_sampler();
        }
    }
    
private:
    /// \brief Create appropriate sampler based on strategy
    void create_sampler() {
        switch (strategy_) {
            case strategy::uniform:
                sampler_ = std::make_unique<uniform_sampler<Real>>(dim_, seed_);
                break;
                
            case strategy::stratified:
                sampler_ = std::make_unique<stratified_sampler<Real>>(dim_, 0, seed_);
                break;
                
            case strategy::antithetic:
                sampler_ = std::make_unique<antithetic_sampler<Real>>(dim_, seed_);
                break;
                
            case strategy::latin_hypercube:
                sampler_ = std::make_unique<latin_hypercube_sampler<Real>>(dim_, seed_);
                break;
                
            case strategy::importance:
                sampler_ = std::make_unique<importance_sampler<Real>>(dim_, seed_);
                break;
                
            case strategy::randomized_qmc:
                sampler_ = std::make_unique<randomized_qmc_sampler<Real>>(dim_, 10, seed_);
                break;
                
            default:
                sampler_ = std::make_unique<uniform_sampler<Real>>(dim_, seed_);
        }
    }
    
    /// \brief Consider switching strategy based on statistics
    void consider_strategy_switch(const running_statistics<Real>& stats) {
        // Simple heuristic: switch to stratified if variance is high
        Real cv = stats.std_dev() / std::abs(stats.mean());
        
        if (cv > Real(1) && strategy_ == strategy::uniform) {
            // High variance, try stratified sampling
            set_strategy(strategy::stratified);
        } else if (cv > Real(2) && strategy_ == strategy::stratified) {
            // Very high variance, try antithetic variates
            set_strategy(strategy::antithetic);
        }
    }
    
    /// \brief Compute z-score for confidence level
    Real compute_z_score(Real confidence) const {
        // Approximate inverse normal CDF
        if (confidence <= Real(0.9)) return Real(1.645);
        if (confidence <= Real(0.95)) return Real(1.96);
        if (confidence <= Real(0.99)) return Real(2.576);
        return Real(3.0);  // Very high confidence
    }
    
    std::size_t dim_;
    strategy strategy_;
    std::uint64_t seed_;
    std::unique_ptr<monte_carlo_sampler<Real>> sampler_;
    
    std::size_t min_samples_;
    std::size_t max_samples_;
    std::size_t batch_size_;
    Real confidence_level_;
};

/// \brief Parallel adaptive Monte Carlo (stub for future implementation)
template <typename Real>
class parallel_adaptive_monte_carlo {
public:
    using result_type = mc_result<Real>;
    
    parallel_adaptive_monte_carlo(std::size_t dim, std::size_t n_threads = 0)
        : dim_(dim), n_threads_(n_threads ? n_threads : std::max(std::size_t(1), std::size_t(std::thread::hardware_concurrency()))) {}
    
    template <typename F>
    result_type integrate(F&& f, const std::vector<Real>& lower,
                         const std::vector<Real>& upper,
                         Real abs_tol, Real rel_tol) {
        // Parallel implementation would go here
        // For now, delegate to single-threaded version
        adaptive_monte_carlo<Real> mc(dim_);
        return mc.integrate(std::forward<F>(f), lower, upper, abs_tol, rel_tol);
    }
    
private:
    std::size_t dim_;
    std::size_t n_threads_;
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_ADAPTIVE_MONTE_CARLO_HPP