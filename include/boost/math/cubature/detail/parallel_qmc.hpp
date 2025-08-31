// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_PARALLEL_QMC_HPP
#define BOOST_MATH_CUBATURE_DETAIL_PARALLEL_QMC_HPP

/// \file parallel_qmc.hpp
/// \brief Parallel QMC integration with deterministic block partitioning
/// \details Implements parallel evaluation of QMC points with reproducible
///          Sobol index assignment and independent scrambling per thread.

#include <boost/math/cubature/detail/parallel_executor.hpp>
#include <boost/math/cubature/detail/sobol_owen.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/policies.hpp>
#include <boost/random/sobol.hpp>
#include <vector>
#include <future>
#include <cstdint>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief QMC work block for parallel processing
template <typename Real>
struct qmc_block {
    std::size_t start_index;     ///< Starting Sobol index
    std::size_t end_index;       ///< Ending Sobol index (exclusive)
    std::size_t thread_id;       ///< Thread processing this block
    std::uint32_t scramble_seed; ///< Seed for Owen scrambling
    Real partial_sum;            ///< Partial sum for this block
    std::size_t evaluations;     ///< Number of evaluations
    
    qmc_block(std::size_t start, std::size_t end, std::size_t tid)
        : start_index(start), end_index(end), thread_id(tid),
          scramble_seed(0), partial_sum(0), evaluations(0) {}
};

/// \brief Parallel QMC integrator with deterministic partitioning
template <typename Real, typename F, std::size_t Dim, typename Policy = default_policy>
class parallel_qmc_integrator {
private:
    const F& f_;
    const hypercube<Real>& box_;
    std::size_t n_points_;
    bool use_scrambling_;
    std::size_t n_replicates_;
    parallel_executor<Real, Policy> executor_;
    Policy policy_;
    
public:
    parallel_qmc_integrator(
        const F& f,
        const hypercube<Real>& box,
        std::size_t n_points,
        bool use_scrambling,
        std::size_t n_replicates,
        const parallel_config& config,
        const Policy& pol = Policy{})
        : f_(f), box_(box),
          n_points_(n_points),
          use_scrambling_(use_scrambling),
          n_replicates_(n_replicates),
          executor_(config, pol),
          policy_(pol) {}
    
    /// \brief Run parallel QMC integration
    result<Real> integrate() {
        std::vector<Real> replicate_values;
        replicate_values.reserve(n_replicates_);
        
        Real total_sum = 0;
        std::size_t total_evals = 0;
        
        // Process each replicate
        for (std::size_t rep = 0; rep < n_replicates_; ++rep) {
            Real replicate_sum = integrate_replicate(rep);
            replicate_values.push_back(replicate_sum);
            total_sum += replicate_sum;
            total_evals += n_points_;
        }
        
        // Prepare result
        result<Real> res;
        res.evaluations = total_evals;
        res.status = status_code::success;
        
        if (use_scrambling_ && n_replicates_ > 1) {
            res.value = total_sum / Real(n_replicates_);
            
            // Compute standard error
            Real mean = res.value;
            Real variance = 0;
            for (const auto& val : replicate_values) {
                Real diff = val - mean;
                variance += diff * diff;
            }
            variance /= Real(n_replicates_ - 1);
            res.error = Real(2) * std::sqrt(variance / Real(n_replicates_));
        } else {
            res.value = total_sum / Real(n_replicates_);
            res.error = 0;  // No error estimate without scrambling
        }
        
        return res;
    }
    
private:
    /// \brief Integrate a single QMC replicate in parallel
    Real integrate_replicate(std::size_t replicate_id) {
        std::size_t n_threads = executor_.num_threads();
        
        // Create deterministic block partitioning
        auto ranges = work_partitioner<int>::create_index_ranges(
            n_points_, n_threads);
        
        std::vector<qmc_block<Real>> blocks;
        blocks.reserve(ranges.size());
        
        // Initialize blocks with deterministic scramble seeds
        for (std::size_t t = 0; t < ranges.size(); ++t) {
            auto range = ranges[t];
            auto start = range.first;
            auto end = range.second;
            qmc_block<Real> block(start, end, t);
            
            if (use_scrambling_) {
                // Deterministic seed based on replicate and thread
                block.scramble_seed = static_cast<std::uint32_t>(
                    replicate_id * 65521 + t * 32749 + 17);
            }
            
            blocks.push_back(block);
        }
        
        // Process blocks in parallel
        auto results = process_blocks_parallel(blocks, replicate_id);
        
        // Tree reduction for deterministic accumulation
        std::vector<Real> partial_sums;
        partial_sums.reserve(results.size());
        
        for (const auto& block : results) {
            partial_sums.push_back(block.partial_sum);
        }
        
        tree_reducer<Real> reducer;
        Real total = reducer.reduce(partial_sums);
        
        // Scale by volume / n_points
        Real volume = 1.0;
        for (std::size_t d = 0; d < Dim; ++d) {
            volume *= (box_.upper[d] - box_.lower[d]);
        }
        
        return (total * volume) / Real(n_points_);
    }
    
    /// \brief Process QMC blocks in parallel
    std::vector<qmc_block<Real>> process_blocks_parallel(
        const std::vector<qmc_block<Real>>& blocks,
        std::size_t replicate_id)
    {
        std::vector<std::future<qmc_block<Real>>> futures;
        futures.reserve(blocks.size());
        
        auto* pool = executor_.get_pool();
        
        for (const auto& block : blocks) {
            auto promise = std::make_shared<std::promise<qmc_block<Real>>>();
            futures.push_back(promise->get_future());
            
            if (pool) {
                boost::asio::post(*pool, 
                    [this, block, replicate_id, promise]() {
                    auto result = evaluate_block(block, replicate_id);
                    promise->set_value(result);
                });
            } else {
                // Sequential fallback
                auto result = evaluate_block(block, replicate_id);
                promise->set_value(result);
            }
        }
        
        // Collect results in deterministic order
        std::vector<qmc_block<Real>> results;
        results.reserve(futures.size());
        
        for (auto& fut : futures) {
            results.push_back(fut.get());
        }
        
        return results;
    }
    
    /// \brief Evaluate a single QMC block
    qmc_block<Real> evaluate_block(
        const qmc_block<Real>& block,
        std::size_t replicate_id)
    {
        qmc_block<Real> result = block;
        
        // Create Sobol engine and skip to block's start index
        boost::random::sobol_engine<std::uint64_t, Dim> sobol(0);
        sobol.discard(block.start_index);
        
        Real sum = 0;
        std::vector<Real> point(Dim);
        
        for (std::size_t i = block.start_index; i < block.end_index; ++i) {
            // Generate Sobol point
            auto sobol_point = sobol();
            
            // Convert to [0,1)^d and apply scrambling if needed
            for (std::size_t d = 0; d < Dim; ++d) {
                Real u = Real(sobol_point[d]) / 
                        Real(std::numeric_limits<std::uint64_t>::max());
                
                if (use_scrambling_) {
                    u = owen_scramble(u, block.scramble_seed + 
                                        static_cast<std::uint32_t>(d));
                }
                
                // Apply tent transform for variance reduction
                u = tent_transform(u);
                
                // Map from [0,1] to [a,b]
                point[d] = box_.lower[d] + u * (box_.upper[d] - box_.lower[d]);
            }
            
            // Evaluate integrand
            Real value = evaluate_integrand(point);
            sum += value;
            result.evaluations++;
        }
        
        result.partial_sum = sum;
        return result;
    }
    
    // SFINAE helper for function signature
    template <typename Func, typename Point, typename = void>
    struct integrand_sig {
        static constexpr int value = 0;
    };
    
    template <typename Func, typename Point>
    struct integrand_sig<Func, Point,
        typename std::enable_if<
            std::is_convertible<decltype(std::declval<Func>()(std::declval<Point>())), Real>::value
        >::type> {
        static constexpr int value = 1; // vector
    };
    
    template <typename Func, typename Point>
    struct integrand_sig<Func, Point,
        typename std::enable_if<
            !std::is_convertible<decltype(std::declval<Func>()(std::declval<Point>())), Real>::value &&
            std::is_convertible<decltype(std::declval<Func>()(std::declval<const Real*>(), std::declval<std::size_t>())), Real>::value
        >::type> {
        static constexpr int value = 2; // pointer + size
    };
    
    template <typename Func, typename Point>
    struct integrand_sig<Func, Point,
        typename std::enable_if<
            !std::is_convertible<decltype(std::declval<Func>()(std::declval<Point>())), Real>::value &&
            !std::is_convertible<decltype(std::declval<Func>()(std::declval<const Real*>(), std::declval<std::size_t>())), Real>::value &&
            std::is_convertible<decltype(std::declval<Func>()(std::declval<const Real*>())), Real>::value
        >::type> {
        static constexpr int value = 3; // pointer only
    };
    
    // Tag dispatch evaluation
    Real evaluate_dispatch(const std::vector<Real>& point, std::integral_constant<int, 0>) const {
        static_assert(sizeof(F) == 0, "Integrand must be callable");
        return Real(0);
    }
    
    Real evaluate_dispatch(const std::vector<Real>& point, std::integral_constant<int, 1>) const {
        return f_(point);
    }
    
    Real evaluate_dispatch(const std::vector<Real>& point, std::integral_constant<int, 2>) const {
        return f_(point.data(), Dim);
    }
    
    Real evaluate_dispatch(const std::vector<Real>& point, std::integral_constant<int, 3>) const {
        return f_(point.data());
    }
    
    /// \brief Evaluate integrand with proper signature detection
    Real evaluate_integrand(const std::vector<Real>& point) const {
        return evaluate_dispatch(point, std::integral_constant<int, integrand_sig<F, decltype(point)>::value>());
    }
    
    /// \brief Tent transform for variance reduction
    static Real tent_transform(Real u) {
        return Real(1) - Real(2) * std::abs(u - Real(0.5));
    }
    
    /// \brief Owen scrambling for randomized QMC
    static Real owen_scramble(Real u, std::uint32_t seed) {
        // Use full nested uniform Owen scrambling from sobol_owen.hpp
        return detail::owen_scramble<Real>(u, seed);
    }
};

/// \brief Factory function for parallel QMC integration
template <typename Real, typename F, std::size_t Dim, typename Policy>
result<Real> integrate_qmc_parallel(
    const F& f,
    const hypercube<Real>& box,
    std::size_t n_points,
    bool use_scrambling,
    std::size_t n_replicates,
    const parallel_config& config,
    const Policy& pol = Policy{})
{
    parallel_qmc_integrator<Real, F, Dim, Policy> integrator(
        f, box, n_points, use_scrambling, n_replicates, config, pol);
    return integrator.integrate();
}

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_PARALLEL_QMC_HPP