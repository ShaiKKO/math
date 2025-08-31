// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_PARALLEL_ADAPTIVE_HPP
#define BOOST_MATH_CUBATURE_DETAIL_PARALLEL_ADAPTIVE_HPP

/// \file parallel_adaptive.hpp
/// \brief Parallel adaptive integration with deterministic region processing
/// \details Implements parallel evaluation of regions in adaptive integration
///          while maintaining bitwise reproducibility across thread counts.

#include <boost/math/cubature/detail/parallel_executor.hpp>
#include <boost/math/cubature/detail/adaptivity.hpp>
#include <boost/math/cubature/detail/genz_malik_evaluator.hpp>
#include <boost/math/cubature/regions.hpp>
#include <boost/math/cubature/policies.hpp>
#include <vector>
#include <queue>
#include <mutex>
#include <atomic>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Parallel batch of regions for evaluation
template <typename Real>
struct region_batch {
    std::vector<region<Real>> regions;
    std::vector<embedded_pair_result<Real>> results;
    std::size_t thread_id;
    
    explicit region_batch(std::size_t capacity = 32) 
        : thread_id(0) {
        regions.reserve(capacity);
        results.reserve(capacity);
    }
};

/// \brief Thread-safe region accumulator with deterministic ordering
template <typename Real, typename Policy = default_policy>
class parallel_region_accumulator {
private:
    policy_accumulator<Real, Policy> integral_;
    policy_accumulator<Real, Policy> error_;
    mutable std::mutex mutex_;
    std::vector<region<Real>> processed_regions_;
    
public:
    explicit parallel_region_accumulator(const Policy& pol = Policy{})
        : integral_(pol), error_(pol) {}
    
    /// \brief Add regions from a batch (thread-safe)
    void add_batch(const region_batch<Real>& batch) {
        std::lock_guard<std::mutex> lock(mutex_);
        
        for (std::size_t i = 0; i < batch.regions.size(); ++i) {
            const auto& reg = batch.regions[i];
            integral_.add(reg.estimate_fine);
            error_.add(reg.error);
            processed_regions_.push_back(reg);
        }
    }
    
    /// \brief Get accumulated integral
    Real get_integral() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return integral_.sum();
    }
    
    /// \brief Get accumulated error
    Real get_error() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return error_.sum();
    }
    
    /// \brief Get all processed regions
    std::vector<region<Real>> get_regions() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return processed_regions_;
    }
};

/// \brief Parallel adaptive integrator
template <typename Real, typename F, typename Policy = default_policy>
class parallel_adaptive_integrator {
private:
    const F& f_;
    const hypercube<Real>& box_;
    Real abs_tol_;
    Real rel_tol_;
    std::size_t max_eval_;
    parallel_executor<Real, Policy> executor_;
    Policy policy_;
    
    // Evaluation state
    std::atomic<std::size_t> total_evals_{0};
    std::atomic<std::size_t> num_regions_{0};
    
public:
    parallel_adaptive_integrator(
        const F& f,
        const hypercube<Real>& box,
        Real abs_tol,
        Real rel_tol,
        std::size_t max_eval,
        const parallel_config& config,
        const Policy& pol = Policy{})
        : f_(f), box_(box),
          abs_tol_(abs_tol), rel_tol_(rel_tol),
          max_eval_(max_eval),
          executor_(config, pol),
          policy_(pol) {}
    
    /// \brief Run parallel adaptive integration
    result<Real> integrate() {
        const std::size_t dim = box_.dimension();
        
        // Initial region evaluation (sequential)
        region<Real> initial_region(dim);
        initial_region.a = box_.lower;
        initial_region.b = box_.upper;
        
        embedded_pair_result<Real> initial_result;
        bool success = evaluate_region(initial_region, initial_result);
        
        if (!success) {
            result<Real> res;
            res.status = status_code::dimension_error;
            res.error = std::numeric_limits<Real>::max();
            return res;
        }
        
        initial_region.estimate_fine = initial_result.estimate_fine;
        initial_region.estimate_coarse = initial_result.estimate_coarse;
        initial_region.error = initial_result.embedded_error;
        initial_region.evaluations = initial_result.evaluations;
        
        total_evals_ = initial_result.evaluations;
        num_regions_ = 1;
        
        // Check initial convergence
        Real current_integral = initial_region.estimate_fine;
        Real current_error = initial_region.error;
        
        if (current_error <= abs_tol_ + rel_tol_ * std::abs(current_integral)) {
            result<Real> res;
            res.value = current_integral;
            res.error = current_error;
            res.evaluations = total_evals_;
            res.status = status_code::success;
            return res;
        }
        
        // Priority queue for regions (by error)
        auto cmp = [](const region<Real>& a, const region<Real>& b) {
            return a.error < b.error;  // Max heap
        };
        std::priority_queue<region<Real>, std::vector<region<Real>>, decltype(cmp)> pq(cmp);
        pq.push(initial_region);
        
        // Parallel accumulator
        parallel_region_accumulator<Real, Policy> accumulator(policy_);
        
        // Main adaptive loop with parallel batch processing
        std::size_t batch_size = executor_.config().min_batch_size;
        
        while (!pq.empty() && total_evals_ < max_eval_) {
            // Collect regions for parallel processing
            std::vector<region<Real>> batch_regions;
            batch_regions.reserve(batch_size);
            
            // Extract up to batch_size regions with largest errors
            while (!pq.empty() && batch_regions.size() < batch_size) {
                batch_regions.push_back(pq.top());
                pq.pop();
            }
            
            // Process regions in parallel
            std::vector<region_batch<Real>> results = 
                process_regions_parallel(batch_regions);
            
            // Add results back to queue and accumulator
            for (const auto& batch : results) {
                accumulator.add_batch(batch);
                
                for (const auto& reg : batch.regions) {
                    pq.push(reg);
                }
            }
            
            // Check convergence
            current_integral = accumulator.get_integral();
            current_error = accumulator.get_error();
            
            if (current_error <= abs_tol_ + rel_tol_ * std::abs(current_integral)) {
                break;
            }
        }
        
        // Prepare final result
        result<Real> res;
        res.value = accumulator.get_integral();
        res.error = accumulator.get_error();
        res.evaluations = total_evals_;
        res.status = (res.error <= abs_tol_ + rel_tol_ * std::abs(res.value)) 
                     ? status_code::success 
                     : status_code::maxeval_reached;
        
        return res;
    }
    
private:
    /// \brief Process regions in parallel
    std::vector<region_batch<Real>> process_regions_parallel(
        const std::vector<region<Real>>& regions)
    {
        std::size_t n_threads = executor_.num_threads();
        
        // Partition regions among threads deterministically
        auto ranges = work_partitioner<int>::create_index_ranges(
            regions.size(), n_threads);
        
        std::vector<std::future<region_batch<Real>>> futures;
        futures.reserve(ranges.size());
        
        auto* pool = executor_.get_pool();
        
        for (std::size_t t = 0; t < ranges.size(); ++t) {
            auto range = ranges[t];
            auto start = range.first;
            auto end = range.second;
            
            auto promise = std::make_shared<std::promise<region_batch<Real>>>();
            futures.push_back(promise->get_future());
            
            // Copy regions for this thread
            std::vector<region<Real>> thread_regions(
                regions.begin() + start, 
                regions.begin() + end);
            
            // Process on thread pool
            if (pool) {
                boost::asio::post(*pool, 
                    [this, thread_regions = std::move(thread_regions), t, promise]() {
                    region_batch<Real> batch;
                    batch.thread_id = t;
                    
                    for (const auto& reg : thread_regions) {
                        // Bisect region
                        auto bisect_result = bisect_region(reg, 
                            select_split_dimension(reg));
                        auto left = bisect_result.first;
                        auto right = bisect_result.second;
                        
                        // Evaluate children
                        embedded_pair_result<Real> left_result, right_result;
                        evaluate_region(left, left_result);
                        evaluate_region(right, right_result);
                        
                        // Update regions
                        left.estimate_fine = left_result.estimate_fine;
                        left.estimate_coarse = left_result.estimate_coarse;
                        left.error = left_result.embedded_error;
                        left.evaluations = left_result.evaluations;
                        
                        right.estimate_fine = right_result.estimate_fine;
                        right.estimate_coarse = right_result.estimate_coarse;
                        right.error = right_result.embedded_error;
                        right.evaluations = right_result.evaluations;
                        
                        batch.regions.push_back(left);
                        batch.regions.push_back(right);
                        
                        total_evals_ += left.evaluations + right.evaluations;
                        num_regions_ += 2;
                    }
                    
                    promise->set_value(std::move(batch));
                });
            } else {
                // Sequential fallback
                region_batch<Real> batch;
                batch.thread_id = t;
                
                for (const auto& reg : thread_regions) {
                    auto bisect_result = bisect_region(reg, 
                        select_split_dimension(reg));
                    auto left = bisect_result.first;
                    auto right = bisect_result.second;
                    
                    embedded_pair_result<Real> left_result, right_result;
                    evaluate_region(left, left_result);
                    evaluate_region(right, right_result);
                    
                    left.estimate_fine = left_result.estimate_fine;
                    left.estimate_coarse = left_result.estimate_coarse;
                    left.error = left_result.embedded_error;
                    left.evaluations = left_result.evaluations;
                    
                    right.estimate_fine = right_result.estimate_fine;
                    right.estimate_coarse = right_result.estimate_coarse;
                    right.error = right_result.embedded_error;
                    right.evaluations = right_result.evaluations;
                    
                    batch.regions.push_back(left);
                    batch.regions.push_back(right);
                    
                    total_evals_ += left.evaluations + right.evaluations;
                    num_regions_ += 2;
                }
                
                promise->set_value(std::move(batch));
            }
        }
        
        // Collect results in deterministic order
        std::vector<region_batch<Real>> results;
        results.reserve(futures.size());
        
        for (auto& fut : futures) {
            results.push_back(fut.get());
        }
        
        return results;
    }
    
    /// \brief Evaluate a single region
    bool evaluate_region(region<Real>& reg, embedded_pair_result<Real>& result) {
        // Dispatch to appropriate evaluator based on dimension
        std::size_t dim = reg.dimension();
        
        if (dim <= 15) {
            // Use Genz-Malik evaluator with dimension dispatch
            switch (dim) {
                case 2:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 2>(
                        f_, reg, result, nullptr);
                case 3:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 3>(
                        f_, reg, result, nullptr);
                case 4:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 4>(
                        f_, reg, result, nullptr);
                case 5:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 5>(
                        f_, reg, result, nullptr);
                case 6:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 6>(
                        f_, reg, result, nullptr);
                case 7:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 7>(
                        f_, reg, result, nullptr);
                case 8:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 8>(
                        f_, reg, result, nullptr);
                case 9:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 9>(
                        f_, reg, result, nullptr);
                case 10:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 10>(
                        f_, reg, result, nullptr);
                case 11:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 11>(
                        f_, reg, result, nullptr);
                case 12:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 12>(
                        f_, reg, result, nullptr);
                case 13:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 13>(
                        f_, reg, result, nullptr);
                case 14:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 14>(
                        f_, reg, result, nullptr);
                case 15:
                    return genz_malik_evaluator<Real>::template evaluate_embedded_pair<F, 15>(
                        f_, reg, result, nullptr);
                default:
                    return false;
            }
        }
        
        return false;  // No rules available
    }
    
    /// \brief Select dimension to split
    std::size_t select_split_dimension(const region<Real>& reg) {
        // Simple strategy: split along longest dimension
        std::size_t best_dim = 0;
        Real max_width = 0;
        
        for (std::size_t i = 0; i < reg.dimension(); ++i) {
            Real width = reg.b[i] - reg.a[i];
            if (width > max_width) {
                max_width = width;
                best_dim = i;
            }
        }
        
        return best_dim;
    }
    
    /// \brief Bisect region along given dimension
    std::pair<region<Real>, region<Real>> bisect_region(
        const region<Real>& parent, 
        std::size_t split_dim)
    {
        region<Real> left = parent;
        region<Real> right = parent;
        
        Real split_point = (parent.a[split_dim] + parent.b[split_dim]) / Real(2);
        
        left.b[split_dim] = split_point;
        right.a[split_dim] = split_point;
        
        // Clear cached values as they're no longer valid
        left.cached_values.reset();
        right.cached_values.reset();
        
        return {left, right};
    }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_PARALLEL_ADAPTIVE_HPP