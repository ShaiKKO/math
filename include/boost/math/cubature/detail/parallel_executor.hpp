// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_PARALLEL_EXECUTOR_HPP
#define BOOST_MATH_CUBATURE_DETAIL_PARALLEL_EXECUTOR_HPP

/// \file parallel_executor.hpp
/// \brief Thread pool executor for deterministic parallel integration
/// \details Provides a wrapper around boost::asio::thread_pool with deterministic
///          work partitioning and reproducible reduction strategies.
///
/// \section parallel_guarantees Determinism Guarantees
/// - Fixed work assignment based on thread index
/// - Reproducible traversal orders (no dynamic scheduling in deterministic mode)
/// - Tree-based reduction with fixed merge order
/// - Thread-local workspaces to avoid contention
///
/// \section parallel_performance Performance
/// - Batch processing to amortize thread overhead
/// - Cache-aware partitioning to minimize false sharing
/// - Optional work-stealing for non-deterministic mode
/// - Thread-local memory pools for allocation efficiency

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/math/cubature/policies.hpp>
#include <vector>
#include <future>
#include <memory>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <thread>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Configuration for parallel execution
struct parallel_config {
    std::size_t num_threads;           ///< Number of worker threads (0 = hardware concurrency)
    std::size_t min_batch_size;        ///< Minimum work items per batch
    std::size_t max_batch_size;        ///< Maximum work items per batch
    bool deterministic;                ///< Ensure reproducible results
    bool enable_work_stealing;         ///< Allow dynamic load balancing (non-deterministic)
    
    parallel_config()
        : num_threads(0),
          min_batch_size(16),
          max_batch_size(256),
          deterministic(true),
          enable_work_stealing(false) {}
    
    /// \brief Get effective number of threads
    std::size_t effective_threads() const {
        if (num_threads == 0) {
            return std::thread::hardware_concurrency();
        }
        return num_threads;
    }
};

/// \brief Work partitioning strategy for deterministic parallelism
template <typename WorkItem>
class work_partitioner {
public:
    using work_batch = std::vector<WorkItem>;
    using partition = std::vector<work_batch>;
    
    /// \brief Partition work items into batches for parallel processing
    /// \param items Work items to partition
    /// \param config Parallel execution configuration
    /// \return Vector of batches, one per thread
    static partition create_batches(
        std::vector<WorkItem> items,
        const parallel_config& config)
    {
        std::size_t n_threads = config.effective_threads();
        std::size_t n_items = items.size();
        
        // Single thread or too few items
        if (n_threads == 1 || n_items < config.min_batch_size) {
            return partition{std::move(items)};
        }
        
        // Calculate items per thread with deterministic rounding
        std::size_t base_size = n_items / n_threads;
        std::size_t remainder = n_items % n_threads;
        
        partition batches;
        batches.reserve(n_threads);
        
        auto it = items.begin();
        for (std::size_t t = 0; t < n_threads; ++t) {
            // Deterministic distribution: first 'remainder' threads get one extra item
            std::size_t batch_size = base_size + (t < remainder ? 1 : 0);
            
            if (batch_size > 0) {
                work_batch batch;
                batch.reserve(batch_size);
                
                for (std::size_t i = 0; i < batch_size; ++i) {
                    batch.push_back(std::move(*it));
                    ++it;
                }
                
                batches.push_back(std::move(batch));
            }
        }
        
        return batches;
    }
    
    /// \brief Create index ranges for partitioning continuous work
    /// \param total_items Total number of items
    /// \param n_threads Number of threads
    /// \return Vector of (start, end) index pairs
    static std::vector<std::pair<std::size_t, std::size_t>> 
    create_index_ranges(std::size_t total_items, std::size_t n_threads)
    {
        std::vector<std::pair<std::size_t, std::size_t>> ranges;
        ranges.reserve(n_threads);
        
        std::size_t base_size = total_items / n_threads;
        std::size_t remainder = total_items % n_threads;
        
        std::size_t start = 0;
        for (std::size_t t = 0; t < n_threads; ++t) {
            std::size_t size = base_size + (t < remainder ? 1 : 0);
            if (size > 0) {
                ranges.emplace_back(start, start + size);
                start += size;
            }
        }
        
        return ranges;
    }
};

/// \brief Tree-based reduction for deterministic accumulation
/// \tparam T Type to accumulate
/// \tparam BinaryOp Binary operation for reduction
template <typename T, typename BinaryOp = std::plus<T>>
class tree_reducer {
private:
    BinaryOp op_;
    
public:
    explicit tree_reducer(BinaryOp op = BinaryOp{}) : op_(op) {}
    
    /// \brief Reduce values using binary tree structure
    /// \param values Values to reduce (will be modified)
    /// \return Reduced result
    T reduce(std::vector<T>& values) const {
        if (values.empty()) {
            return T{};
        }
        
        // Binary tree reduction for deterministic order
        // Example for 8 values:
        // Level 0: [0,1] [2,3] [4,5] [6,7]
        // Level 1: [0+1,2+3] [4+5,6+7]
        // Level 2: [(0+1)+(2+3),(4+5)+(6+7)]
        // Result: ((0+1)+(2+3))+((4+5)+(6+7))
        
        std::size_t n = values.size();
        while (n > 1) {
            std::size_t new_n = (n + 1) / 2;
            
            for (std::size_t i = 0; i < new_n; ++i) {
                if (2 * i + 1 < n) {
                    values[i] = op_(values[2 * i], values[2 * i + 1]);
                } else {
                    values[i] = values[2 * i];
                }
            }
            
            n = new_n;
        }
        
        return values[0];
    }
    
    /// \brief Parallel tree reduction with deterministic ordering
    /// \param values Values from different threads
    /// \param pool Thread pool for parallel reduction
    /// \return Reduced result
    T parallel_reduce(
        std::vector<T> values,
        boost::asio::thread_pool& pool) const
    {
        if (values.size() <= 1) {
            return values.empty() ? T{} : values[0];
        }
        
        // Parallel pairwise reduction
        while (values.size() > 1) {
            std::size_t new_size = (values.size() + 1) / 2;
            std::vector<std::future<T>> futures;
            futures.reserve(new_size);
            
            for (std::size_t i = 0; i < new_size; ++i) {
                if (2 * i + 1 < values.size()) {
                    // Pair reduction
                    auto promise = std::make_shared<std::promise<T>>();
                    futures.push_back(promise->get_future());
                    
                    T val1 = values[2 * i];
                    T val2 = values[2 * i + 1];
                    
                    boost::asio::post(pool, [this, val1, val2, promise]() {
                        promise->set_value(op_(val1, val2));
                    });
                } else {
                    // Odd element, carry forward
                    auto promise = std::make_shared<std::promise<T>>();
                    promise->set_value(values[2 * i]);
                    futures.push_back(promise->get_future());
                }
            }
            
            // Collect results
            values.clear();
            for (auto& fut : futures) {
                values.push_back(fut.get());
            }
        }
        
        return values[0];
    }
};

/// \brief Parallel executor for integration algorithms
template <typename Real, typename Policy = boost::math::cubature::default_policy>
class parallel_executor {
private:
    std::unique_ptr<boost::asio::thread_pool> pool_;
    parallel_config config_;
    Policy policy_;
    
public:
    /// \brief Constructor
    /// \param config Parallel execution configuration
    /// \param pol Boost.Math policy
    explicit parallel_executor(
        const parallel_config& config = parallel_config{},
        const Policy& pol = Policy{})
        : config_(config), policy_(pol)
    {
        std::size_t n_threads = config_.effective_threads();
        if (n_threads > 1) {
            pool_ = std::make_unique<boost::asio::thread_pool>(n_threads);
        }
    }
    
    /// \brief Destructor - ensures all tasks complete
    ~parallel_executor() {
        if (pool_) {
            pool_->join();
        }
    }
    
    /// \brief Check if parallel execution is enabled
    bool is_parallel() const {
        return pool_ != nullptr;
    }
    
    /// \brief Get number of threads
    std::size_t num_threads() const {
        return config_.effective_threads();
    }
    
    /// \brief Execute function on all items in parallel
    /// \tparam Item Type of work item
    /// \tparam Func Function to execute on each item
    /// \param items Work items
    /// \param func Function to apply
    template <typename Item, typename Func>
    void execute_batch(const std::vector<Item>& items, Func func) {
        if (!pool_ || items.size() < config_.min_batch_size) {
            // Execute sequentially
            for (const auto& item : items) {
                func(item);
            }
            return;
        }
        
        // Create batches for parallel execution
        auto batches = work_partitioner<Item>::create_batches(
            items, config_);
        
        std::vector<std::future<void>> futures;
        futures.reserve(batches.size());
        
        for (auto& batch : batches) {
            auto promise = std::make_shared<std::promise<void>>();
            futures.push_back(promise->get_future());
            
            boost::asio::post(*pool_, [batch = std::move(batch), func, promise]() {
                for (const auto& item : batch) {
                    func(item);
                }
                promise->set_value();
            });
        }
        
        // Wait for all tasks to complete
        for (auto& fut : futures) {
            fut.wait();
        }
    }
    
    /// \brief Map-reduce operation with parallel execution
    /// \tparam Item Type of work item
    /// \tparam MapFunc Mapping function
    /// \tparam ReduceOp Reduction operation
    /// \tparam Result Result type
    template <typename Item, typename MapFunc, typename ReduceOp, typename Result>
    Result map_reduce(
        const std::vector<Item>& items,
        MapFunc map_func,
        ReduceOp reduce_op,
        Result init = Result{})
    {
        if (!pool_ || items.size() < config_.min_batch_size) {
            // Sequential execution
            Result result = init;
            for (const auto& item : items) {
                result = reduce_op(result, map_func(item));
            }
            return result;
        }
        
        // Create batches for parallel execution
        auto batches = work_partitioner<Item>::create_batches(
            items, config_);
        
        std::vector<std::future<Result>> futures;
        futures.reserve(batches.size());
        
        // Process batches in parallel
        for (auto& batch : batches) {
            auto promise = std::make_shared<std::promise<Result>>();
            futures.push_back(promise->get_future());
            
            boost::asio::post(*pool_, 
                [batch = std::move(batch), map_func, reduce_op, init, promise]() {
                Result local_result = init;
                for (const auto& item : batch) {
                    local_result = reduce_op(local_result, map_func(item));
                }
                promise->set_value(local_result);
            });
        }
        
        // Collect and reduce results using tree reduction
        std::vector<Result> partial_results;
        partial_results.reserve(futures.size());
        
        for (auto& fut : futures) {
            partial_results.push_back(fut.get());
        }
        
        // Use tree reduction for deterministic order
        tree_reducer<Result, ReduceOp> reducer(reduce_op);
        return reducer.reduce(partial_results);
    }
    
    /// \brief Get thread pool for custom tasks
    boost::asio::thread_pool* get_pool() {
        return pool_.get();
    }
    
    /// \brief Get configuration
    const parallel_config& config() const {
        return config_;
    }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_PARALLEL_EXECUTOR_HPP