// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_MEMORY_POOL_HPP
#define BOOST_MATH_CUBATURE_DETAIL_MEMORY_POOL_HPP

#include <vector>
#include <memory>
#include <cstddef>
#include <new>
#include <utility>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Memory pool for efficient region allocation during adaptive integration
/// \details Implementation to reduce allocation overhead when managing
///          thousands of regions. Follows Boost.Pool design principles.
template <typename T>
class memory_pool {
public:
  /// \brief Configuration for pool behavior
  struct pool_config {
    std::size_t initial_capacity;     // Initial number of objects to pre-allocate
    std::size_t growth_factor;        // Multiplicative growth factor (e.g., 2)
    std::size_t max_capacity;         // Maximum pool size (0 = unlimited)
    bool allow_fallback;              // Use standard allocation if pool exhausted
    
    pool_config()
      : initial_capacity(1024),
        growth_factor(2),
        max_capacity(0),
        allow_fallback(true) {}
  };
  
private:
  /// \brief Single chunk of contiguous memory
  struct memory_chunk {
    std::unique_ptr<std::aligned_storage_t<sizeof(T), alignof(T)>[]> storage;
    std::size_t capacity;
    std::size_t used;
    
    explicit memory_chunk(std::size_t cap)
      : storage(std::make_unique<std::aligned_storage_t<sizeof(T), alignof(T)>[]>(cap)),
        capacity(cap),
        used(0) {}
    
    /// \brief Check if chunk has space for n objects
    bool has_space(std::size_t n = 1) const noexcept {
      return used + n <= capacity;
    }
    
    /// \brief Get pointer to next available slot
    T* get_next() noexcept {
      if (!has_space()) return nullptr;
      return reinterpret_cast<T*>(&storage[used++]);
    }
    
    /// \brief Reset chunk for reuse
    void reset() noexcept {
      used = 0;
    }
  };
  
  pool_config config_;
  std::vector<std::unique_ptr<memory_chunk>> chunks_;
  std::vector<T*> free_list_;           // List of deallocated objects for reuse
  std::size_t total_allocated_;         // Total objects allocated
  std::size_t total_deallocated_;       // Total objects deallocated
  std::size_t current_chunk_index_;     // Index of current allocation chunk
  
public:
  /// \brief Construct pool with given configuration
  explicit memory_pool(const pool_config& cfg = pool_config())
    : config_(cfg),
      total_allocated_(0),
      total_deallocated_(0),
      current_chunk_index_(0)
  {
    // Pre-allocate initial chunk
    if (config_.initial_capacity > 0) {
      chunks_.emplace_back(std::make_unique<memory_chunk>(config_.initial_capacity));
    }
  }
  
  /// \brief Destructor - no automatic cleanup of allocated objects
  /// \details User is responsible for calling destroy on all allocated objects
  ~memory_pool() = default;
  
  // Disable copy operations
  memory_pool(const memory_pool&) = delete;
  memory_pool& operator=(const memory_pool&) = delete;
  
  // Enable move operations
  memory_pool(memory_pool&&) noexcept = default;
  memory_pool& operator=(memory_pool&&) noexcept = default;
  
  /// \brief Allocate space for one object (uninitialized)
  T* allocate() {
    // First check free list for reusable objects
    if (!free_list_.empty()) {
      T* ptr = free_list_.back();
      free_list_.pop_back();
      total_deallocated_--;
      return ptr;
    }
    
    // Find a chunk with available space
    for (std::size_t i = current_chunk_index_; i < chunks_.size(); ++i) {
      if (chunks_[i]->has_space()) {
        current_chunk_index_ = i;
        total_allocated_++;
        return chunks_[i]->get_next();
      }
    }
    
    // All chunks full, need to allocate new chunk
    if (!grow_pool()) {
      // Pool growth failed, use fallback if allowed
      if (config_.allow_fallback) {
        return static_cast<T*>(::operator new(sizeof(T)));
      }
      return nullptr;
    }
    
    // Retry allocation with new chunk
    total_allocated_++;
    return chunks_[current_chunk_index_]->get_next();
  }
  
  /// \brief Construct object in pool with given arguments
  template <typename... Args>
  T* construct(Args&&... args) {
    T* ptr = allocate();
    if (ptr) {
      try {
        new(ptr) T(std::forward<Args>(args)...);
      } catch (...) {
        deallocate(ptr);
        throw;
      }
    }
    return ptr;
  }
  
  /// \brief Deallocate object (does not call destructor)
  void deallocate(T* ptr) noexcept {
    if (!ptr) return;
    
    // Check if pointer belongs to our pool
    bool in_pool = false;
    for (const auto& chunk : chunks_) {
      auto base = reinterpret_cast<T*>(chunk->storage.get());
      if (ptr >= base && ptr < base + chunk->capacity) {
        in_pool = true;
        break;
      }
    }
    
    if (in_pool) {
      // Add to free list for reuse
      free_list_.push_back(ptr);
      total_deallocated_++;
    } else if (config_.allow_fallback) {
      // Was allocated with fallback, use standard delete
      ::operator delete(ptr);
    }
  }
  
  /// \brief Destroy object and return memory to pool
  void destroy(T* ptr) noexcept {
    if (ptr) {
      ptr->~T();
      deallocate(ptr);
    }
  }
  
  /// \brief Reset pool, invalidating all allocations
  /// \warning All allocated objects become invalid after reset
  void reset() noexcept {
    for (auto& chunk : chunks_) {
      chunk->reset();
    }
    free_list_.clear();
    total_allocated_ = 0;
    total_deallocated_ = 0;
    current_chunk_index_ = 0;
  }
  
  /// \brief Get pool statistics
  struct pool_stats {
    std::size_t chunks_count;
    std::size_t total_capacity;
    std::size_t allocated_count;
    std::size_t free_list_size;
    std::size_t bytes_allocated;
    
    double utilization() const noexcept {
      return total_capacity > 0 ? 
             static_cast<double>(allocated_count - free_list_size) / total_capacity : 
             0.0;
    }
  };
  
  pool_stats get_stats() const noexcept {
    pool_stats stats{};
    stats.chunks_count = chunks_.size();
    stats.allocated_count = total_allocated_ - total_deallocated_;
    stats.free_list_size = free_list_.size();
    
    for (const auto& chunk : chunks_) {
      stats.total_capacity += chunk->capacity;
    }
    
    stats.bytes_allocated = stats.total_capacity * sizeof(T);
    return stats;
  }
  
private:
  /// \brief Grow pool by adding a new chunk
  bool grow_pool() {
    if (chunks_.empty()) {
      // First allocation
      std::size_t initial_size = config_.initial_capacity > 0 ? 
                                 config_.initial_capacity : 1024;
      chunks_.emplace_back(std::make_unique<memory_chunk>(initial_size));
      current_chunk_index_ = 0;
      return true;
    }
    
    // Calculate new chunk size
    std::size_t last_size = chunks_.back()->capacity;
    std::size_t new_size = last_size * config_.growth_factor;
    
    // Check against maximum capacity
    if (config_.max_capacity > 0) {
      std::size_t current_total = 0;
      for (const auto& chunk : chunks_) {
        current_total += chunk->capacity;
      }
      
      if (current_total >= config_.max_capacity) {
        return false;  // Pool at maximum capacity
      }
      
      if (current_total + new_size > config_.max_capacity) {
        new_size = config_.max_capacity - current_total;
      }
    }
    
    if (new_size == 0) {
      return false;
    }
    
    try {
      chunks_.emplace_back(std::make_unique<memory_chunk>(new_size));
      current_chunk_index_ = chunks_.size() - 1;
      return true;
    } catch (const std::bad_alloc&) {
      return false;
    }
  }
};

/// \brief Scoped pool allocator for automatic cleanup
template <typename T>
class scoped_pool_allocator {
  memory_pool<T>& pool_;
  std::vector<T*> allocations_;
  
public:
  explicit scoped_pool_allocator(memory_pool<T>& pool) 
    : pool_(pool) {
    allocations_.reserve(1000);  // Pre-size for efficiency
  }
  
  ~scoped_pool_allocator() {
    // Destroy all allocations in reverse order
    for (auto it = allocations_.rbegin(); it != allocations_.rend(); ++it) {
      pool_.destroy(*it);
    }
  }
  
  // Disable copy
  scoped_pool_allocator(const scoped_pool_allocator&) = delete;
  scoped_pool_allocator& operator=(const scoped_pool_allocator&) = delete;
  
  template <typename... Args>
  T* construct(Args&&... args) {
    T* ptr = pool_.construct(std::forward<Args>(args)...);
    if (ptr) {
      allocations_.push_back(ptr);
    }
    return ptr;
  }
  
  void destroy(T* ptr) {
    auto it = std::find(allocations_.begin(), allocations_.end(), ptr);
    if (it != allocations_.end()) {
      pool_.destroy(ptr);
      allocations_.erase(it);
    }
  }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_MEMORY_POOL_HPP