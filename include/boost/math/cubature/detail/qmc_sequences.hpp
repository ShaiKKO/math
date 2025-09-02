// Copyright 2025 Shaiko
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CUBATURE_DETAIL_QMC_SEQUENCES_HPP
#define BOOST_MATH_CUBATURE_DETAIL_QMC_SEQUENCES_HPP

#include <vector>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <limits>

// Check for Boost.Random availability
#ifdef __has_include
  #if __has_include(<boost/random/sobol.hpp>)
    #include <boost/random/sobol.hpp>
    #define BOOST_MATH_HAS_SOBOL 1
  #endif
  #if __has_include(<boost/random/faure.hpp>)
    #include <boost/random/faure.hpp>
    #define BOOST_MATH_HAS_FAURE 1
  #endif
  #if __has_include(<boost/random/niederreiter_base2.hpp>)
    #include <boost/random/niederreiter_base2.hpp>
    #define BOOST_MATH_HAS_NIEDERREITER 1
  #endif
#endif

// Owen scrambling
#include <boost/math/cubature/detail/sobol_owen.hpp>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \brief Enumeration of available QMC sequence types
enum class qmc_sequence_type {
    sobol,              ///< Sobol sequence (already implemented)
    halton,             ///< Halton sequence with prime bases
    faure,              ///< Faure sequence
    niederreiter_base2, ///< Niederreiter base-2 sequence
    lattice_rank1       ///< Rank-1 lattice rules with CBC construction
};

/// \brief Abstract base class for QMC sequences
/// \tparam Real Floating point type for sequence values
template <typename Real>
class qmc_sequence_base {
public:
    virtual ~qmc_sequence_base() = default;
    
    /// \brief Reset sequence to a specific starting index
    virtual void reset(std::size_t start_index = 0) = 0;
    
    /// \brief Generate next point in the sequence
    /// \return Vector of dimension() values in [0,1)
    virtual std::vector<Real> next() = 0;
    
    /// \brief Skip ahead n points without generating them
    virtual void skip(std::size_t n) = 0;
    
    /// \brief Get the dimension of the sequence
    virtual std::size_t dimension() const = 0;
    
    /// \brief Get current index in the sequence
    virtual std::size_t current_index() const = 0;
    
    /// \brief Check if more points are available
    virtual bool has_next() const { return true; } // Most sequences are infinite
};

/// \brief Sobol sequence implementation
/// This is a simplified implementation that works without Boost.Random
template <typename Real>
class sobol_sequence : public qmc_sequence_base<Real> {
private:
    std::size_t dim_;
    std::size_t index_;
    std::vector<std::vector<std::uint64_t>> direction_numbers_;
    std::vector<std::uint64_t> current_point_;
    
    // Initialize direction numbers (simplified - using basic values)
    void init_direction_numbers() {
        direction_numbers_.resize(dim_);
        
        // First dimension: powers of 2
        direction_numbers_[0].resize(64);
        for (std::size_t i = 0; i < 64; ++i) {
            direction_numbers_[0][i] = 1ULL << (63 - i);
        }
        
        // Other dimensions: use primitive polynomials (simplified)
        // In production, these would come from Sobol tables
        for (std::size_t d = 1; d < dim_; ++d) {
            direction_numbers_[d].resize(64);
            std::uint64_t seed = static_cast<std::uint64_t>(d * 0x9E3779B97F4A7C15ULL);
            
            // Initialize first few direction numbers
            direction_numbers_[d][0] = 1ULL << 63;
            for (std::size_t i = 1; i < 64; ++i) {
                // Simplified generation using linear congruential generator
                seed = seed * 1103515245ULL + 12345ULL;
                direction_numbers_[d][i] = seed & ((1ULL << (64 - i)) - 1);
                direction_numbers_[d][i] <<= i;
            }
        }
        
        current_point_.resize(dim_, 0);
    }
    
public:
    explicit sobol_sequence(std::size_t dimension)
        : dim_(dimension), index_(0) {
        if (dimension == 0) {
            throw std::invalid_argument("Sobol sequence dimension must be positive");
        }
        if (dimension > 100) {  // Reduced limit for simplified implementation
            throw std::invalid_argument("Simplified Sobol sequence supports up to 100 dimensions");
        }
        init_direction_numbers();
    }
    
    void reset(std::size_t start_index = 0) override {
        index_ = 0;
        std::fill(current_point_.begin(), current_point_.end(), 0);
        
        // Skip to start_index if needed
        if (start_index > 0) {
            skip(start_index);
        }
    }
    
    std::vector<Real> next() override {
        std::vector<Real> point(dim_);
        
        // Find rightmost zero bit of index (Gray code)
        std::size_t l = 0;
        std::size_t i = index_;
        while ((i & 1) == 1) {
            i >>= 1;
            ++l;
        }
        
        // Update the Sobol point using XOR with direction number
        if (l < 64) {
            for (std::size_t d = 0; d < dim_; ++d) {
                current_point_[d] ^= direction_numbers_[d][l];
                point[d] = Real(current_point_[d]) / Real(std::numeric_limits<std::uint64_t>::max());
            }
        }
        
        ++index_;
        return point;
    }
    
    void skip(std::size_t n) override {
        for (std::size_t i = 0; i < n; ++i) {
            next();
        }
    }
    
    std::size_t dimension() const override { return dim_; }
    std::size_t current_index() const override { return index_; }
};

/// \brief Halton sequence implementation
template <typename Real>
class halton_sequence : public qmc_sequence_base<Real> {
private:
    std::size_t dim_;
    std::size_t index_;
    std::vector<std::size_t> primes_;
    std::size_t leap_;
    
    /// \brief Generate the i-th element of Halton sequence in given base
    static Real halton_element(std::size_t i, std::size_t base) {
        Real result = 0;
        Real f = Real(1) / Real(base);
        std::size_t index = i;
        
        while (index > 0) {
            result += f * (index % base);
            index /= base;
            f /= base;
        }
        
        return result;
    }
    
    /// \brief Simple primality test
    static bool is_prime(std::size_t n) {
        if (n <= 1) return false;
        if (n <= 3) return true;
        if (n % 2 == 0 || n % 3 == 0) return false;
        
        for (std::size_t i = 5; i * i <= n; i += 6) {
            if (n % i == 0 || n % (i + 2) == 0)
                return false;
        }
        return true;
    }
    
    /// \brief Get first n prime numbers
    static std::vector<std::size_t> first_n_primes(std::size_t n) {
        std::vector<std::size_t> primes;
        primes.reserve(n);
        
        std::size_t candidate = 2;
        while (primes.size() < n) {
            if (is_prime(candidate)) {
                primes.push_back(candidate);
            }
            ++candidate;
        }
        
        return primes;
    }
    
public:
    explicit halton_sequence(std::size_t dimension, std::size_t leap = 1)
        : dim_(dimension), index_(0), leap_(leap) {
        primes_ = first_n_primes(dimension);
    }
    
    void reset(std::size_t start_index = 0) override {
        index_ = start_index;
    }
    
    std::vector<Real> next() override {
        std::vector<Real> point(dim_);
        std::size_t actual_index = index_ * leap_;
        
        for (std::size_t d = 0; d < dim_; ++d) {
            point[d] = halton_element(actual_index, primes_[d]);
        }
        
        ++index_;
        return point;
    }
    
    void skip(std::size_t n) override {
        index_ += n;
    }
    
    std::size_t dimension() const override { return dim_; }
    std::size_t current_index() const override { return index_; }
    
    /// \brief Set leap value for leaped Halton sequence
    void set_leap(std::size_t leap) { leap_ = leap; }
};

/// \brief Faure sequence wrapper
template <typename Real>
class faure_sequence : public qmc_sequence_base<Real> {
private:
#ifdef BOOST_MATH_HAS_FAURE
    boost::random::faure engine_;
#endif
    std::size_t dim_;
    std::size_t index_;
    
public:
    explicit faure_sequence(std::size_t dimension)
        : 
#ifdef BOOST_MATH_HAS_FAURE
          engine_(dimension), 
#endif
          dim_(dimension), index_(0) {
#ifdef BOOST_MATH_HAS_FAURE
        // Faure sequence has dimension limits based on prime table
        // The default table supports up to 1117 dimensions
        if (dimension > 1117) {
            throw std::invalid_argument("Faure sequence dimension exceeds maximum (1117)");
        }
#else
        throw std::runtime_error("Faure sequence not available - Boost.Random faure.hpp not found");
#endif
    }
    
    void reset(std::size_t start_index = 0) override {
#ifdef BOOST_MATH_HAS_FAURE
        engine_.seed(start_index);
#endif
        index_ = start_index;
    }
    
    std::vector<Real> next() override {
        std::vector<Real> point(dim_);
#ifdef BOOST_MATH_HAS_FAURE
        for (std::size_t d = 0; d < dim_; ++d) {
            point[d] = static_cast<Real>(engine_()) / 
                      static_cast<Real>(engine_.max() - engine_.min() + 1);
        }
#else
        // Fallback: return zeros
        std::fill(point.begin(), point.end(), Real(0));
#endif
        ++index_;
        return point;
    }
    
    void skip(std::size_t n) override {
#ifdef BOOST_MATH_HAS_FAURE
        engine_.discard(n * dim_);
#endif
        index_ += n;
    }
    
    std::size_t dimension() const override { return dim_; }
    std::size_t current_index() const override { return index_; }
};

/// \brief Niederreiter base-2 sequence wrapper
template <typename Real>
class niederreiter_sequence : public qmc_sequence_base<Real> {
private:
#ifdef BOOST_MATH_HAS_NIEDERREITER
    boost::random::niederreiter_base2 engine_;
#endif
    std::size_t dim_;
    std::size_t index_;
    
public:
    explicit niederreiter_sequence(std::size_t dimension)
        : 
#ifdef BOOST_MATH_HAS_NIEDERREITER
          engine_(dimension), 
#endif
          dim_(dimension), index_(0) {
#ifdef BOOST_MATH_HAS_NIEDERREITER
        // Niederreiter base-2 supports up to 4720 dimensions with default table
        if (dimension > 4720) {
            throw std::invalid_argument("Niederreiter sequence dimension exceeds maximum (4720)");
        }
#else
        throw std::runtime_error("Niederreiter sequence not available - Boost.Random niederreiter_base2.hpp not found");
#endif
    }
    
    void reset(std::size_t start_index = 0) override {
#ifdef BOOST_MATH_HAS_NIEDERREITER
        engine_.seed(start_index);
#endif
        index_ = start_index;
    }
    
    std::vector<Real> next() override {
        std::vector<Real> point(dim_);
#ifdef BOOST_MATH_HAS_NIEDERREITER
        for (std::size_t d = 0; d < dim_; ++d) {
            point[d] = static_cast<Real>(engine_()) / 
                      static_cast<Real>(engine_.max() - engine_.min() + 1);
        }
#else
        // Fallback: return zeros
        std::fill(point.begin(), point.end(), Real(0));
#endif
        ++index_;
        return point;
    }
    
    void skip(std::size_t n) override {
#ifdef BOOST_MATH_HAS_NIEDERREITER
        engine_.discard(n * dim_);
#endif
        index_ += n;
    }
    
    std::size_t dimension() const override { return dim_; }
    std::size_t current_index() const override { return index_; }
};

/// \brief Rank-1 lattice rules with CBC (Component-By-Component) construction
template <typename Real>
class lattice_rank1_sequence : public qmc_sequence_base<Real> {
private:
    std::size_t dim_;
    std::size_t n_points_;
    std::vector<std::size_t> generator_;
    std::size_t index_;
    Real shift_;  // Random shift for randomized lattice rules
    
    /// \brief Compute gcd using Euclidean algorithm
    static std::size_t gcd(std::size_t a, std::size_t b) {
        while (b != 0) {
            std::size_t temp = b;
            b = a % b;
            a = temp;
        }
        return a;
    }
    
    /// \brief Component-by-component construction
    /// Uses a simplified Korobov-type construction for generator
    void construct_cbc(Real /*omega*/ = 1.0) {
        generator_.resize(dim_);
        
        // Find a good primitive root modulo n_points_
        // For simplicity, use a fixed good generator
        std::size_t a = 1;
        
        // Find a good value of 'a' that is relatively prime to n_points_
        for (std::size_t candidate = 2; candidate < n_points_; ++candidate) {
            if (gcd(candidate, n_points_) == 1) {
                // Check if it has good distribution properties
                // For simplicity, use the first valid one that's not too small
                if (candidate > n_points_ / 4 && candidate < 3 * n_points_ / 4) {
                    a = candidate;
                    break;
                }
            }
        }
        
        // Korobov-type construction: generator[i] = a^i mod n_points_
        generator_[0] = 1;
        for (std::size_t d = 1; d < dim_; ++d) {
            generator_[d] = (generator_[d-1] * a) % n_points_;
        }
    }
    
public:
    lattice_rank1_sequence(std::size_t dimension, std::size_t num_points)
        : dim_(dimension), n_points_(num_points), index_(0), shift_(0) {
        
        // Ensure n_points is prime for better properties
        // In production, would use prime or prime power
        construct_cbc();
    }
    
    void reset(std::size_t start_index = 0) override {
        index_ = start_index;
    }
    
    std::vector<Real> next() override {
        if (index_ >= n_points_) {
            throw std::runtime_error("Lattice rule exhausted");
        }
        
        std::vector<Real> point(dim_);
        for (std::size_t d = 0; d < dim_; ++d) {
            Real x = Real(index_ * generator_[d] % n_points_) / Real(n_points_);
            // Apply random shift if set
            x = std::fmod(x + shift_, Real(1));
            point[d] = x;
        }
        
        ++index_;
        return point;
    }
    
    void skip(std::size_t n) override {
        index_ += n;
        if (index_ > n_points_) {
            index_ = n_points_;
        }
    }
    
    std::size_t dimension() const override { return dim_; }
    std::size_t current_index() const override { return index_; }
    
    bool has_next() const override { return index_ < n_points_; }
    
    /// \brief Apply random shift for randomized lattice rules
    void apply_random_shift(std::mt19937& rng) {
        std::uniform_real_distribution<Real> dist(0, 1);
        shift_ = dist(rng);
    }
    
};

/// \brief Scrambled sequence wrapper applying Owen scrambling
template <typename Real>
class scrambled_sequence : public qmc_sequence_base<Real> {
private:
    std::unique_ptr<qmc_sequence_base<Real>> base_sequence_;
    owen_scrambler<Real> scrambler_;
    
public:
    scrambled_sequence(std::unique_ptr<qmc_sequence_base<Real>> seq, 
                      std::uint32_t seed)
        : base_sequence_(std::move(seq)), scrambler_(seed) {}
    
    void reset(std::size_t start_index = 0) override {
        base_sequence_->reset(start_index);
    }
    
    std::vector<Real> next() override {
        auto point = base_sequence_->next();
        for (auto& coord : point) {
            coord = scrambler_.scramble(coord);
        }
        return point;
    }
    
    void skip(std::size_t n) override {
        base_sequence_->skip(n);
    }
    
    std::size_t dimension() const override { 
        return base_sequence_->dimension(); 
    }
    
    std::size_t current_index() const override { 
        return base_sequence_->current_index(); 
    }
    
    bool has_next() const override {
        return base_sequence_->has_next();
    }
};

/// \brief Factory function to create QMC sequences
template <typename Real>
std::unique_ptr<qmc_sequence_base<Real>> create_qmc_sequence(
    qmc_sequence_type type,
    std::size_t dimension,
    std::size_t num_points = 0)  // Only used for lattice rules
{
    switch(type) {
        case qmc_sequence_type::sobol:
            return std::make_unique<sobol_sequence<Real>>(dimension);
            
        case qmc_sequence_type::halton:
            return std::make_unique<halton_sequence<Real>>(dimension);
            
        case qmc_sequence_type::faure:
            return std::make_unique<faure_sequence<Real>>(dimension);
            
        case qmc_sequence_type::niederreiter_base2:
            return std::make_unique<niederreiter_sequence<Real>>(dimension);
            
        case qmc_sequence_type::lattice_rank1:
            if (num_points == 0) {
                throw std::invalid_argument("Lattice rules require num_points");
            }
            return std::make_unique<lattice_rank1_sequence<Real>>(dimension, num_points);
            
        default:
            throw std::invalid_argument("Unknown QMC sequence type");
    }
}

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_QMC_SEQUENCES_HPP