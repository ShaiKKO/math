// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_RANDOM_SOBOL_HPP
#define BOOST_RANDOM_SOBOL_HPP

/// \file sobol.hpp
/// \brief Sobol low-discrepancy sequence generator
/// \details Implements the Sobol sequence, a quasi-random low-discrepancy
///          sequence for quasi-Monte Carlo integration. Based on the
///          algorithms from Bratley & Fox (1988) and Joe & Kuo (2008).
///
/// \references
/// - I.M. Sobol, "On the distribution of points in a cube and the approximate
///   evaluation of integrals", USSR Comp. Math. Math. Phys. 7 (1967), 86-112
/// - P. Bratley and B.L. Fox, "Algorithm 659: Implementing Sobol's quasirandom
///   sequence generator", ACM Trans. Math. Software 14 (1988), 88-100
/// - S. Joe and F.Y. Kuo, "Constructing Sobol sequences with better two-dimensional
///   projections", SIAM J. Sci. Comput. 30 (2008), 2635-2654

#include <vector>
#include <array>
#include <cstdint>
#include <limits>
#include <algorithm>

namespace boost { namespace random {

/// \brief Direction numbers for Sobol sequence (up to dimension 21201)
/// \details These are the primitive polynomials and initial direction
///          numbers from Joe & Kuo (2008) for dimensions up to 21201.
///          For this implementation, we include only dimensions up to 20.
namespace detail {

struct sobol_direction_numbers {
    static constexpr std::size_t max_dim = 20;
    static constexpr std::size_t max_bits = 32;
    
    /// Primitive polynomials for each dimension (in binary representation)
    static constexpr std::uint32_t primitive_polynomials[max_dim] = {
        1,      // dim 1: x (degree 1)
        3,      // dim 2: x + 1 (degree 1)
        7,      // dim 3: x^2 + x + 1 (degree 2)
        11,     // dim 4: x^3 + x + 1 (degree 3)
        13,     // dim 5: x^3 + x^2 + 1 (degree 3)
        19,     // dim 6: x^4 + x + 1 (degree 4)
        25,     // dim 7: x^4 + x^3 + 1 (degree 4)
        37,     // dim 8: x^5 + x^2 + 1 (degree 5)
        59,     // dim 9: x^5 + x^4 + x^3 + x + 1 (degree 5)
        47,     // dim 10: x^5 + x^3 + x^2 + x + 1 (degree 5)
        61,     // dim 11: x^5 + x^4 + x^3 + x^2 + 1 (degree 5)
        55,     // dim 12: x^5 + x^4 + x^2 + x + 1 (degree 5)
        41,     // dim 13: x^5 + x^3 + 1 (degree 5)
        67,     // dim 14: x^6 + x + 1 (degree 6)
        97,     // dim 15: x^6 + x^5 + 1 (degree 6)
        91,     // dim 16: x^6 + x^4 + x^3 + x + 1 (degree 6)
        109,    // dim 17: x^6 + x^5 + x^3 + x^2 + 1 (degree 6)
        103,    // dim 18: x^6 + x^5 + x^2 + x + 1 (degree 6)
        115,    // dim 19: x^6 + x^5 + x^4 + x + 1 (degree 6)
        131     // dim 20: x^7 + x + 1 (degree 7)
    };
    
    /// Degree of primitive polynomial for each dimension
    static constexpr std::uint32_t polynomial_degrees[max_dim] = {
        1, 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7
    };
    
    /// Initial direction numbers m_i for each dimension
    /// First index is dimension-1, second index is i-1
    static constexpr std::uint32_t initial_direction_numbers[max_dim][8] = {
        {1, 0, 0, 0, 0, 0, 0, 0},  // dim 1
        {1, 0, 0, 0, 0, 0, 0, 0},  // dim 2
        {1, 1, 0, 0, 0, 0, 0, 0},  // dim 3
        {1, 3, 7, 0, 0, 0, 0, 0},  // dim 4
        {1, 1, 5, 0, 0, 0, 0, 0},  // dim 5
        {1, 3, 1, 1, 0, 0, 0, 0},  // dim 6
        {1, 1, 3, 7, 0, 0, 0, 0},  // dim 7
        {1, 3, 3, 9, 9, 0, 0, 0},  // dim 8
        {1, 3, 7, 13, 3, 0, 0, 0}, // dim 9
        {1, 1, 5, 11, 27, 0, 0, 0}, // dim 10
        {1, 3, 5, 1, 15, 0, 0, 0},  // dim 11
        {1, 1, 7, 3, 29, 0, 0, 0},  // dim 12
        {1, 3, 7, 7, 21, 0, 0, 0},  // dim 13
        {1, 1, 1, 9, 23, 37, 0, 0}, // dim 14
        {1, 3, 3, 5, 19, 33, 0, 0}, // dim 15
        {1, 1, 3, 13, 11, 7, 0, 0}, // dim 16
        {1, 1, 7, 13, 25, 5, 0, 0}, // dim 17
        {1, 3, 5, 11, 7, 11, 0, 0}, // dim 18
        {1, 1, 1, 3, 13, 39, 0, 0}, // dim 19
        {1, 3, 1, 15, 17, 63, 13, 0} // dim 20
    };
};

} // namespace detail

/// \brief Sobol quasi-random sequence generator
/// \tparam UIntType Unsigned integer type for output (typically uint32_t or uint64_t)
/// \tparam Dimension Number of dimensions (compile-time constant)
template <typename UIntType, std::size_t Dimension>
class sobol_engine {
public:
    static_assert(Dimension > 0 && Dimension <= detail::sobol_direction_numbers::max_dim,
                  "Dimension must be between 1 and 20");
    
    using result_type = std::vector<UIntType>;
    
private:
    static constexpr std::size_t bits = std::numeric_limits<UIntType>::digits;
    static constexpr UIntType norm = static_cast<UIntType>(1) << (bits - 1);
    
    std::uint32_t count_;  ///< Sequence counter
    std::array<std::uint32_t, Dimension> x_; ///< Current state vector
    std::array<std::array<std::uint32_t, 32>, Dimension> v_; ///< Direction numbers
    
    /// \brief Initialize direction numbers for all dimensions
    void initialize_direction_numbers() {
        // First dimension is special (van der Corput sequence)
        for (std::size_t i = 0; i < 32; ++i) {
            v_[0][i] = 1 << (31 - i);
        }
        
        // Initialize other dimensions using primitive polynomials
        for (std::size_t d = 1; d < Dimension; ++d) {
            std::uint32_t poly = detail::sobol_direction_numbers::primitive_polynomials[d];
            std::uint32_t degree = detail::sobol_direction_numbers::polynomial_degrees[d];
            
            // Use initial direction numbers
            for (std::size_t i = 0; i < degree; ++i) {
                v_[d][i] = detail::sobol_direction_numbers::initial_direction_numbers[d][i] 
                          << (31 - i);
            }
            
            // Generate remaining direction numbers using recurrence
            for (std::size_t i = degree; i < 32; ++i) {
                std::uint32_t newv = v_[d][i - degree];
                newv ^= (newv >> degree);
                
                for (std::size_t k = 1; k < degree; ++k) {
                    if ((poly >> (degree - k)) & 1) {
                        newv ^= v_[d][i - k];
                    }
                }
                
                v_[d][i] = newv;
            }
        }
    }
    
    /// \brief Find position of rightmost zero bit
    static std::uint32_t rightmost_zero_bit(std::uint32_t n) {
        std::uint32_t pos = 0;
        while (n & 1) {
            n >>= 1;
            ++pos;
        }
        return pos;
    }
    
public:
    /// \brief Constructor
    /// \param seed Starting index in sequence (default 0)
    explicit sobol_engine(std::size_t seed = 0) : count_(0) {
        initialize_direction_numbers();
        std::fill(x_.begin(), x_.end(), 0);
        
        // Skip to seed position
        for (std::size_t i = 0; i < seed; ++i) {
            (*this)();
        }
    }
    
    /// \brief Generate next point in Sobol sequence
    /// \return Vector of dimension values in [0, 2^bits - 1]
    result_type operator()() {
        if (count_ == 0) {
            // First point is always zero
            count_++;
            result_type result(Dimension, 0);
            return result;
        }
        
        // Find rightmost zero bit of count
        std::uint32_t c = rightmost_zero_bit(count_ - 1);
        
        // Update state using Gray code
        for (std::size_t d = 0; d < Dimension; ++d) {
            x_[d] ^= v_[d][c];
        }
        
        // Convert to output type
        result_type result(Dimension);
        for (std::size_t d = 0; d < Dimension; ++d) {
            if constexpr (sizeof(UIntType) == sizeof(std::uint32_t)) {
                result[d] = x_[d];
            } else {
                // Scale up to 64-bit if needed
                result[d] = static_cast<UIntType>(x_[d]) << 32;
            }
        }
        
        count_++;
        return result;
    }
    
    /// \brief Reset sequence to beginning
    void reset() {
        count_ = 0;
        std::fill(x_.begin(), x_.end(), 0);
    }
    
    /// \brief Discard n values
    void discard(std::size_t n) {
        for (std::size_t i = 0; i < n; ++i) {
            (*this)();
        }
    }
    
    /// \brief Get current sequence index
    std::uint32_t count() const { return count_; }
};

/// \brief Helper to generate normalized Sobol points in [0,1]^d
template <typename Real, std::size_t Dimension>
class sobol_generator {
private:
    sobol_engine<std::uint64_t, Dimension> engine_;
    static constexpr Real normalizer = Real(1) / (Real(1) << 32);
    
public:
    explicit sobol_generator(std::size_t seed = 0) : engine_(seed) {}
    
    /// \brief Generate next point in [0,1]^d
    std::vector<Real> operator()() {
        auto raw = engine_();
        std::vector<Real> point(Dimension);
        
        for (std::size_t d = 0; d < Dimension; ++d) {
            // Convert to [0,1] using upper 32 bits for better uniformity
            point[d] = static_cast<Real>(raw[d] >> 32) * normalizer;
        }
        
        return point;
    }
    
    void reset() { engine_.reset(); }
    void discard(std::size_t n) { engine_.discard(n); }
};

}} // namespace boost::random

#endif // BOOST_RANDOM_SOBOL_HPP