// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_SOBOL_OWEN_HPP
#define BOOST_MATH_CUBATURE_DETAIL_SOBOL_OWEN_HPP

// STL headers first per Boost conventions
#include <cstdint>
#include <random>

namespace boost { namespace math { namespace cubature { namespace detail {

/// \file sobol_owen.hpp
/// \brief Owen scrambling for randomized quasi-Monte Carlo
/// \details Implements digital scrambling for Sobol sequences
///          Following ALGORITHMS.md section 3.2

/// \brief Owen digital scrambling for a single value
/// \param u Value in [0,1) to scramble
/// \param seed Scrambling seed (deterministic)
/// \return Scrambled value in [0,1)
template <typename Real>
Real owen_scramble(Real u, std::uint32_t seed) {
    static constexpr int max_digits = 32;
    
    // Convert to integer representation for bit manipulation
    const std::uint64_t max_val = 1ULL << max_digits;
    std::uint32_t u_int = static_cast<std::uint32_t>(u * max_val);
    
    // Extract binary digits
    std::uint32_t digits[max_digits] = {0};
    for (int i = 0; i < max_digits; ++i) {
        digits[i] = (u_int >> (max_digits - 1 - i)) & 1;
    }
    
    // Initialize RNG with seed
    std::mt19937 rng(seed);
    std::uniform_int_distribution<std::uint32_t> dist(0, 1);
    
    // Apply nested uniform scrambling
    std::uint32_t scrambled_digits[max_digits] = {0};
    
    for (int i = 0; i < max_digits; ++i) {
        // Build context from parent digits
        std::uint32_t context = 0;
        for (int j = 0; j < i; ++j) {
            context = (context << 1) | scrambled_digits[j];
        }
        
        // Generate deterministic permutation based on seed and context
        // This ensures each digit's scrambling depends on all parent digits
        std::uint32_t perm_seed = seed + i * 65521 + context * 32749;
        rng.seed(perm_seed);
        std::uint32_t perm = dist(rng);
        
        // Apply permutation (for base-2, this is XOR)
        scrambled_digits[i] = digits[i] ^ perm;
    }
    
    // Reconstruct scrambled value
    std::uint32_t scrambled_int = 0;
    for (int i = 0; i < max_digits; ++i) {
        scrambled_int = (scrambled_int << 1) | scrambled_digits[i];
    }
    
    return Real(scrambled_int) / Real(max_val);
}

/// \brief Full nested uniform Owen scrambling (production version)
/// \details Implements the proper nested uniform scrambling algorithm
///          where each digit is scrambled based on parent digits
template <typename Real>
class owen_scrambler {
private:
    std::mt19937 rng_;
    static constexpr int base_ = 2;  // Binary scrambling for Sobol
    static constexpr int max_digits_ = 32;
    
public:
    explicit owen_scrambler(std::uint32_t seed) : rng_(seed) {}
    
    Real scramble(Real u) {
        // Extract binary digits
        std::uint32_t digits[max_digits_] = {0};
        std::uint32_t u_int = static_cast<std::uint32_t>(u * (1ULL << max_digits_));
        
        for (int i = 0; i < max_digits_; ++i) {
            digits[i] = (u_int >> (max_digits_ - 1 - i)) & 1;
        }
        
        // Apply nested uniform scrambling
        std::uniform_int_distribution<std::uint32_t> dist(0, 1);
        std::uint32_t scrambled_digits[max_digits_] = {0};
        
        for (int i = 0; i < max_digits_; ++i) {
            // Scrambling depends on parent digits
            std::uint32_t context = 0;
            for (int j = 0; j < i; ++j) {
                context = (context << 1) | scrambled_digits[j];
            }
            
            // Generate scrambling permutation based on context
            rng_.seed(rng_() + context);
            std::uint32_t perm = dist(rng_);
            
            // Apply permutation (for base-2, this is just XOR)
            scrambled_digits[i] = digits[i] ^ perm;
        }
        
        // Reconstruct scrambled value
        std::uint32_t scrambled_int = 0;
        for (int i = 0; i < max_digits_; ++i) {
            scrambled_int = (scrambled_int << 1) | scrambled_digits[i];
        }
        
        return Real(scrambled_int) / Real(1ULL << max_digits_);
    }
};

}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_SOBOL_OWEN_HPP