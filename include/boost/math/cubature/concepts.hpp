// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_CONCEPTS_HPP
#define BOOST_MATH_CUBATURE_CONCEPTS_HPP

#include <type_traits>
#include <cstddef>
#include <vector>

#if __cplusplus >= 202002L
  #include <span>
  #ifdef __cpp_lib_span
    #define BOOST_MATH_CUBATURE_HAVE_SPAN 1
  #endif
#elif defined(__has_include)
  #if __has_include(<span>)
    #include <span>
    #ifdef __cpp_lib_span
      #define BOOST_MATH_CUBATURE_HAVE_SPAN 1
    #endif
  #endif
#endif

namespace boost { namespace math { namespace cubature {

/// \file concepts.hpp
/// \brief Concepts and traits for scalar and vector-valued integrands
/// \details Provides C++20 concepts and C++17 fallback traits for detecting
///          valid integrand signatures. Supports multiple calling conventions
///          including raw pointers, vectors, and spans.
/// 
/// \section concepts_usage Usage
/// \code
/// // Scalar integrand examples
/// auto f1 = [](const Real* x, std::size_t n) { return x[0] * x[0]; };
/// auto f2 = [](const std::vector<Real>& x) { return sin(x[0]); };
/// auto f3 = [](std::span<const Real> x) { return exp(x[0]); };  // C++20
/// 
/// // Vector integrand examples  
/// auto g1 = [](const Real* x, Real* out, std::size_t m) {
///     out[0] = x[0]; out[1] = x[0]*x[0]; 
/// };
/// auto g2 = [](const std::vector<Real>& x, Real* out, std::size_t m) {
///     for(std::size_t i = 0; i < m; ++i) out[i] = pow(x[0], i);
/// };
/// \endcode
/// 
/// \section concepts_detection Detection
/// \code
/// // C++20 with concepts
/// if constexpr (ScalarIntegrand<F, Real>) { /* scalar path */ }
/// if constexpr (VectorIntegrand<F, Real>) { /* vector path */ }
/// 
/// // C++17 with traits
/// if constexpr (is_scalar_integrand<F, Real>::value) { /* scalar */ }
/// if constexpr (is_vector_integrand<F, Real>::value) { /* vector */ }
/// \endcode

#if __cpp_concepts

/// \brief Helper concepts for multiple accepted call signatures
/// \details Internal concepts that check individual signatures
namespace detail {
#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  /// \brief Check if F accepts span<const Real> and returns Real
  template <class F, class Real>
  concept accepts_span_scalar = requires(F f, std::span<const Real> x) {
    { f(x) } -> std::convertible_to<Real>;
  };

  /// \brief Check if F accepts span for vector output
  template <class F, class Real>
  concept accepts_span_vector = requires(F f, std::span<const Real> x, Real* out, std::size_t m) {
    { f(x, out, m) } -> std::same_as<void>;
  };
#endif

  /// \brief Check if F accepts vector<Real> and returns Real
  template <class F, class Real>
  concept accepts_vec_scalar = requires(F f, const std::vector<Real>& x) {
    { f(x) } -> std::convertible_to<Real>;
  };

  /// \brief Check if F accepts vector for vector output
  template <class F, class Real>
  concept accepts_vec_vector = requires(F f, const std::vector<Real>& x, Real* out, std::size_t m) {
    { f(x, out, m) } -> std::same_as<void>;
  };

  /// \brief Check if F accepts (const Real*, size_t) and returns Real
  template <class F, class Real>
  concept accepts_ptr_scalar = requires(F f, const Real* p, std::size_t n) {
    { f(p, n) } -> std::convertible_to<Real>;
  };

  /// \brief Check if F accepts raw pointers for vector output
  template <class F, class Real>
  concept accepts_ptr_vector = requires(F f, const Real* p, Real* out, std::size_t m) {
    { f(p, out, m) } -> std::same_as<void>;
  };
}

/// \brief Concept for scalar-valued integrands
/// \tparam F Integrand type
/// \tparam Real Floating point type
/// \details Accepts any of: f(span), f(vector), f(ptr, size)
template <class F, class Real>
concept ScalarIntegrand =
#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  detail::accepts_span_scalar<F, Real> ||
#endif
  detail::accepts_vec_scalar<F, Real> ||
  detail::accepts_ptr_scalar<F, Real>;

/// \brief Concept for vector-valued integrands
/// \tparam F Integrand type
/// \tparam Real Floating point type
/// \details Accepts any of: f(span, out, m), f(vector, out, m), f(ptr, out, m)
template <class F, class Real>
concept VectorIntegrand =
#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  detail::accepts_span_vector<F, Real> ||
#endif
  detail::accepts_vec_vector<F, Real> ||
  detail::accepts_ptr_vector<F, Real>;

#else // __cpp_concepts

/// \brief C++17 fallback: detection traits for the same accepted signatures
/// \details SFINAE-based trait detection for pre-C++20 compilers
namespace detail {
  template <class, class, class = void>
  struct has_vec_scalar : std::false_type {};
  template <class F, class Real>
  struct has_vec_scalar<F, Real, decltype(void(std::declval<F&>()(std::declval<const std::vector<Real>&>())))> : std::true_type {};

  template <class, class, class = void>
  struct has_vec_vector : std::false_type {};
  template <class F, class Real>
  struct has_vec_vector<F, Real, decltype(void(std::declval<F&>()(std::declval<const std::vector<Real>&>(), (Real*)nullptr, std::declval<std::size_t>())))> : std::true_type {};

  template <class, class, class = void>
  struct has_ptr_scalar : std::false_type {};
  template <class F, class Real>
  struct has_ptr_scalar<F, Real, decltype(void(std::declval<F&>()((const Real*)nullptr, std::declval<std::size_t>())))> : std::true_type {};

  template <class, class, class = void>
  struct has_ptr_vector : std::false_type {};
  template <class F, class Real>
  struct has_ptr_vector<F, Real, decltype(void(std::declval<F&>()((const Real*)nullptr, (Real*)nullptr, std::declval<std::size_t>())))> : std::true_type {};

#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  template <class, class, class = void>
  struct has_span_scalar : std::false_type {};
  template <class F, class Real>
  struct has_span_scalar<F, Real, decltype(void(std::declval<F&>()(std::declval<std::span<const Real>>())))> : std::true_type {};

  template <class, class, class = void>
  struct has_span_vector : std::false_type {};
  template <class F, class Real>
  struct has_span_vector<F, Real, decltype(void(std::declval<F&>()(std::declval<std::span<const Real>>(), (Real*)nullptr, std::declval<std::size_t>())))> : std::true_type {};
#endif
}

/// \brief C++17 trait for scalar integrand detection
/// \tparam F Integrand type
/// \tparam Real Floating point type
template <class F, class Real>
struct is_scalar_integrand : std::integral_constant<bool,
#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  detail::has_span_scalar<F, Real>::value ||
#endif
  detail::has_vec_scalar<F, Real>::value ||
  detail::has_ptr_scalar<F, Real>::value> {};

/// \brief C++17 trait for vector integrand detection
/// \tparam F Integrand type
/// \tparam Real Floating point type
template <class F, class Real>
struct is_vector_integrand : std::integral_constant<bool,
#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  detail::has_span_vector<F, Real>::value ||
#endif
  detail::has_vec_vector<F, Real>::value ||
  detail::has_ptr_vector<F, Real>::value> {};

#endif // __cpp_concepts

}}} // namespace boost::math::cubature

#endif // BOOST_MATH_CUBATURE_CONCEPTS_HPP

