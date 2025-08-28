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
/// \brief Concepts/traits for scalar and vector-valued integrands.

#if __cpp_concepts

// Helper concepts for multiple accepted call signatures
namespace detail {
#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  template <class F, class Real>
  concept accepts_span_scalar = requires(F f, std::span<const Real> x) {
    { f(x) } -> std::convertible_to<Real>;
  };

  template <class F, class Real>
  concept accepts_span_vector = requires(F f, std::span<const Real> x, Real* out, std::size_t m) {
    { f(x, out, m) } -> std::same_as<void>;
  };
#endif

  template <class F, class Real>
  concept accepts_vec_scalar = requires(F f, const std::vector<Real>& x) {
    { f(x) } -> std::convertible_to<Real>;
  };

  template <class F, class Real>
  concept accepts_vec_vector = requires(F f, const std::vector<Real>& x, Real* out, std::size_t m) {
    { f(x, out, m) } -> std::same_as<void>;
  };

  template <class F, class Real>
  concept accepts_ptr_scalar = requires(F f, const Real* p, std::size_t n) {
    { f(p, n) } -> std::convertible_to<Real>;
  };

  template <class F, class Real>
  concept accepts_ptr_vector = requires(F f, const Real* p, std::size_t n, Real* out, std::size_t m) {
    { f(p, n, out, m) } -> std::same_as<void>;
  };
}

template <class F, class Real>
concept ScalarIntegrand =
#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  detail::accepts_span_scalar<F, Real> ||
#endif
  detail::accepts_vec_scalar<F, Real> ||
  detail::accepts_ptr_scalar<F, Real>;

template <class F, class Real>
concept VectorIntegrand =
#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  detail::accepts_span_vector<F, Real> ||
#endif
  detail::accepts_vec_vector<F, Real> ||
  detail::accepts_ptr_vector<F, Real>;

#else // __cpp_concepts

// C++17 fallback: detection traits for the same accepted signatures
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
  struct has_ptr_vector<F, Real, decltype(void(std::declval<F&>()((const Real*)nullptr, std::declval<std::size_t>(), (Real*)nullptr, std::declval<std::size_t>())))> : std::true_type {};

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

template <class F, class Real>
struct is_scalar_integrand : std::integral_constant<bool,
#ifdef BOOST_MATH_CUBATURE_HAVE_SPAN
  detail::has_span_scalar<F, Real>::value ||
#endif
  detail::has_vec_scalar<F, Real>::value ||
  detail::has_ptr_scalar<F, Real>::value> {};

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

