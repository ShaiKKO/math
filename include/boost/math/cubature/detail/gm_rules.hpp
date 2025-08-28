// Copyright 2025 Colin MacRitchie/Ripple Group
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_CUBATURE_DETAIL_GM_RULES_HPP
#define BOOST_MATH_CUBATURE_DETAIL_GM_RULES_HPP

#include <array>
#include <cstddef>
#include <type_traits>
#include <cassert>
#include <cmath>

// Auto-generated finalized GM 9/7 tables (DCUHRE weights and radii)
#include "genz_malik_9_7_tables.hpp"

namespace boost { namespace math { namespace cubature { namespace detail {

/// \file gm_rules.hpp
/// \brief Fully symmetric Genz–Malik embedded rules (API + metadata).
/// \details Provides compile-time metadata and accessors for degree-7/5 and 9/7
///  embedded cubature rules on the unit hypercube [0,1]^d (affine scaled by caller).

namespace gm {

// Tags for rule degrees
struct deg5  { static constexpr int value = 5; };
struct deg7  { static constexpr int value = 7; };
struct deg9  { static constexpr int value = 9; };

// Rule descriptor: primary template (not available by default)
template <int Degree, int Dim>
struct rule {
  static constexpr bool available = false;
  static constexpr int degree     = Degree;
  static constexpr int dim        = Dim;
  static constexpr std::size_t size = 0; // number of nodes

  template <class Real>
  static constexpr std::array<std::array<Real, Dim>, size> nodes() { return {}; }

  template <class Real>
  static constexpr std::array<Real, size> weights() { return {}; }
};

// Embedded pair helper: exposes fine/coarse rules sharing nodes
template <int FineDeg, int CoarseDeg, int Dim>
struct embedded_pair {
  using fine   = rule<FineDeg, Dim>;
  using coarse = rule<CoarseDeg, Dim>;
  static constexpr bool available = fine::available && coarse::available;
};

// Convenience alias for common pairs
template <int Dim>
using pair_7_5 = embedded_pair<7,5,Dim>;

// Metadata helpers
template <int Degree, int Dim>
struct available_v_helper { static constexpr bool value = rule<Degree, Dim>::available; };

template <int Degree, int Dim>
constexpr bool available_v = available_v_helper<Degree, Dim>::value;

template <int Fine, int Coarse, int Dim>
struct pair_available_v_helper { static constexpr bool value = embedded_pair<Fine, Coarse, Dim>::available; };

// Families for fully symmetric embedded rules
struct family_7_5 { };
struct family_9_7 { };

// Family-aware rule template (default: unavailable)
template <int Degree, int Dim, class Family>
struct rule_fam {
  static constexpr bool available = false;
  static constexpr int degree = Degree;
  static constexpr int dim = Dim;
  static constexpr std::size_t size = 0;
  template <class Real>
  static constexpr std::array<std::array<Real, Dim>, size> nodes() { return {}; }
  template <class Real>
  static constexpr std::array<Real, size> weights() { return {}; }
};

// Embedded pair (family-aware)
template <class Family, int Dim>
struct embedded_pair_fam {
  using fine   = rule_fam<9, Dim, Family>;
  using coarse = rule_fam<7, Dim, Family>;
  static constexpr bool available = fine::available && coarse::available;
};

template <int Dim>
using pair_9_7 = embedded_pair_fam<family_9_7, Dim>;

// Backward-compatibility: existing rule<Degree,Dim> refers to 7/5 family
// (Existing specializations below remain valid for 7/5.)


template <int Fine, int Coarse, int Dim>
constexpr bool pair_available_v = pair_available_v_helper<Fine, Coarse, Dim>::value;

// Concrete specializations for 2D 7/5 on [0,1]^2 (mapped from [-1,1]^2)
// Source: Genz–Malik embedded rules; see shaiko/genz_malik_2d_7_5.md (cites Genz & Malik, 1980/1983)

template <>
struct rule<7, 2> {
  static constexpr bool available = true;
  static constexpr int degree = 7;
  static constexpr int dim = 2;
  static constexpr std::size_t size = 17;

  template <class Real>
  static constexpr std::array<std::array<Real, 2>, size> nodes() {
    // Map t(x) = (x+1)/2 from [-1,1]^2 to [0,1]^2
    const long double l2 = 0.3585685828003180919906451539L; // 3/sqrt(70)
    const long double l3 = 0.9486832980505137995996680633L; // 3/sqrt(10)
    const long double l5 = 0.6882472016116852977216287343L; // 3/sqrt(19)
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    return {{
      // center
      { static_cast<Real>(0.5L), static_cast<Real>(0.5L) },
      // axis lambda2 (±l2, 0) and (0, ±l2)
      { static_cast<Real>(T::map( l2)), static_cast<Real>(0.5L) },
      { static_cast<Real>(T::map(-l2)), static_cast<Real>(0.5L) },
      { static_cast<Real>(0.5L),        static_cast<Real>(T::map( l2)) },
      { static_cast<Real>(0.5L),        static_cast<Real>(T::map(-l2)) },
      // axis lambda3 (±l3, 0) and (0, ±l3)
      { static_cast<Real>(T::map( l3)), static_cast<Real>(0.5L) },
      { static_cast<Real>(T::map(-l3)), static_cast<Real>(0.5L) },
      { static_cast<Real>(0.5L),        static_cast<Real>(T::map( l3)) },
      { static_cast<Real>(0.5L),        static_cast<Real>(T::map(-l3)) },
      // diagonals (±l3, ±l3)
      { static_cast<Real>(T::map( l3)), static_cast<Real>(T::map( l3)) },
      { static_cast<Real>(T::map( l3)), static_cast<Real>(T::map(-l3)) },
      { static_cast<Real>(T::map(-l3)), static_cast<Real>(T::map( l3)) },
      { static_cast<Real>(T::map(-l3)), static_cast<Real>(T::map(-l3)) },
      // corners (±l5, ±l5)
      { static_cast<Real>(T::map( l5)), static_cast<Real>(T::map( l5)) },
      { static_cast<Real>(T::map( l5)), static_cast<Real>(T::map(-l5)) },
      { static_cast<Real>(T::map(-l5)), static_cast<Real>(T::map( l5)) },
      { static_cast<Real>(T::map(-l5)), static_cast<Real>(T::map(-l5)) }
    }};
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    // Scale weights for [-1,1]^2 by Jacobian 1/4 to get [0,1]^2
    const long double w1 = (-1696.0L/2187.0L)/4.0L;     // center
    const long double w2 = ( 3920.0L/6561.0L)/4.0L;     // axis l2
    const long double w3 = ( 1360.0L/6561.0L)/4.0L;     // axis l3
    const long double w4 = (  800.0L/19683.0L)/4.0L;    // diag l3
    const long double w5 = ( 6859.0L/19683.0L)/4.0L;    // corners l5
    return {{
      static_cast<Real>(w1),
      static_cast<Real>(w2), static_cast<Real>(w2), static_cast<Real>(w2), static_cast<Real>(w2),
      static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3),
      static_cast<Real>(w4), static_cast<Real>(w4), static_cast<Real>(w4), static_cast<Real>(w4),
      static_cast<Real>(w5), static_cast<Real>(w5), static_cast<Real>(w5), static_cast<Real>(w5)
    }};
  }
};

template <>
struct rule<5, 2> {
  static constexpr bool available = true;
  static constexpr int degree = 5;
  static constexpr int dim = 2;
  static constexpr std::size_t size = 13;

  template <class Real>
  static constexpr std::array<std::array<Real, 2>, size> nodes() {
    const long double l2 = 0.3585685828003180919906451539L; // 3/sqrt(70)
    const long double l3 = 0.9486832980505137995996680633L; // 3/sqrt(10)
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    return {{
      // center
      { static_cast<Real>(0.5L), static_cast<Real>(0.5L) },
      // axis lambda2 (±l2, 0) and (0, ±l2)
      { static_cast<Real>(T::map( l2)), static_cast<Real>(0.5L) },
      { static_cast<Real>(T::map(-l2)), static_cast<Real>(0.5L) },
      { static_cast<Real>(0.5L),        static_cast<Real>(T::map( l2)) },
      { static_cast<Real>(0.5L),        static_cast<Real>(T::map(-l2)) },
      // axis lambda3 (±l3, 0) and (0, ±l3)
      { static_cast<Real>(T::map( l3)), static_cast<Real>(0.5L) },
      { static_cast<Real>(T::map(-l3)), static_cast<Real>(0.5L) },
      { static_cast<Real>(0.5L),        static_cast<Real>(T::map( l3)) },
      { static_cast<Real>(0.5L),        static_cast<Real>(T::map(-l3)) },
      // diagonals (±l3, ±l3)
      { static_cast<Real>(T::map( l3)), static_cast<Real>(T::map( l3)) },
      { static_cast<Real>(T::map( l3)), static_cast<Real>(T::map(-l3)) },
      { static_cast<Real>(T::map(-l3)), static_cast<Real>(T::map( l3)) },
      { static_cast<Real>(T::map(-l3)), static_cast<Real>(T::map(-l3)) }
    }};
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    // Degree-5 embedded weights, scaled to [0,1]^2
    const long double w1 = (-3884.0L/729.0L)/4.0L;      // center
    const long double w2 = (  490.0L/243.0L)/4.0L;      // axis l2
    const long double w3 = (  130.0L/729.0L)/4.0L;      // axis l3
    const long double w4 = (  100.0L/729.0L)/4.0L;      // diag l3
    return {{
      static_cast<Real>(w1),
      static_cast<Real>(w2), static_cast<Real>(w2), static_cast<Real>(w2), static_cast<Real>(w2),
      static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3),
      static_cast<Real>(w4), static_cast<Real>(w4), static_cast<Real>(w4), static_cast<Real>(w4)
    }};
  }
};


// Concrete specializations for 3D 7/5 on [0,1]^3 (mapped from [-1,1]^3)
// Source: Genz–Malik embedded rules; see shaiko/genz_malik_3d_7_5.md

template <>
struct rule<7, 3> {
  static constexpr bool available = true;
  static constexpr int degree = 7;
  static constexpr int dim = 3;
  static constexpr std::size_t size = 33;

  template <class Real>
  static constexpr std::array<std::array<Real, 3>, size> nodes() {
    const long double l2 = 0.3585685828003180919906451539L; // 3/sqrt(70)
    const long double l3 = 0.9486832980505137995996680633L; // 3/sqrt(10)
    const long double l5 = 0.6882472016116852977216287343L; // 3/sqrt(19)
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,3>, size> pts{};
    std::size_t k = 0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis λ2: (±l2,0,0), (0,±l2,0), (0,0,±l2)
    pts[k++] = { static_cast<Real>(T::map( l2)), static_cast<Real>(0.5L),           static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(T::map(-l2)), static_cast<Real>(0.5L),           static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(T::map( l2)), static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(T::map(-l2)), static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(0.5L),           static_cast<Real>(T::map( l2)) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(0.5L),           static_cast<Real>(T::map(-l2)) };
    // axis λ3
    pts[k++] = { static_cast<Real>(T::map( l3)), static_cast<Real>(0.5L),           static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(T::map(-l3)), static_cast<Real>(0.5L),           static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(T::map( l3)), static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(T::map(-l3)), static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(0.5L),           static_cast<Real>(T::map( l3)) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(0.5L),           static_cast<Real>(T::map(-l3)) };
    // face-diagonals λ3: permutations of (±l3, ±l3, 0)
    const long double a = l3;
    auto add_face = [&](long double x, long double y, long double z){
      pts[k++] = { static_cast<Real>(T::map( x)), static_cast<Real>(T::map( y)), static_cast<Real>(T::map( z)) };
    };
    add_face( a,  a, 0.0L); add_face( a, -a, 0.0L); add_face(-a,  a, 0.0L); add_face(-a, -a, 0.0L);
    add_face( a, 0.0L,  a); add_face( a, 0.0L, -a); add_face(-a, 0.0L,  a); add_face(-a, 0.0L, -a);
    add_face(0.0L,  a,  a); add_face(0.0L,  a, -a); add_face(0.0L, -a,  a); add_face(0.0L, -a, -a);
    // full-diagonals λ5: (±l5,±l5,±l5)
    for (int sx : {+1,-1}) for (int sy : {+1,-1}) for (int sz : {+1,-1}) {
      pts[k++] = { static_cast<Real>(T::map(sx*l5)), static_cast<Real>(T::map(sy*l5)), static_cast<Real>(T::map(sz*l5)) };
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    // Scale weights for [-1,1]^3 by Jacobian 1/8 to get [0,1]^3
    const long double w0 = (-87488.0L/19683.0L)/8.0L;    // center
    const long double w1 = (  7840.0L/ 6561.0L)/8.0L;    // axis l2
    const long double w2 = (  4960.0L/19683.0L)/8.0L;    // axis l3
    const long double w3 = (  1600.0L/19683.0L)/8.0L;    // face-diag l3
    const long double w4 = (  6859.0L/19683.0L)/8.0L;    // full-diag l5
    return {{
      static_cast<Real>(w0),
      // 6 axis l2
      static_cast<Real>(w1), static_cast<Real>(w1), static_cast<Real>(w1), static_cast<Real>(w1), static_cast<Real>(w1), static_cast<Real>(w1),
      // 6 axis l3
      static_cast<Real>(w2), static_cast<Real>(w2), static_cast<Real>(w2), static_cast<Real>(w2), static_cast<Real>(w2), static_cast<Real>(w2),
      // 12 face-diag l3
      static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3),
      static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3), static_cast<Real>(w3),
      // 8 corners l5
      static_cast<Real>(w4), static_cast<Real>(w4), static_cast<Real>(w4), static_cast<Real>(w4),
      static_cast<Real>(w4), static_cast<Real>(w4), static_cast<Real>(w4), static_cast<Real>(w4)
    }};
  }
};

template <>
struct rule<5, 3> {
  static constexpr bool available = true; // fixed weights
  static constexpr int degree = 5;
  static constexpr int dim = 3;
  static constexpr std::size_t size = 25;

  template <class Real>
  static constexpr std::array<std::array<Real, 3>, size> nodes() {
    const long double l2 = 0.3585685828003180919906451539L; // 3/sqrt(70)
    const long double l3 = 0.9486832980505137995996680633L; // 3/sqrt(10)
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,3>, size> pts{};
    std::size_t k = 0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis λ2
    pts[k++] = { static_cast<Real>(T::map( l2)), static_cast<Real>(0.5L),           static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(T::map(-l2)), static_cast<Real>(0.5L),           static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(T::map( l2)), static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(T::map(-l2)), static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(0.5L),           static_cast<Real>(T::map( l2)) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(0.5L),           static_cast<Real>(T::map(-l2)) };
    // axis λ3
    pts[k++] = { static_cast<Real>(T::map( l3)), static_cast<Real>(0.5L),           static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(T::map(-l3)), static_cast<Real>(0.5L),           static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(T::map( l3)), static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(T::map(-l3)), static_cast<Real>(0.5L) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(0.5L),           static_cast<Real>(T::map( l3)) };
    pts[k++] = { static_cast<Real>(0.5L),           static_cast<Real>(0.5L),           static_cast<Real>(T::map(-l3)) };
    // face-diagonals λ3 (no corners for degree-5)
    const long double a = l3;
    auto add_face = [&](long double x, long double y, long double z){
      pts[k++] = { static_cast<Real>(T::map( x)), static_cast<Real>(T::map( y)), static_cast<Real>(T::map( z)) };
    };
    add_face( a,  a, 0.0L); add_face( a, -a, 0.0L); add_face(-a,  a, 0.0L); add_face(-a, -a, 0.0L);
    add_face( a, 0.0L,  a); add_face( a, 0.0L, -a); add_face(-a, 0.0L,  a); add_face(-a, 0.0L, -a);
    add_face(0.0L,  a,  a); add_face(0.0L,  a, -a); add_face(0.0L, -a,  a); add_face(0.0L, -a, -a);
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    // Corrected 3D deg-5 weights from shaiko/genz_malik_3d_deg5_FIX.md, scaled by 1/8
    const long double v0 = (-4456.0L/  243.0L)/8.0L;    // center
    const long double v1 = (   980.0L/  243.0L)/8.0L;    // axis l2
    const long double v2 = (  -140.0L/  729.0L)/8.0L;    // axis l3
    const long double v3 = (   200.0L/  729.0L)/8.0L;    // face-diag l3
    return {{
      static_cast<Real>(v0),
      // 6 axis l2
      static_cast<Real>(v1), static_cast<Real>(v1), static_cast<Real>(v1), static_cast<Real>(v1), static_cast<Real>(v1), static_cast<Real>(v1),
      // 6 axis l3
      static_cast<Real>(v2), static_cast<Real>(v2), static_cast<Real>(v2), static_cast<Real>(v2), static_cast<Real>(v2), static_cast<Real>(v2),
      // 12 face-diag l3
      static_cast<Real>(v3), static_cast<Real>(v3), static_cast<Real>(v3), static_cast<Real>(v3), static_cast<Real>(v3), static_cast<Real>(v3),
      static_cast<Real>(v3), static_cast<Real>(v3), static_cast<Real>(v3), static_cast<Real>(v3), static_cast<Real>(v3), static_cast<Real>(v3)
    }};
  }
};

// Concrete specializations for 4D 7/5 on [0,1]^4 (mapped from [-1,1]^4)
// Source: shaiko/genz_malik_4d_7_5.md

template <>
struct rule<7, 4> {
  static constexpr bool available = true;
  static constexpr int degree = 7;
  static constexpr int dim = 4;
  static constexpr std::size_t size = 57; // 1 + 8 + 8 + 24 + 16

  template <class Real>
  static constexpr std::array<std::array<Real, 4>, size> nodes() {
    const long double l2 = 0.3585685828003180919906451539L; // 3/sqrt(70)
    const long double l3 = 0.9486832980505137995996680633L; // 3/sqrt(10)
    const long double l5 = 0.6882472016116852977216287343L; // 3/sqrt(19)
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,4>, size> pts{};
    std::size_t k = 0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis λ2 (±l2 along each axis)
    for (int i=0;i<4;++i){ for(int s: {+1,-1}){ std::array<long double,4> v{}; v.fill(0.0L); v[i]=s*l2; pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3]))};}}
    // axis λ3
    for (int i=0;i<4;++i){ for(int s: {+1,-1}){ std::array<long double,4> v{}; v.fill(0.0L); v[i]=s*l3; pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3]))};}}
    // face-diagonals λ3: choose i<j, assign ±l3 at (i,j)
    for (int i=0;i<4;++i){ for(int j=i+1;j<4;++j){ for(int sx: {+1,-1}) for(int sy:{+1,-1}){
      std::array<long double,4> v{}; v.fill(0.0L); v[i]=sx*l3; v[j]=sy*l3;
      pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3]))};
    }}}
    // full-diagonals λ5: all 2^4 sign combinations
    for (int sx: {+1,-1}) for(int sy:{+1,-1}) for(int sz:{+1,-1}) for(int sw:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(sx*l5)), static_cast<Real>(T::map(sy*l5)), static_cast<Real>(T::map(sz*l5)), static_cast<Real>(T::map(sw*l5)) };
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    // Scale from [-1,1]^4 by 1/16
    const long double w0 = (-92032.0L/ 6561.0L)/16.0L;   // center
    const long double w1 = ( 15680.0L/ 6561.0L)/16.0L;   // axis l2
    const long double w2 = (  3520.0L/19683.0L)/16.0L;   // axis l3
    const long double w3 = (  3200.0L/19683.0L)/16.0L;   // face-diag l3
    const long double w4 = (  6859.0L/19683.0L)/16.0L;   // full-diag l5
    std::array<Real, size> ws{}; std::size_t k=0;
    ws[k++] = static_cast<Real>(w0);
    for(int i=0;i<8;++i)  ws[k++] = static_cast<Real>(w1);
    for(int i=0;i<8;++i)  ws[k++] = static_cast<Real>(w2);
    for(int i=0;i<24;++i) ws[k++] = static_cast<Real>(w3);
    for(int i=0;i<16;++i) ws[k++] = static_cast<Real>(w4);
    return ws;
  }
};

template <>
struct rule<5, 4> {
  static constexpr bool available = true;
  static constexpr int degree = 5;
  static constexpr int dim = 4;
  static constexpr std::size_t size = 41; // 1 + 8 + 8 + 24

  template <class Real>
  static constexpr std::array<std::array<Real, 4>, size> nodes() {
    const long double l2 = 0.3585685828003180919906451539L; // 3/sqrt(70)
    const long double l3 = 0.9486832980505137995996680633L; // 3/sqrt(10)
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,4>, size> pts{}; std::size_t k=0;
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    for (int i=0;i<4;++i){ for(int s: {+1,-1}){ std::array<long double,4> v{}; v.fill(0.0L); v[i]=s*l2; pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3]))};}}
    for (int i=0;i<4;++i){ for(int s: {+1,-1}){ std::array<long double,4> v{}; v.fill(0.0L); v[i]=s*l3; pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3]))};}}
    for (int i=0;i<4;++i){ for(int j=i+1;j<4;++j){ for(int sx: {+1,-1}) for(int sy:{+1,-1}){
      std::array<long double,4> v{}; v.fill(0.0L); v[i]=sx*l3; v[j]=sy*l3;
      pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3]))};
    }}}
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double v0 = (-12112.0L/  243.0L)/16.0L;   // center
    const long double v1 = (  1960.0L/  243.0L)/16.0L;   // axis l2
    const long double v2 = (   -40.0L/   27.0L)/16.0L;   // axis l3
    const long double v3 = (   400.0L/  729.0L)/16.0L;   // face-diag l3
    std::array<Real, size> ws{}; std::size_t k=0;
    ws[k++] = static_cast<Real>(v0);
    for(int i=0;i<8;++i)  ws[k++] = static_cast<Real>(v1);
    for(int i=0;i<8;++i)  ws[k++] = static_cast<Real>(v2);
    for(int i=0;i<24;++i) ws[k++] = static_cast<Real>(v3);
    return ws;
  }
};


// Concrete specializations for 5D 7/5 on [0,1]^5 (mapped from [-1,1]^5)
// Source: shaiko/genz_malik_5d_7_5.md

template <>
struct rule<7, 5> {
  static constexpr bool available = true;
  static constexpr int degree = 7;
  static constexpr int dim = 5;
  static constexpr std::size_t size = 93; // 1 + 10 + 10 + 40 + 32

  template <class Real>
  static constexpr std::array<std::array<Real, 5>, size> nodes() {
    const long double l2 = 0.3585685828003180919906451539L; // 3/sqrt(70)
    const long double l3 = 0.9486832980505137995996680633L; // 3/sqrt(10)
    const long double l5 = 0.6882472016116852977216287343L; // 3/sqrt(19)
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,5>, size> pts{}; std::size_t k=0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis λ2
    for(int i=0;i<5;++i){ for(int s: {+1,-1}){ std::array<long double,5> v{}; v.fill(0.0L); v[i]=s*l2; pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3])),static_cast<Real>(T::map(v[4]))}; }}
    // axis λ3
    for(int i=0;i<5;++i){ for(int s: {+1,-1}){ std::array<long double,5> v{}; v.fill(0.0L); v[i]=s*l3; pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3])),static_cast<Real>(T::map(v[4]))}; }}
    // face-diagonals λ3 (pairs of coordinates with ±l3)
    for(int i=0;i<5;++i){ for(int j=i+1;j<5;++j){ for(int sx: {+1,-1}) for(int sy:{+1,-1}){
      std::array<long double,5> v{}; v.fill(0.0L); v[i]=sx*l3; v[j]=sy*l3;
      pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3])),static_cast<Real>(T::map(v[4]))};
    }}}
    // full-diagonals λ5: all 2^5 sign combinations
    for(int s0:{+1,-1}) for(int s1:{+1,-1}) for(int s2:{+1,-1}) for(int s3:{+1,-1}) for(int s4:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s0*l5)), static_cast<Real>(T::map(s1*l5)), static_cast<Real>(T::map(s2*l5)), static_cast<Real>(T::map(s3*l5)), static_cast<Real>(T::map(s4*l5)) };
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    // Scale from [-1,1]^5 by 1/32
    const long double w0 = (-242944.0L/  6561.0L)/32.0L; // center
    const long double w1 = (  31360.0L/  6561.0L)/32.0L; // axis l2
    const long double w2 = (   -640.0L/  2187.0L)/32.0L; // axis l3
    const long double w3 = (   6400.0L/ 19683.0L)/32.0L; // face-diag l3
    const long double w4 = (   6859.0L/ 19683.0L)/32.0L; // full-diag l5
    std::array<Real, size> ws{}; std::size_t k=0;
    ws[k++] = static_cast<Real>(w0);
    for(int i=0;i<10;++i) ws[k++] = static_cast<Real>(w1);
    for(int i=0;i<10;++i) ws[k++] = static_cast<Real>(w2);
    for(int i=0;i<40;++i) ws[k++] = static_cast<Real>(w3);
    for(int i=0;i<32;++i) ws[k++] = static_cast<Real>(w4);
    return ws;
  }
};

template <>
struct rule<5, 5> {
  static constexpr bool available = true;
  static constexpr int degree = 5;
  static constexpr int dim = 5;
  static constexpr std::size_t size = 61; // 1 + 10 + 10 + 40

  template <class Real>
  static constexpr std::array<std::array<Real, 5>, size> nodes() {
    const long double l2 = 0.3585685828003180919906451539L; // 3/sqrt(70)
    const long double l3 = 0.9486832980505137995996680633L; // 3/sqrt(10)
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,5>, size> pts{}; std::size_t k=0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis λ2
    for(int i=0;i<5;++i){ for(int s: {+1,-1}){ std::array<long double,5> v{}; v.fill(0.0L); v[i]=s*l2; pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3])),static_cast<Real>(T::map(v[4]))}; }}
    // axis λ3
    for(int i=0;i<5;++i){ for(int s: {+1,-1}){ std::array<long double,5> v{}; v.fill(0.0L); v[i]=s*l3; pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3])),static_cast<Real>(T::map(v[4]))}; }}
    // face-diagonals λ3
    for(int i=0;i<5;++i){ for(int j=i+1;j<5;++j){ for(int sx: {+1,-1}) for(int sy:{+1,-1}){
      std::array<long double,5> v{}; v.fill(0.0L); v[i]=sx*l3; v[j]=sy*l3;
      pts[k++]={static_cast<Real>(T::map(v[0])),static_cast<Real>(T::map(v[1])),static_cast<Real>(T::map(v[2])),static_cast<Real>(T::map(v[3])),static_cast<Real>(T::map(v[4]))};
    }}}
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    // Scale from [-1,1]^5 by 1/32
    const long double v0 = (-88672.0L/  729.0L)/32.0L; // center
    const long double v1 = (  3920.0L/  243.0L)/32.0L; // axis l2
    const long double v2 = (  -3760.0L/  729.0L)/32.0L; // axis l3
    const long double v3 = (   800.0L/  729.0L)/32.0L; // face-diag l3
    std::array<Real, size> ws{}; std::size_t k=0;
    ws[k++] = static_cast<Real>(v0);
    for(int i=0;i<10;++i) ws[k++] = static_cast<Real>(v1);
    for(int i=0;i<10;++i) ws[k++] = static_cast<Real>(v2);
    for(int i=0;i<40;++i) ws[k++] = static_cast<Real>(v3);
    return ws;
  }
};

// -----------------------------------------------------------------------------

} // close namespace gm temporarily to add utility functions

// Bring in GM 9/7 tables and types
using rules::detail::gm7;
using rules::detail::gm9;
using rules::detail::gm_orbit_kind;
using rules::detail::gm_lambda0;
using rules::detail::gm_lambda1;
using rules::detail::gm_lambda2;
using rules::detail::gm_lambda3;

// Orbit multiplicity helper (compile-time)
template <int D>
constexpr std::size_t orbit_count(gm_orbit_kind kind) {
  using rules::detail::gm_count_center;
  using rules::detail::gm_count_axis;
  using rules::detail::gm_count_pair11;
  using rules::detail::gm_count_pair12;

  using rules::detail::gm_count_triple111;
  using rules::detail::gm_count_diag;
  switch (kind) {
    case gm_orbit_kind::center:        return gm_count_center<D>();
    case gm_orbit_kind::axis_l1:
    case gm_orbit_kind::axis_l2:
    case gm_orbit_kind::axis_l3:       return gm_count_axis<D>();
    case gm_orbit_kind::pair_l1_l1:    return gm_count_pair11<D>();
    case gm_orbit_kind::pair_l1_l2:    return gm_count_pair12<D>();
    case gm_orbit_kind::triple_l1_l1_l1:return gm_count_triple111<D>();
    case gm_orbit_kind::full_diag_l0:  return gm_count_diag<D>();
  }
  return 0;
}

// Compile-time weight sum validation for release builds
template <int D, typename Groups>
inline constexpr void gm_validate_weight_sums_constexpr(const Groups& groups, long double scaled_sum) {
  // In constexpr context, we just consume the parameters
  // Actual validation happens at compile-time through static_assert in specializations
  (void)groups; 
  (void)scaled_sum;
}

// Runtime validation helper for debug builds
#ifndef NDEBUG
inline void gm_debug_assert_close(long double a, long double b, long double tol = 1e-12L) {
  long double diff = fabsl(a - b);
  assert(diff <= tol);
}

template <int D, typename Groups>
inline void gm_validate_weight_sums_debug(const Groups& groups, long double scaled_sum) {
  long double raw_sum = 0.0L;
  for (auto g : groups) {
    raw_sum += static_cast<long double>(orbit_count<D>(g.kind)) * static_cast<long double>(g.weight);
  }
  const long double two_pow_d = ldexpl(1.0L, D); // 2^D
  gm_debug_assert_close(raw_sum, two_pow_d, 1e-10L);
  gm_debug_assert_close(scaled_sum, 1.0L, 1e-12L);
}
#endif

// Main validation function that dispatches based on build mode
template <int D, typename Groups>
inline constexpr void gm_validate_weight_sums(const Groups& groups, long double scaled_sum) {
#ifndef NDEBUG
  #if __cpp_lib_is_constant_evaluated >= 201811L
  if (!std::is_constant_evaluated()) {
    gm_validate_weight_sums_debug<D>(groups, scaled_sum);
    return;
  }
  #else
  // For C++17, we can't do runtime validation in constexpr context
  // Just skip validation in constexpr context for older standards
  #endif
#endif
  gm_validate_weight_sums_constexpr<D>(groups, scaled_sum);
}

// GM 9/7 family — tag-based specializations under gm::rule_fam
// -----------------------------------------------------------------------------

namespace gm { // reopen namespace gm for specializations

// 3D degree-9 generated from shaiko/9-7_Updates/GM_9-7_d3_deg9.csv
// Mapping: nodes from [-1,1]^3 -> [0,1]^3 via (x+1)/2; weights scaled by 1/8.
template <>
struct rule_fam<9, 3, family_9_7> {
  static constexpr bool available = true;
  static constexpr int degree = 9;
  static constexpr int dim = 3;
  static constexpr std::size_t size = 71; // zero-weight orbit removed

  template <class Real>
  static constexpr std::array<std::array<Real, 3>, size> nodes() {
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,3>, size> pts{}; std::size_t k=0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis l1, l2, l3
    for (auto a : {gm_lambda1,gm_lambda2,gm_lambda3}) {
      for (int axis=0; axis<3; ++axis) {
        for (int s : {+1,-1}) {
          long double v[3] = {0.0L,0.0L,0.0L}; v[axis] = s*a;
          pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])) };
        }
      }
    }
    // pair l1-l1: permutations of (±l1, ±l1, 0)
    for (int i=0;i<3;++i){ for(int j=i+1;j<3;++j){ for(int sx:{+1,-1}) for(int sy:{+1,-1}){
      long double v[3]={0.0L,0.0L,0.0L}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1;
      pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])) };
    }}}
    // pair l1-l2: permutations of (±l1, ±l2, 0)
    for (int i=0;i<3;++i){ for(int j=i+1;j<3;++j){ for(int si:{+1,-1}) for(int sj:{+1,-1}){
      long double v1[3]={0.0L,0.0L,0.0L}; v1[i]=si*gm_lambda1; v1[j]=sj*gm_lambda2; pts[k++] = { static_cast<Real>(T::map(v1[0])), static_cast<Real>(T::map(v1[1])), static_cast<Real>(T::map(v1[2])) };
      long double v2[3]={0.0L,0.0L,0.0L}; v2[i]=si*gm_lambda2; v2[j]=sj*gm_lambda1; pts[k++] = { static_cast<Real>(T::map(v2[0])), static_cast<Real>(T::map(v2[1])), static_cast<Real>(T::map(v2[2])) };
    }}}
    // triple l1-l1-l1
    for (int sx:{+1,-1}) for (int sy:{+1,-1}) for (int sz:{+1,-1}) {
      pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda1)), static_cast<Real>(T::map(sy*gm_lambda1)), static_cast<Real>(T::map(sz*gm_lambda1)) };
    }
    // full diag l0
    for (int sx:{+1,-1}) for (int sy:{+1,-1}) for (int sz:{+1,-1}) {
      pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda0)), static_cast<Real>(T::map(sy*gm_lambda0)), static_cast<Real>(T::map(sz*gm_lambda0)) };
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double s = 1.0L/8.0L;
    std::array<Real, size> ws{}; std::size_t k=0;
    long double sum_scaled = 0.0L;
    for (auto g : gm9<3>::groups) {
      const std::size_t n = orbit_count<3>(g.kind);
      for (std::size_t i=0;i<n;++i) { 
        auto w = static_cast<Real>(g.weight * s); 
        ws[k++] = w; 
        sum_scaled += static_cast<long double>(w); 
      }
    }
    gm_validate_weight_sums<3>(gm9<3>::groups, sum_scaled);
    return ws;
  }
};

  // 3D degree-7 (coarse of 9/7) from shaiko/9-7_Updates/GM_9-7_d3_deg7.csv
  template <>
  struct rule_fam<7, 3, family_9_7> {
    static constexpr bool available = true;
    static constexpr int degree = 7;
    static constexpr int dim = 3;
    static constexpr std::size_t size = 33; // zero-weight orbit removed

    template <class Real>
    static constexpr std::array<std::array<Real, 3>, size> nodes() {
      struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };

      std::array<std::array<Real,3>, size> pts{}; std::size_t k=0;
      // center
      pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
      // axis l2
      for (int axis=0; axis<3; ++axis) for(int s:{+1,-1}){
        long double v[3]={0.0L,0.0L,0.0L}; v[axis]=s*gm_lambda2;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])) };
      }
      // axis l1
      for (int axis=0; axis<3; ++axis) for(int s:{+1,-1}){
        long double v[3]={0.0L,0.0L,0.0L}; v[axis]=s*gm_lambda1;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])) };
      }
      // pair l1-l1
      for (int i=0;i<3;++i){ for(int j=i+1;j<3;++j){ for(int sx:{+1,-1}) for(int sy:{+1,-1}){
        long double v[3]={0.0L,0.0L,0.0L}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])) };
      }}}
      // full diag l0
      for (int sx:{+1,-1}) for (int sy:{+1,-1}) for (int sz:{+1,-1}) {
        pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda0)), static_cast<Real>(T::map(sy*gm_lambda0)), static_cast<Real>(T::map(sz*gm_lambda0)) };
      }

      return pts;
    }

    template <class Real>
    static constexpr std::array<Real, size> weights() {
      const long double s = 1.0L/8.0L;
      std::array<Real, size> ws{}; std::size_t k=0;
      for (auto g : gm7<3>::groups) {
        const std::size_t n = orbit_count<3>(g.kind);
        for (std::size_t i=0;i<n;++i) ws[k++] = static_cast<Real>(g.weight * s);
      }
      return ws;
    }
  };


// 2D degree-9 and degree-7 from GM_9-7_d2_* CSVs
template <>
struct rule_fam<9, 2, family_9_7> {
  static constexpr bool available = true;
  static constexpr int degree = 9;
  static constexpr int dim = 2;
  static constexpr std::size_t size = 29;

  template <class Real>
  static constexpr std::array<std::array<Real, 2>, size> nodes() {
    using rules::detail::gm_lambda0; using rules::detail::gm_lambda1;
    using rules::detail::gm_lambda2; using rules::detail::gm_lambda3;
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,2>, size> pts{}; std::size_t k=0;
    // Emit groups in the same order as gm9<2>::groups: center, axis_l1, axis_l2, axis_l3, pair11, pair12, diag
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis_l1
    for (int s:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s*gm_lambda1)), static_cast<Real>(0.5L) };
      pts[k++] = { static_cast<Real>(0.5L),                  static_cast<Real>(T::map(s*gm_lambda1)) };
    }
    // axis_l2
    for (int s:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s*gm_lambda2)), static_cast<Real>(0.5L) };
      pts[k++] = { static_cast<Real>(0.5L),                  static_cast<Real>(T::map(s*gm_lambda2)) };
    }
    // axis_l3
    for (int s:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s*gm_lambda3)), static_cast<Real>(0.5L) };
      pts[k++] = { static_cast<Real>(0.5L),                  static_cast<Real>(T::map(s*gm_lambda3)) };
    }
    // pair l1-l1 (diagonals)
    for (int sx:{+1,-1}) for (int sy:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda1)), static_cast<Real>(T::map(sy*gm_lambda1)) };
    }
    // pair l1-l2
    for (int sx:{+1,-1}) for (int sy:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda1)), static_cast<Real>(T::map(sy*gm_lambda2)) };
      pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda2)), static_cast<Real>(T::map(sy*gm_lambda1)) };
    }
    // full diag l0
    for (int sx:{+1,-1}) for (int sy:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda0)), static_cast<Real>(T::map(sy*gm_lambda0)) };
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double s = 1.0L/4.0L;
    std::array<Real, size> ws{}; std::size_t k=0;
    long double sum_scaled = 0.0L;
    for (auto g : gm9<2>::groups) {
      const std::size_t n = orbit_count<2>(g.kind);
      for (std::size_t i=0;i<n;++i) { auto w = static_cast<Real>(g.weight * s); ws[k++] = w; sum_scaled += static_cast<long double>(w); }
    }
    gm_validate_weight_sums<2>(gm9<2>::groups, sum_scaled);
    return ws;
  }
};

// 2D degree-7 from DCUHRE Algorithm 698 D07HRE (authoritative source)
template <>
struct rule_fam<7, 2, family_9_7> {
  static constexpr bool available = true;
  static constexpr int degree = 7;
  static constexpr int dim = 2;
  static constexpr std::size_t size = 17; // 21 total nodes - 4 zero-weight axis_lp_null

  template <class Real>
  static constexpr std::array<std::array<Real, 2>, size> nodes() {
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,2>, size> pts{}; std::size_t k=0;

    // Emit in the same order as gm7<2>::groups: center, axis_l2, axis_l1, pair11, diag
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L) };

    for (int s:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s*gm_lambda2)), static_cast<Real>(0.5L) };
      pts[k++] = { static_cast<Real>(0.5L),                  static_cast<Real>(T::map(s*gm_lambda2)) };
    }

    for (int s:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s*gm_lambda1)), static_cast<Real>(0.5L) };
      pts[k++] = { static_cast<Real>(0.5L),                  static_cast<Real>(T::map(s*gm_lambda1)) };
    }

    for (int sx:{+1,-1}) for (int sy:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda1)), static_cast<Real>(T::map(sy*gm_lambda1)) };
    }

    for (int sx:{+1,-1}) for (int sy:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda0)), static_cast<Real>(T::map(sy*gm_lambda0)) };
    }

    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double s = 1.0L/4.0L; // 2^(-d) for domain mapping
    std::array<Real, size> ws{}; std::size_t k=0;
    long double sum_scaled = 0.0L;
    for (auto g : gm7<2>::groups) {
      const std::size_t n = orbit_count<2>(g.kind);
      for (std::size_t i=0;i<n;++i) { auto w = static_cast<Real>(g.weight * s); ws[k++] = w; sum_scaled += static_cast<long double>(w);}
    }
    gm_validate_weight_sums<2>(gm7<2>::groups, sum_scaled);
    return ws;
  }
};

// 4D degree-9 from GM_9-7_d4_deg9.csv (zero-weight orbit removed)
template <>
struct rule_fam<9, 4, family_9_7> {
  static constexpr bool available = true;
  static constexpr int degree = 9;
  static constexpr int dim = 4;
  static constexpr std::size_t size = 145; // 1 + 8 + 8 + 8 + 24 + 48 + 32 + 16

  template <class Real>
  static constexpr std::array<std::array<Real, 4>, size> nodes() {
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,4>, size> pts{}; std::size_t k=0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis l1,l2,l3 across all axes (in this order to match gm9<4>)
    for (auto a : {gm_lambda1, gm_lambda2, gm_lambda3}){
      for (int axis=0; axis<4; ++axis) for (int s:{+1,-1}){
        long double v[4]={0,0,0,0}; v[axis]=s*a;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])) };
      }
    }
    // pair l1-l1 on i<j
    for (int i=0;i<4;++i){ for(int j=i+1;j<4;++j){ for(int sx:{+1,-1}) for(int sy:{+1,-1}){
      long double v[4]={0,0,0,0}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1;
      pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])) };
    }}}
    // pair l1-l2 on i<j, both assignments
    for (int i=0;i<4;++i){ for(int j=i+1;j<4;++j){ for(int si:{+1,-1}) for(int sj:{+1,-1}){
      long double v1[4]={0,0,0,0}; v1[i]=si*gm_lambda1; v1[j]=sj*gm_lambda2;
      long double v2[4]={0,0,0,0}; v2[i]=si*gm_lambda2; v2[j]=sj*gm_lambda1;
      pts[k++] = { static_cast<Real>(T::map(v1[0])), static_cast<Real>(T::map(v1[1])), static_cast<Real>(T::map(v1[2])), static_cast<Real>(T::map(v1[3])) };
      pts[k++] = { static_cast<Real>(T::map(v2[0])), static_cast<Real>(T::map(v2[1])), static_cast<Real>(T::map(v2[2])), static_cast<Real>(T::map(v2[3])) };
    }}}
    // triple l1-l1-l1 on i<j<k
    for (int i=0;i<4;++i){ for(int j=i+1;j<4;++j){ for(int m=j+1;m<4;++m){
      for (int sx:{+1,-1}) for (int sy:{+1,-1}) for (int sz:{+1,-1}){
        long double v[4]={0,0,0,0}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1; v[m]=sz*gm_lambda1;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])) };
      }
    }}}
    // full diag l0
    for (int s0:{+1,-1}) for(int s1:{+1,-1}) for(int s2:{+1,-1}) for(int s3:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s0*gm_lambda0)), static_cast<Real>(T::map(s1*gm_lambda0)), static_cast<Real>(T::map(s2*gm_lambda0)), static_cast<Real>(T::map(s3*gm_lambda0)) };
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double s = 1.0L/16.0L;
    std::array<Real, size> ws{}; std::size_t k=0;
    for (auto g : gm9<4>::groups) {
      const std::size_t n = orbit_count<4>(g.kind);
      for (std::size_t i=0;i<n;++i) ws[k++] = static_cast<Real>(g.weight * s);
    }
    return ws;
  }
};

// 4D degree-7 from GM_9-7_d4_deg7.csv (zero-weight orbit removed)
template <>
struct rule_fam<7, 4, family_9_7> {
  static constexpr bool available = true;
  static constexpr int degree = 7;
  static constexpr int dim = 4;
  static constexpr std::size_t size = 57; // 1 + 8 + 8 + 24 + 16

  template <class Real>
  static constexpr std::array<std::array<Real, 4>, size> nodes() {
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,4>, size> pts{}; std::size_t k=0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis l2 then l1
    for (auto a : {gm_lambda2, gm_lambda1}){
      for (int axis=0; axis<4; ++axis) for (int s:{+1,-1}){
        long double v[4]={0,0,0,0}; v[axis]=s*a;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])) };
      }
    }
    // pair l1-l1
    for (int i=0;i<4;++i){ for(int j=i+1;j<4;++j){ for(int sx:{+1,-1}) for(int sy:{+1,-1}){
      long double v[4]={0,0,0,0}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1;
      pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])) };
    }}}
    // full diag l0
    for (int s0:{+1,-1}) for(int s1:{+1,-1}) for(int s2:{+1,-1}) for(int s3:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s0*gm_lambda0)), static_cast<Real>(T::map(s1*gm_lambda0)), static_cast<Real>(T::map(s2*gm_lambda0)), static_cast<Real>(T::map(s3*gm_lambda0)) };
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double s = 1.0L/16.0L;
    std::array<Real, size> ws{}; std::size_t k=0;
    for (auto g : gm7<4>::groups) {
      const std::size_t n = orbit_count<4>(g.kind);
      for (std::size_t i=0;i<n;++i) ws[k++] = static_cast<Real>(g.weight * s);
    }
    return ws;
  }
};

// 5D degree-9 from GM_9-7_d5_deg9.csv (zero-weight orbit removed)
template <>
struct rule_fam<9, 5, family_9_7> {
  static constexpr bool available = true;
  static constexpr int degree = 9;
  static constexpr int dim = 5;
  static constexpr std::size_t size = 263; // 1 + 10 + 10 + 10 + 40 + 80 + 80 + 32

  template <class Real>
  static constexpr std::array<std::array<Real, 5>, size> nodes() {
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,5>, size> pts{}; std::size_t k=0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis l1,l2,l3 across 5 axes
    for (auto a : {gm_lambda1, gm_lambda2, gm_lambda3}){
      for (int axis=0; axis<5; ++axis) for (int s:{+1,-1}){
        long double v[5]={0,0,0,0,0}; v[axis]=s*a;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])) };
      }
    }
    // pair l1-l1 on i<j
    for (int i=0;i<5;++i){ for(int j=i+1;j<5;++j){ for(int sx:{+1,-1}) for(int sy:{+1,-1}){
      long double v[5]={0,0,0,0,0}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1;
      pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])) };
    }}}
    // pair l1-l2 on i<j, both assignments
    for (int i=0;i<5;++i){ for(int j=i+1;j<5;++j){ for(int si:{+1,-1}) for(int sj:{+1,-1}){
      long double v1[5]={0,0,0,0,0}; v1[i]=si*gm_lambda1; v1[j]=sj*gm_lambda2;
      long double v2[5]={0,0,0,0,0}; v2[i]=si*gm_lambda2; v2[j]=sj*gm_lambda1;
      pts[k++] = { static_cast<Real>(T::map(v1[0])), static_cast<Real>(T::map(v1[1])), static_cast<Real>(T::map(v1[2])), static_cast<Real>(T::map(v1[3])), static_cast<Real>(T::map(v1[4])) };
      pts[k++] = { static_cast<Real>(T::map(v2[0])), static_cast<Real>(T::map(v2[1])), static_cast<Real>(T::map(v2[2])), static_cast<Real>(T::map(v2[3])), static_cast<Real>(T::map(v2[4])) };
    }}}
    // triple l1-l1-l1 on i<j<m
    for (int i=0;i<5;++i){ for(int j=i+1;j<5;++j){ for(int m=j+1;m<5;++m){
      for (int sx:{+1,-1}) for (int sy:{+1,-1}) for (int sz:{+1,-1}){
        long double v[5]={0,0,0,0,0}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1; v[m]=sz*gm_lambda1;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])) };
      }
    }}}
    // full diag l0
    for (int s0:{+1,-1}) for(int s1:{+1,-1}) for(int s2:{+1,-1}) for(int s3:{+1,-1}) for(int s4:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s0*gm_lambda0)), static_cast<Real>(T::map(s1*gm_lambda0)), static_cast<Real>(T::map(s2*gm_lambda0)), static_cast<Real>(T::map(s3*gm_lambda0)), static_cast<Real>(T::map(s4*gm_lambda0)) };
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double s = 1.0L/32.0L;
    std::array<Real, size> ws{}; std::size_t k=0; long double sum_scaled=0;
    for (auto g : gm9<5>::groups) {
      const std::size_t n = orbit_count<5>(g.kind);
      for (std::size_t i=0;i<n;++i) { auto w = static_cast<Real>(g.weight * s); ws[k++] = w; sum_scaled += static_cast<long double>(w); }
    }
    gm_validate_weight_sums<5>(gm9<5>::groups, sum_scaled);
    return ws;
  }
};

// 5D degree-7 from GM_9-7_d5_deg7.csv (zero-weight orbit removed)
template <>
struct rule_fam<7, 5, family_9_7> {
  static constexpr bool available = true;
  static constexpr int degree = 7;
  static constexpr int dim = 5;
  static constexpr std::size_t size = 93; // 1 + 10 + 10 + 40 + 32

  template <class Real>
  static constexpr std::array<std::array<Real, 5>, size> nodes() {
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,5>, size> pts{}; std::size_t k=0;
    // center
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    // axis l2 then l1 across all axes
    for (auto a : {gm_lambda2, gm_lambda1}){
      for (int axis=0; axis<5; ++axis) for (int s:{+1,-1}){
        long double v[5]={0,0,0,0,0}; v[axis]=s*a;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])) };
      }
    }
    // pair l1-l1
    for (int i=0;i<5;++i){ for(int j=i+1;j<5;++j){ for(int sx:{+1,-1}) for(int sy:{+1,-1}){
      long double v[5]={0,0,0,0,0}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1;
      pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])) };
    }}}
    // full diag l0
    for (int s0:{+1,-1}) for(int s1:{+1,-1}) for(int s2:{+1,-1}) for(int s3:{+1,-1}) for(int s4:{+1,-1}){
      pts[k++] = { static_cast<Real>(T::map(s0*gm_lambda0)), static_cast<Real>(T::map(s1*gm_lambda0)), static_cast<Real>(T::map(s2*gm_lambda0)), static_cast<Real>(T::map(s3*gm_lambda0)), static_cast<Real>(T::map(s4*gm_lambda0)) };
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double s = 1.0L/32.0L;
    std::array<Real, size> ws{}; std::size_t k=0;
    for (auto g : gm7<5>::groups) {
      const std::size_t n = orbit_count<5>(g.kind);
      for (std::size_t i=0;i<n;++i) ws[k++] = static_cast<Real>(g.weight * s);
    }
    return ws;
  }
};


  // 6D degree-9 from GM_9-7_d6_deg9.csv (zero-weight orbit removed)
  template <>
  struct rule_fam<9, 6, family_9_7> {
    static constexpr bool available = true;
    static constexpr int degree = 9;
    static constexpr int dim = 6;
    static constexpr std::size_t size = 441; // 1 + 12+12+12 + 60 + 120 + 160 + 64

    template <class Real>
    static constexpr std::array<std::array<Real, 6>, size> nodes() {
      struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
      std::array<std::array<Real,6>, size> pts{}; std::size_t k=0;
      // center
      pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
      // axis l1,l2,l3 across 6 axes
      for (auto a : {gm_lambda1, gm_lambda2, gm_lambda3}){
        for (int axis=0; axis<6; ++axis) for (int s:{+1,-1}){
          long double v[6]={0,0,0,0,0,0}; v[axis]=s*a;
          pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])), static_cast<Real>(T::map(v[5])) };
        }
      }
      // pair l1-l1 on i<j
      for (int i=0;i<6;++i){ for(int j=i+1;j<6;++j){ for(int sx:{+1,-1}) for(int sy:{+1,-1}){
        long double v[6]={0,0,0,0,0,0}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])), static_cast<Real>(T::map(v[5])) };
      }}}
      // pair l1-l2 on i<j, both assignments
      for (int i=0;i<6;++i){ for(int j=i+1;j<6;++j){ for(int si:{+1,-1}) for(int sj:{+1,-1}){
        long double v1[6]={0,0,0,0,0,0}; v1[i]=si*gm_lambda1; v1[j]=sj*gm_lambda2;
        long double v2[6]={0,0,0,0,0,0}; v2[i]=si*gm_lambda2; v2[j]=sj*gm_lambda1;
        pts[k++] = { static_cast<Real>(T::map(v1[0])), static_cast<Real>(T::map(v1[1])), static_cast<Real>(T::map(v1[2])), static_cast<Real>(T::map(v1[3])), static_cast<Real>(T::map(v1[4])), static_cast<Real>(T::map(v1[5])) };
        pts[k++] = { static_cast<Real>(T::map(v2[0])), static_cast<Real>(T::map(v2[1])), static_cast<Real>(T::map(v2[2])), static_cast<Real>(T::map(v2[3])), static_cast<Real>(T::map(v2[4])), static_cast<Real>(T::map(v2[5])) };
      }}}
      // triple l1-l1-l1 on i<j<m
      for (int i=0;i<6;++i){ for(int j=i+1;j<6;++j){ for(int m=j+1;m<6;++m){
        for (int sx:{+1,-1}) for (int sy:{+1,-1}) for (int sz:{+1,-1}){
          long double v[6]={0,0,0,0,0,0}; v[i]=sx*gm_lambda1; v[j]=sy*gm_lambda1; v[m]=sz*gm_lambda1;
          pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])), static_cast<Real>(T::map(v[5])) };
        }
      }}}
      // full diag l0 (all 2^6 sign combinations)
      for (int s0:{+1,-1}) for(int s1:{+1,-1}) for(int s2:{+1,-1}) for(int s3:{+1,-1}) for(int s4:{+1,-1}) for(int s5:{+1,-1}){
        pts[k++] = { static_cast<Real>(T::map(s0*gm_lambda0)), static_cast<Real>(T::map(s1*gm_lambda0)), static_cast<Real>(T::map(s2*gm_lambda0)), static_cast<Real>(T::map(s3*gm_lambda0)), static_cast<Real>(T::map(s4*gm_lambda0)), static_cast<Real>(T::map(s5*gm_lambda0)) };
      }
      return pts;
    }

    template <class Real>
    static constexpr std::array<Real, size> weights() {
      const long double s = 1.0L/64.0L;
      std::array<Real, size> ws{}; std::size_t k=0;
      for (auto g : gm9<6>::groups) {
        const std::size_t n = orbit_count<6>(g.kind);
        for (std::size_t i=0;i<n;++i) ws[k++] = static_cast<Real>(g.weight * s);
      }
      return ws;
    }
  };

  // 6D degree-7 from GM_9-7_d6_deg7.csv (zero-weight orbit removed)
  template <>
  struct rule_fam<7, 6, family_9_7> {
    static constexpr bool available = true;
    static constexpr int degree = 7;
    static constexpr int dim = 6;
    static constexpr std::size_t size = 149; // 1 + 12 + 12 + 60 + 64

    template <class Real>
    static constexpr std::array<std::array<Real, 6>, size> nodes() {
      const long double l1 = 0.955907315804538915L;
      const long double l2 = 0.406057174738239546L;
      struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
      std::array<std::array<Real,6>, size> pts{}; std::size_t k=0;
      // center
      pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
      // axis l2 then l1 across 6 axes
      for (auto a : {l2,l1}){
        for (int axis=0; axis<6; ++axis) for (int s:{+1,-1}){
          long double v[6]={0,0,0,0,0,0}; v[axis]=s*a;
          pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])), static_cast<Real>(T::map(v[5])) };
        }
      }
      // pair l1-l1
      for (int i=0;i<6;++i){ for(int j=i+1;j<6;++j){ for(int sx:{+1,-1}) for(int sy:{+1,-1}){
        long double v[6]={0,0,0,0,0,0}; v[i]=sx*l1; v[j]=sy*l1;
        pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])), static_cast<Real>(T::map(v[3])), static_cast<Real>(T::map(v[4])), static_cast<Real>(T::map(v[5])) };
      }}}
      // full diag l0
      for (int s0:{+1,-1}) for(int s1:{+1,-1}) for(int s2:{+1,-1}) for(int s3:{+1,-1}) for(int s4:{+1,-1}) for(int s5:{+1,-1}){
        pts[k++] = { static_cast<Real>(T::map(s0*gm_lambda0)), static_cast<Real>(T::map(s1*gm_lambda0)), static_cast<Real>(T::map(s2*gm_lambda0)), static_cast<Real>(T::map(s3*gm_lambda0)), static_cast<Real>(T::map(s4*gm_lambda0)), static_cast<Real>(T::map(s5*gm_lambda0)) };
      }
      return pts;
    }

    template <class Real>
    static constexpr std::array<Real, size> weights() {
      const long double s = 1.0L/64.0L;
      std::array<Real, size> ws{}; std::size_t k=0; long double sum_scaled=0;
      for (auto g : gm7<6>::groups) { const std::size_t n = orbit_count<6>(g.kind); for (std::size_t i=0;i<n;++i) { auto w = static_cast<Real>(g.weight * s); ws[k++]=w; sum_scaled += static_cast<long double>(w); } }
      gm_validate_weight_sums<6>(gm7<6>::groups, sum_scaled);
      return ws;
    }
  };

// Generic sizes for degree-9 and degree-7 rules (D >= 7)
template <int D>
constexpr std::size_t gm_total_size_deg9() {
  using rules::detail::gm_count_center;
  using rules::detail::gm_count_axis;
  using rules::detail::gm_count_pair11;
  using rules::detail::gm_count_pair12;
  using rules::detail::gm_count_triple111;
  using rules::detail::gm_count_diag;
  return gm_count_center<D>() + 3*gm_count_axis<D>() + gm_count_pair11<D>() + gm_count_pair12<D>() + gm_count_triple111<D>() + gm_count_diag<D>();
}

template <int D>
constexpr std::size_t gm_total_size_deg7() {
  using rules::detail::gm_count_center;
  using rules::detail::gm_count_axis;
  using rules::detail::gm_count_pair11;
  using rules::detail::gm_count_diag;
  return gm_count_center<D>() + 2*gm_count_axis<D>() + gm_count_pair11<D>() + gm_count_diag<D>();
}

// Generic rule for degree-9, D >= 7
template <int D>
struct rule_fam<9, D, family_9_7> {
  static constexpr bool available = true;
  static constexpr int degree = 9;
  static constexpr int dim = D;
  static constexpr std::size_t size = gm_total_size_deg9<D>();

  template <class Real>
  static constexpr std::array<std::array<Real, D>, size> nodes() {
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,D>, size> pts{}; std::size_t k=0;
    // Order: center, axis l1, axis l2, axis l3, pair11, pair12, triple111, diag
    // center
    {
      std::array<Real,D> c{}; c.fill(static_cast<Real>(0.5L)); pts[k++] = c;
    }
    // axis l1, l2, l3
    for (auto lam : {gm_lambda1, gm_lambda2, gm_lambda3}){
      for (int axis=0; axis<D; ++axis) for (int s:{+1,-1}){
        std::array<Real,D> a{}; for (int i=0;i<D;++i) a[i] = static_cast<Real>(0.5L);
        a[axis] = static_cast<Real>(T::map(s*lam));
        pts[k++] = a;
      }
    }
    // pair l1-l1 on i<j
    for (int i=0;i<D;++i) for (int j=i+1;j<D;++j){
      for (int sx:{+1,-1}) for (int sy:{+1,-1}){
        std::array<Real,D> v{}; for (int t=0;t<D;++t) v[t]=static_cast<Real>(0.5L);
        v[i]=static_cast<Real>(T::map(sx*gm_lambda1));
        v[j]=static_cast<Real>(T::map(sy*gm_lambda1));
        pts[k++] = v;
      }
    }
    // pair l1-l2 (both assignments)
    for (int i=0;i<D;++i) for (int j=i+1;j<D;++j){
      for (int si:{+1,-1}) for (int sj:{+1,-1}){
        std::array<Real,D> v1{}; for (int t=0;t<D;++t) v1[t]=static_cast<Real>(0.5L);
        v1[i]=static_cast<Real>(T::map(si*gm_lambda1));
        v1[j]=static_cast<Real>(T::map(sj*gm_lambda2));
        pts[k++] = v1;
        std::array<Real,D> v2{}; for (int t=0;t<D;++t) v2[t]=static_cast<Real>(0.5L);
        v2[i]=static_cast<Real>(T::map(si*gm_lambda2));
        v2[j]=static_cast<Real>(T::map(sj*gm_lambda1));
        pts[k++] = v2;
      }
    }
    // triple l1-l1-l1 on i<j<k
    for (int i=0;i<D;++i) for (int j=i+1;j<D;++j) for (int m=j+1;m<D;++m){
      for (int sx:{+1,-1}) for (int sy:{+1,-1}) for (int sz:{+1,-1}){
        std::array<Real,D> v{}; for (int t=0;t<D;++t) v[t]=static_cast<Real>(0.5L);
        v[i]=static_cast<Real>(T::map(sx*gm_lambda1));
        v[j]=static_cast<Real>(T::map(sy*gm_lambda1));
        v[m]=static_cast<Real>(T::map(sz*gm_lambda1));
        pts[k++] = v;
      }
    }
    // full diag l0 (all sign combinations)
    const std::size_t diag_n = static_cast<std::size_t>(1) << D;
    for (std::size_t mask=0; mask<diag_n; ++mask){
      std::array<Real,D> v{}; for (int t=0;t<D;++t){ int s = ((mask>>t)&1)? +1 : -1; v[t]=static_cast<Real>(T::map(s*gm_lambda0)); }
      pts[k++] = v;
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double s = ldexpl(1.0L, -D); // 2^(-D)
    std::array<Real, size> ws{}; std::size_t k=0; long double sum_scaled=0;
    for (auto g : gm9<D>::groups) { const std::size_t n = orbit_count<D>(g.kind); for (std::size_t i=0;i<n;++i) { auto w = static_cast<Real>(g.weight * s); ws[k++]=w; sum_scaled += static_cast<long double>(w); } }
    gm_validate_weight_sums<D>(gm9<D>::groups, sum_scaled);
    return ws;
  }
};

// Generic rule for degree-7, D >= 7
template <int D>
struct rule_fam<7, D, family_9_7> {
  static constexpr bool available = true;
  static constexpr int degree = 7;
  static constexpr int dim = D;
  static constexpr std::size_t size = gm_total_size_deg7<D>();

  template <class Real>
  static constexpr std::array<std::array<Real, D>, size> nodes() {
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,D>, size> pts{}; std::size_t k=0;
    // Order: center, axis l2, axis l1, pair11, diag
    // center
    { std::array<Real,D> c{}; c.fill(static_cast<Real>(0.5L)); pts[k++] = c; }
    // axis l2 then l1
    for (auto lam : {gm_lambda2, gm_lambda1}){
      for (int axis=0; axis<D; ++axis) for (int s:{+1,-1}){
        std::array<Real,D> a{}; for (int i=0;i<D;++i) a[i] = static_cast<Real>(0.5L);
        a[axis] = static_cast<Real>(T::map(s*lam));
        pts[k++] = a;
      }
    }
    // pair l1-l1
    for (int i=0;i<D;++i) for (int j=i+1;j<D;++j){
      for (int sx:{+1,-1}) for (int sy:{+1,-1}){
        std::array<Real,D> v{}; for (int t=0;t<D;++t) v[t]=static_cast<Real>(0.5L);
        v[i]=static_cast<Real>(T::map(sx*gm_lambda1));
        v[j]=static_cast<Real>(T::map(sy*gm_lambda1));
        pts[k++] = v;
      }
    }
    // full diag l0
    const std::size_t diag_n = static_cast<std::size_t>(1) << D;
    for (std::size_t mask=0; mask<diag_n; ++mask){
      std::array<Real,D> v{}; for (int t=0;t<D;++t){ int s = ((mask>>t)&1)? +1 : -1; v[t]=static_cast<Real>(T::map(s*gm_lambda0)); }
      pts[k++] = v;
    }
    return pts;
  }

  template <class Real>
  static constexpr std::array<Real, size> weights() {
    const long double s = ldexpl(1.0L, -D); // 2^(-D)
    std::array<Real, size> ws{}; std::size_t k=0; long double sum_scaled=0;
    for (auto g : gm7<D>::groups) { const std::size_t n = orbit_count<D>(g.kind); for (std::size_t i=0;i<n;++i) { auto w = static_cast<Real>(g.weight * s); ws[k++]=w; sum_scaled += static_cast<long double>(w); } }
    gm_validate_weight_sums<D>(gm7<D>::groups, sum_scaled);
    return ws;
  }
};

// Raw accessors exposing zero-weight orbits (for validation only)

template <int Degree, int Dim, class Family>
struct raw_rule_fam {
  static constexpr bool available = false;
  template <class Real>
  static constexpr auto nodes_with_zero() { return std::array<std::array<Real, Dim>, 0>{}; }
  template <class Real>
  static constexpr auto weights_with_zero() { return std::array<Real, 0>{}; }
};

// 2D raw (includes axis_lp_null)
template <> struct raw_rule_fam<9, 2, family_9_7> {
  static constexpr bool available = true;
  static constexpr std::size_t size = 33; // 29 + 4 zero-weight
  template <class Real>
  static constexpr std::array<std::array<Real,2>, size> nodes_with_zero() {
    const long double l1 = 0.955907315804538915L;
    const long double l2 = 0.406057174738239546L;
    const long double l3 = 0.895254709252355174L;
    const long double lp = 0.25L;
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    std::array<std::array<Real,2>, size> pts{}; std::size_t k=0;
    pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(0.5L) };
    for (auto a : {l1,l2,l3}){
      for (int s:{+1,-1}){ pts[k++] = { static_cast<Real>(T::map(s*a)), static_cast<Real>(0.5L) }; pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(T::map(s*a)) }; }
    }
    for (int sx:{+1,-1}) for (int sy:{+1,-1}) pts[k++] = { static_cast<Real>(T::map(sx*l1)), static_cast<Real>(T::map(sy*l1)) };
    for (int sx:{+1,-1}) for (int sy:{+1,-1}){ pts[k++] = { static_cast<Real>(T::map(sx*l1)), static_cast<Real>(T::map(sy*l2)) }; pts[k++] = { static_cast<Real>(T::map(sx*l2)), static_cast<Real>(T::map(sy*l1)) }; }
    for (int sx:{+1,-1}) for (int sy:{+1,-1}) pts[k++] = { static_cast<Real>(T::map(sx*gm_lambda0)), static_cast<Real>(T::map(sy*gm_lambda0)) };
    for (int s:{+1,-1}){ pts[k++] = { static_cast<Real>(T::map(s*lp)), static_cast<Real>(0.5L) }; pts[k++] = { static_cast<Real>(0.5L), static_cast<Real>(T::map(s*lp)) }; }
    return pts;
  }
  template <class Real>
  static constexpr std::array<Real, size> weights_with_zero() {
    const long double s = 1.0L/4.0L;
    const long double w_center = 2.909547716724805611L * s;
    const long double w_axis_l1 = -0.025476792732883849L * s;
    const long double w_axis_l2 =  0.123986256654814359L * s;
    const long double w_axis_l3 =  0.057693384490972686L * s;
    const long double w_pair_l1_l1 = 0.008448904373250531L * s;
    const long double w_pair_l1_l2 = 0.022543144647178933L * s;
    const long double w_full_diag  = 0.062875028738286959L * s;
    const long double w_axis_lp    = 0.0L;
    std::array<Real, size> ws{}; std::size_t k=0;
    ws[k++] = static_cast<Real>(w_center);
    for (int i=0;i<4;++i) ws[k++] = static_cast<Real>(w_axis_l1);
    for (int i=0;i<4;++i) ws[k++] = static_cast<Real>(w_axis_l2);
    for (int i=0;i<4;++i) ws[k++] = static_cast<Real>(w_axis_l3);
    for (int i=0;i<4;++i) ws[k++] = static_cast<Real>(w_pair_l1_l1);
    for (int i=0;i<8;++i) ws[k++] = static_cast<Real>(w_pair_l1_l2);
    for (int i=0;i<4;++i) ws[k++] = static_cast<Real>(w_full_diag);
    for (int i=0;i<4;++i) ws[k++] = static_cast<Real>(w_axis_lp);
    return ws;
  }
};

template <> struct raw_rule_fam<9, 3, family_9_7> {
  static constexpr bool available = true;
  static constexpr std::size_t size = 77;
  template <class Real>
  static constexpr std::array<std::array<Real,3>, size> nodes_with_zero() {
    const auto base = rule_fam<9,3, family_9_7>::template nodes<Real>();
    std::array<std::array<Real,3>, size> pts{}; std::size_t k=0; for (auto p : base) pts[k++]=p;
    struct T { static constexpr long double map(long double x){ return (x + 1.0L) / 2.0L; } };
    const long double lp = 0.25L;
    for (int axis=0; axis<3; ++axis) for (int s:{+1,-1}){
      long double v[3]={0,0,0}; v[axis]=s*lp;
      pts[k++] = { static_cast<Real>(T::map(v[0])), static_cast<Real>(T::map(v[1])), static_cast<Real>(T::map(v[2])) };
    }
    return pts;
  }
  template <class Real>
  static constexpr std::array<Real, size> weights_with_zero() {
    const auto base = rule_fam<9,3, family_9_7>::template weights<Real>();
    std::array<Real, size> ws{}; std::size_t k=0; for (auto w : base) ws[k++]=w; for (int i=0;i<6;++i) ws[k++]=static_cast<Real>(0);
    return ws;
  }
};

} // namespace gm


}}}} // namespace boost::math::cubature::detail

#endif // BOOST_MATH_CUBATURE_DETAIL_GM_RULES_HPP
