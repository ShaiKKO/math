// Copyright 2025 Math/Ripple
// 5D Genz-Malik 7/5 embedded pair raw rules
// Auto-generated from genz_malik_5d_7_5_nodes_deg7.csv and genz_malik_5d_7_5_nodes_deg5.csv

#ifndef BOOST_MATH_CUBATURE_DETAIL_GM_RULES_5D_HPP
#define BOOST_MATH_CUBATURE_DETAIL_GM_RULES_5D_HPP

#include <array>
#include <cstddef>

// Degree-7 raw rule for 5D (93 nodes)
template <> struct raw_rule_fam<7, 5, family_7_5> {
  static constexpr bool available = true;
  static constexpr std::size_t size = 93;
  
  template <class Real>
  static constexpr std::array<std::array<Real,5>, size> nodes_with_zero() {
    return {{
      // Group: center (1 nodes)
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},

      // Group: axis λ2 (10 nodes)
      {{static_cast<Real>(-0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.358568582800318059L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.358568582800318059L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.358568582800318059L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.358568582800318059L)}},

      // Group: axis λ3 (10 nodes)
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},

      // Group: face-diagonal λ3 (40 nodes)
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.948683298050513768L)}},

      // Group: full-diagonal λ5 (32 nodes)
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(-0.688247201611685178L)}},
      {{static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L), static_cast<Real>(0.688247201611685178L)}}
    }};
  }
  
  template <class Real>
  static constexpr std::array<Real, size> weights() {
    return {{
      // Group: center
      static_cast<Real>(-37.028501752781586731L),
      // Group: axis λ2
      static_cast<Real>(4.779759183051363713L),
      static_cast<Real>(4.779759183051363713L),
      static_cast<Real>(4.779759183051363713L),
      static_cast<Real>(4.779759183051363713L),
      static_cast<Real>(4.779759183051363713L),
      static_cast<Real>(4.779759183051363713L),
      static_cast<Real>(4.779759183051363713L),
      static_cast<Real>(4.779759183051363713L),
      static_cast<Real>(4.779759183051363713L),
      static_cast<Real>(4.779759183051363713L),
      // Group: axis λ3
      static_cast<Real>(-0.292638317329675379L),
      static_cast<Real>(-0.292638317329675379L),
      static_cast<Real>(-0.292638317329675379L),
      static_cast<Real>(-0.292638317329675379L),
      static_cast<Real>(-0.292638317329675379L),
      static_cast<Real>(-0.292638317329675379L),
      static_cast<Real>(-0.292638317329675379L),
      static_cast<Real>(-0.292638317329675379L),
      static_cast<Real>(-0.292638317329675379L),
      static_cast<Real>(-0.292638317329675379L),
      // Group: face-diagonal λ3
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      static_cast<Real>(0.325153685921861502L),
      // Group: full-diagonal λ5
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L),
      static_cast<Real>(0.348473301834069993L)
    }};
  }
};

// Degree-5 raw rule for 5D (61 nodes)
template <> struct raw_rule_fam<5, 5, family_7_5> {
  static constexpr bool available = true;
  static constexpr std::size_t size = 61;
  
  template <class Real>
  static constexpr std::array<std::array<Real,5>, size> nodes_with_zero() {
    return {{
      // Group: center (1 nodes)
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},

      // Group: axis λ2 (10 nodes)
      {{static_cast<Real>(-0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.358568582800318059L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.358568582800318059L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.358568582800318059L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.358568582800318059L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.358568582800318059L)}},

      // Group: axis λ3 (10 nodes)
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},

      // Group: face-diagonal λ3 (40 nodes)
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.948683298050513768L), static_cast<Real>(0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(-0.948683298050513768L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.948683298050513768L), static_cast<Real>(0.948683298050513768L)}}
    }};
  }
  
  template <class Real>
  static constexpr std::array<Real, size> weights() {
    return {{
      // Group: center
      static_cast<Real>(-121.635116598079562777L),
      // Group: axis λ2
      static_cast<Real>(16.131687242798353310L),
      static_cast<Real>(16.131687242798353310L),
      static_cast<Real>(16.131687242798353310L),
      static_cast<Real>(16.131687242798353310L),
      static_cast<Real>(16.131687242798353310L),
      static_cast<Real>(16.131687242798353310L),
      static_cast<Real>(16.131687242798353310L),
      static_cast<Real>(16.131687242798353310L),
      static_cast<Real>(16.131687242798353310L),
      static_cast<Real>(16.131687242798353310L),
      // Group: axis λ3
      static_cast<Real>(-5.157750342935528032L),
      static_cast<Real>(-5.157750342935528032L),
      static_cast<Real>(-5.157750342935528032L),
      static_cast<Real>(-5.157750342935528032L),
      static_cast<Real>(-5.157750342935528032L),
      static_cast<Real>(-5.157750342935528032L),
      static_cast<Real>(-5.157750342935528032L),
      static_cast<Real>(-5.157750342935528032L),
      static_cast<Real>(-5.157750342935528032L),
      static_cast<Real>(-5.157750342935528032L),
      // Group: face-diagonal λ3
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L),
      static_cast<Real>(1.097393689986282617L)
    }};
  }
};

#endif // BOOST_MATH_CUBATURE_DETAIL_GM_RULES_5D_HPP
