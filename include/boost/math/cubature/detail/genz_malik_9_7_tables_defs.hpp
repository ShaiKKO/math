// Copyright 2025 Colin MacRitchie/Ripple Group  
// Out-of-class definitions for C++14 compatibility

#ifndef BOOST_MATH_CUBATURE_DETAIL_GENZ_MALIK_9_7_TABLES_DEFS_HPP
#define BOOST_MATH_CUBATURE_DETAIL_GENZ_MALIK_9_7_TABLES_DEFS_HPP

#include "genz_malik_9_7_tables.hpp"

namespace boost { namespace math { namespace cubature { namespace rules { namespace detail {

// C++14 requires out-of-class definitions for static constexpr arrays
// No template<> needed for member definitions

constexpr std::array<gm_group, 5> gm7<2>::groups;
constexpr std::array<gm_group, 5> gm7<3>::groups;
constexpr std::array<gm_group, 5> gm7<4>::groups;
constexpr std::array<gm_group, 5> gm7<5>::groups;
constexpr std::array<gm_group, 5> gm7<6>::groups;
constexpr std::array<gm_group, 5> gm7<7>::groups;
constexpr std::array<gm_group, 5> gm7<8>::groups;
constexpr std::array<gm_group, 5> gm7<9>::groups;
constexpr std::array<gm_group, 5> gm7<10>::groups;
constexpr std::array<gm_group, 5> gm7<11>::groups;
constexpr std::array<gm_group, 5> gm7<12>::groups;
constexpr std::array<gm_group, 5> gm7<13>::groups;
constexpr std::array<gm_group, 5> gm7<14>::groups;
constexpr std::array<gm_group, 5> gm7<15>::groups;

// gm9 arrays have 7 groups for dim 2, but 8 for dim 3+
constexpr std::array<gm_group, 7> gm9<2>::groups;
constexpr std::array<gm_group, 8> gm9<3>::groups;
constexpr std::array<gm_group, 8> gm9<4>::groups;
constexpr std::array<gm_group, 8> gm9<5>::groups;
constexpr std::array<gm_group, 8> gm9<6>::groups;
constexpr std::array<gm_group, 8> gm9<7>::groups;
constexpr std::array<gm_group, 8> gm9<8>::groups;
constexpr std::array<gm_group, 8> gm9<9>::groups;
constexpr std::array<gm_group, 8> gm9<10>::groups;
constexpr std::array<gm_group, 8> gm9<11>::groups;
constexpr std::array<gm_group, 8> gm9<12>::groups;
constexpr std::array<gm_group, 8> gm9<13>::groups;
constexpr std::array<gm_group, 8> gm9<14>::groups;
constexpr std::array<gm_group, 8> gm9<15>::groups;

}}}}} // namespace

#endif // BOOST_MATH_CUBATURE_DETAIL_GENZ_MALIK_9_7_TABLES_DEFS_HPP