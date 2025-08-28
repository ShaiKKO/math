// Auto-generated: Genz–Malik 9/7 group weights for d=2..15 (DCUHRE-finalized)
#pragma once

#include <array>
#include <cstddef>

namespace boost { namespace math { namespace cubature { namespace rules { namespace detail {

enum class gm_orbit_kind { center, axis_l1, axis_l2, axis_l3, pair_l1_l1, pair_l1_l2, triple_l1_l1_l1, full_diag_l0 };

struct gm_group { gm_orbit_kind kind; long double weight; };

// Radii (dimension-independent):
static constexpr long double gm_lambda0 = 0.686075797561756291400285L; // sqrt(0.4707)
static constexpr long double gm_lambda1 = 0.955907315804539012385721L;
static constexpr long double gm_lambda2 = 0.406057174738239735599607L;
static constexpr long double gm_lambda3 = 0.895254709252356257641541L; // deg-9 only
static constexpr long double gm_lambda_p7 = 0.750000000000000000000000L; // 0.75 (excluded: zero-weight)
static constexpr long double gm_lambda_p9 = 0.250000000000000000000000L; // 0.25 (excluded: zero-weight)

// Orbit multiplicities on [-1,1]^D (used by expanders; compile-time functions):
template<std::size_t D> constexpr std::size_t gm_count_center()  { return 1; }
template<std::size_t D> constexpr std::size_t gm_count_axis()    { return 2*D; }
template<std::size_t D> constexpr std::size_t gm_count_pair11()  { return 2*D*(D-1); }
template<std::size_t D> constexpr std::size_t gm_count_pair12()  { return 4*D*(D-1); }
template<std::size_t D> constexpr std::size_t gm_count_triple111(){ return (D<3)?0:(4*D*(D-1)*(D-2))/3; }
template<std::size_t D> constexpr std::size_t gm_count_diag()    { return (std::size_t)1<<D; }

// ---- Degree-7 group tables (excluding axis_lp_null) ----
template<std::size_t D> struct gm7;

template<> struct gm7<2> {
  static constexpr std::size_t dim = 2;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -0.340707998597401099923398L },
    gm_group{ gm_orbit_kind::axis_l2, 0.494375921281729342873450L },
    gm_group{ gm_orbit_kind::axis_l1, 0.196822032692783977478245L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 0.038835733349496776856891L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<3> {
  static constexpr std::size_t dim = 3;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -3.135521946296881266398450L },
    gm_group{ gm_orbit_kind::axis_l2, 0.988751842563458685746901L },
    gm_group{ gm_orbit_kind::axis_l1, 0.238301131987580847528926L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 0.077671466698993553713781L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<4> {
  static constexpr std::size_t dim = 4;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -10.557884057205972236189957L },
    gm_group{ gm_orbit_kind::axis_l2, 1.977503685126917371493801L },
    gm_group{ gm_orbit_kind::axis_l1, 0.165916397179187480202727L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 0.155342933397987107427563L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<5> {
  static constexpr std::size_t dim = 5;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -28.446704976452467019745527L },
    gm_group{ gm_orbit_kind::axis_l2, 3.955007370253834742987603L },
    gm_group{ gm_orbit_kind::axis_l1, -0.289538939233573469304797L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 0.310685866795974214855125L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<6> {
  static constexpr std::size_t dim = 6;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -69.069796742618185415381275L },
    gm_group{ gm_orbit_kind::axis_l2, 7.910014740507669485975205L },
    gm_group{ gm_orbit_kind::axis_l1, -1.821821345651043798030095L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 0.621371733591948429710251L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<7> {
  static constexpr std::size_t dim = 7;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -157.521393195927286144860984L },
    gm_group{ gm_orbit_kind::axis_l2, 15.820029481015338971950411L },
    gm_group{ gm_orbit_kind::axis_l1, -6.129129625669881314901194L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 1.242743467183896859420502L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<8> {
  static constexpr std::size_t dim = 8;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -343.864438075765228042554821L },
    gm_group{ gm_orbit_kind::axis_l2, 31.640058962030677943900822L },
    gm_group{ gm_orbit_kind::axis_l1, -17.229233120075350067484396L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 2.485486934367793718841004L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<9> {
  static constexpr std::size_t dim = 9;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -725.488284044409417840047314L },
    gm_group{ gm_orbit_kind::axis_l2, 63.280117924061355887801643L },
    gm_group{ gm_orbit_kind::axis_l1, -44.400413977621875010332807L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 4.970973868735587437682008L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<10> {
  static constexpr std::size_t dim = 10;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -1486.727592924692059688513912L },
    gm_group{ gm_orbit_kind::axis_l2, 126.560235848122711775603287L },
    gm_group{ gm_orbit_kind::axis_l1, -108.684723430186099771393645L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 9.941947737471174875364016L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<11> {
  static constexpr std::size_t dim = 11;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -2965.421653621361168390954263L },
    gm_group{ gm_orbit_kind::axis_l2, 253.120471696245423551206573L },
    gm_group{ gm_orbit_kind::axis_l1, -257.137237810256899044243353L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 19.883895474942349750728031L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<12> {
  static constexpr std::size_t dim = 12;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -5755.705078987137636803937155L },
    gm_group{ gm_orbit_kind::axis_l2, 506.240943392490847102413147L },
    gm_group{ gm_orbit_kind::axis_l1, -593.810057520283197091398832L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 39.767790949884699501456063L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<13> {
  static constexpr std::size_t dim = 13;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -10842.991373864028277640283065L },
    gm_group{ gm_orbit_kind::axis_l2, 1012.481886784981694204826294L },
    gm_group{ gm_orbit_kind::axis_l1, -1346.691278840105192188621916L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 79.535581899769399002912125L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<14> {
  static constexpr std::size_t dim = 14;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -19712.860524309407371322086636L },
    gm_group{ gm_orbit_kind::axis_l2, 2024.963773569963388409652587L },
    gm_group{ gm_orbit_kind::axis_l1, -3011.524885279287980388892334L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 159.071163799538798005824251L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

template<> struct gm7<15> {
  static constexpr std::size_t dim = 15;
  static constexpr std::array<gm_group, 5> groups = {{
    gm_group{ gm_orbit_kind::center, -34206.907291385205990680620278L },
    gm_group{ gm_orbit_kind::axis_l2, 4049.927547139926776819305174L },
    gm_group{ gm_orbit_kind::axis_l1, -6659.334425756731152801081672L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 318.142327599077596011648502L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.355143312325340177772264L },
  }};
};

// ---- Degree-9 group tables (excluding axis_lp_null) ----
template<std::size_t D> struct gm9;

template<> struct gm9<2> {
  static constexpr std::size_t dim = 2;
  static constexpr std::array<gm_group, 7> groups = {{
    gm_group{ gm_orbit_kind::center, -0.361809133100747141654321L },
    gm_group{ gm_orbit_kind::axis_l1, -0.101907170931537513430354L },
    gm_group{ gm_orbit_kind::axis_l2, 0.495945026619249009960727L },
    gm_group{ gm_orbit_kind::axis_l3, 0.230773537963893934291651L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 0.033795617493002071594637L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 0.090172578588715681519611L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<3> {
  static constexpr std::size_t dim = 3;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, -1.623073849898153095369111L },
    gm_group{ gm_orbit_kind::axis_l1, -0.643019022963093756396076L },
    gm_group{ gm_orbit_kind::axis_l2, 0.631199738883635293843009L },
    gm_group{ gm_orbit_kind::axis_l3, 0.461547075927787868583302L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 0.039257183372578001728462L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 0.180345157177431363039223L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 0.014167025806713070730406L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<4> {
  static constexpr std::size_t dim = 4;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, -2.072151298277507124090396L },
    gm_group{ gm_orbit_kind::axis_l1, -2.051111201672520406019644L },
    gm_group{ gm_orbit_kind::axis_l2, 0.541018849057545135529127L },
    gm_group{ gm_orbit_kind::axis_l3, 0.923094151855575737166604L },
    gm_group{ gm_orbit_kind::pair_l1_l1, 0.021846263518303720535301L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 0.360690314354862726078445L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 0.028334051613426141460811L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<5> {
  static constexpr std::size_t dim = 5;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 3.696160518491999003279411L },
    gm_group{ gm_orbit_kind::axis_l1, -5.405696301930297466807783L },
    gm_group{ gm_orbit_kind::axis_l2, -0.360723559304360633255526L },
    gm_group{ gm_orbit_kind::axis_l3, 1.846188303711151474333207L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -0.069643679417097124772643L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 0.721380628709725452156890L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 0.056668103226852282921623L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<6> {
  static constexpr std::size_t dim = 6;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 33.151498239467218219062365L },
    gm_group{ gm_orbit_kind::axis_l1, -12.964995575216289979779574L },
    gm_group{ gm_orbit_kind::axis_l2, -3.606969633447623075138612L },
    gm_group{ gm_orbit_kind::axis_l3, 3.692376607422302948666415L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -0.365959771741603381231776L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 1.442761257419450904313780L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 0.113336206453704565843245L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<7> {
  static constexpr std::size_t dim = 7;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 136.164473525419991228806164L },
    gm_group{ gm_orbit_kind::axis_l1, -29.330507441514333525141202L },
    gm_group{ gm_orbit_kind::axis_l2, -12.984984296573049767532344L },
    gm_group{ gm_orbit_kind::axis_l3, 7.384753214844605897332829L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -1.185264369298025025836533L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 2.885522514838901808627560L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 0.226672412907409131686490L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<8> {
  static constexpr std::size_t dim = 8;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 445.111387820330774663340053L },
    gm_group{ gm_orbit_kind::axis_l1, -63.648668161932901127954592L },
    gm_group{ gm_orbit_kind::axis_l2, -37.512058652501706769574929L },
    gm_group{ gm_orbit_kind::axis_l3, 14.769506429689211794665658L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -3.277218390225686578419027L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 5.771045029677803617255120L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 0.453344825814818263372981L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<9> {
  static constexpr std::size_t dim = 9;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 1294.653113319645406772897573L },
    gm_group{ gm_orbit_kind::axis_l1, -133.645884275155724304269712L },
    gm_group{ gm_orbit_kind::axis_l2, -98.108297423714628008170339L },
    gm_group{ gm_orbit_kind::axis_l3, 29.539012859378423589331317L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -8.367816083710646210329976L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 11.542090059355607234510241L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 0.906689651629636526745961L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<10> {
  static constexpr std::size_t dim = 10;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 3501.390779851188890079818737L },
    gm_group{ gm_orbit_kind::axis_l1, -272.735347239854200491292789L },
    gm_group{ gm_orbit_kind::axis_l2, -242.384955084851684954381641L },
    gm_group{ gm_orbit_kind::axis_l3, 59.078025718756847178662634L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -20.362390773939838527643799L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 23.084180118711214469020481L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 1.813379303259273053491923L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<11> {
  static constexpr std::size_t dim = 11;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 9004.384352981886287654991202L },
    gm_group{ gm_orbit_kind::axis_l1, -541.850817432719720320156925L },
    gm_group{ gm_orbit_kind::axis_l2, -577.106630644548227784845207L },
    gm_group{ gm_orbit_kind::axis_l3, 118.156051437513694357325268L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -47.978298760916769269255289L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 46.168360237422428938040962L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 3.626758606518546106983846L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<12> {
  static constexpr std::size_t dim = 12;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 22308.813528529917561443561415L },
    gm_group{ gm_orbit_kind::axis_l1, -1047.447811919313710459585778L },
    gm_group{ gm_orbit_kind::axis_l2, -1338.886702238786171321854263L },
    gm_group{ gm_orbit_kind::axis_l3, 236.312102875027388714650536L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -110.463631947907722966445960L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 92.336720474844857876081925L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 7.253517213037092213967691L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<13> {
  static constexpr std::size_t dim = 13;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 53695.338898797787562016540905L },
    gm_group{ gm_orbit_kind::axis_l1, -1964.359840242079222845973883L },
    gm_group{ gm_orbit_kind::axis_l2, -3047.120286376951774148036226L },
    gm_group{ gm_orbit_kind::axis_l3, 472.624205750054777429301071L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -249.941332747963814788762687L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 184.673440949689715752163850L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 14.507034426074184427935383L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<14> {
  static constexpr std::size_t dim = 14;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 126269.233323465617985169471932L },
    gm_group{ gm_orbit_kind::axis_l1, -3551.591837882468574122069353L },
    gm_group{ gm_orbit_kind::axis_l2, -6832.934336552662411304727850L },
    gm_group{ gm_orbit_kind::axis_l3, 945.248411500109554858602142L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -557.910803200224367289266905L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 369.346881899379431504327700L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 29.014068852148368855870766L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

template<> struct gm9<15> {
  static constexpr std::size_t dim = 15;
  static constexpr std::array<gm_group, 8> groups = {{
    gm_group{ gm_orbit_kind::center, 291277.616281825223756672899803L },
    gm_group{ gm_orbit_kind::axis_l1, -6116.815439744370454257415755L },
    gm_group{ gm_orbit_kind::axis_l2, -15143.256200702842548626766498L },
    gm_group{ gm_orbit_kind::axis_l3, 1890.496823000219109717204284L },
    gm_group{ gm_orbit_kind::pair_l1_l1, -1231.877881809042210002016875L },
    gm_group{ gm_orbit_kind::pair_l1_l2, 738.693763798758863008655399L },
    gm_group{ gm_orbit_kind::triple_l1_l1_l1, 58.028137704296737711741532L },
    gm_group{ gm_orbit_kind::full_diag_l0, 0.251500114953147919957697L },
  }};
};

}}}}} // namespaces
