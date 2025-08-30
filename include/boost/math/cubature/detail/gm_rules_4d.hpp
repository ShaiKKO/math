// Copyright 2025 Math/Ripple
// 4D Genz-Malik 9/7 embedded pair raw rules
// Auto-generated from GM_9-7_d4_deg9.csv and GM_9-7_d4_deg7.csv

// Degree-9 raw rule for 4D (153 nodes)
template <> struct raw_rule_fam<9, 4, family_9_7> {
  static constexpr bool available = true;
  static constexpr std::size_t size = 153;
  
  template <class Real>
  static constexpr std::array<std::array<Real,4>, size> nodes_with_zero() {
    return {{
      // Group 1: center (1 node)
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      
      // Group 2: axis_l1 (8 nodes)
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      
      // Group 3: axis_l2 (8 nodes)  
      {{static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L)}},
      
      // Group 4: axis_l3 (8 nodes)
      {{static_cast<Real>(-0.895254709252355174L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.895254709252355174L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.895254709252355174L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.895254709252355174L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.895254709252355174L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.895254709252355174L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.895254709252355174L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.895254709252355174L)}},
      
      // Group 5: axis_l4 (8 nodes)
      {{static_cast<Real>(-0.25L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.25L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.25L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.25L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.25L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.25L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.25L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.25L)}},
      
      // Group 6: face_ll (24 nodes)
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      
      // Group 7: face_ls (48 nodes) - l1 with l2
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L)}},
      {{static_cast<Real>(-0.406057174738239546L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.406057174738239546L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L)}},
      {{static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.406057174738239546L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.406057174738239546L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.406057174738239546L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.406057174738239546L)}},
      {{static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(-0.406057174738239546L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.0L), static_cast<Real>(0.406057174738239546L), static_cast<Real>(0.955907315804538915L)}},
      
      // Group 8: volume_lll (32 nodes)
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(-0.955907315804538915L)}},
      {{static_cast<Real>(0.0L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L), static_cast<Real>(0.955907315804538915L)}},
      
      // Group 9: volume_ssss (16 nodes)
      {{static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L)}},
      {{static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L)}},
      {{static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L)}},
      {{static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L)}},
      {{static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L)}},
      {{static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L)}},
      {{static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L)}},
      {{static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L)}},
      {{static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L)}},
      {{static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L)}},
      {{static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L)}},
      {{static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L)}},
      {{static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(-0.686075797561756295L)}},
      {{static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L), static_cast<Real>(0.686075797561756295L)}},
      {{static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(-0.686075797561756295L)}},
      {{static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L), static_cast<Real>(0.686075797561756295L)}},
    }};
  }
  
  template <class Real>
  static constexpr std::array<Real, size> weights_with_zero() {
    std::array<Real, size> ws{};
    
    // Group 1: center
    ws[0] = static_cast<Real>(14.870490543857641796L);
    
    // Group 2: axis_l1 (8 nodes)
    for (int i = 1; i <= 8; ++i) ws[i] = static_cast<Real>(-0.128194450104532093L);
    
    // Group 3: axis_l2 (8 nodes)
    for (int i = 9; i <= 16; ++i) ws[i] = static_cast<Real>(0.033813678066098640L);
    
    // Group 4: axis_l3 (8 nodes)
    for (int i = 17; i <= 24; ++i) ws[i] = static_cast<Real>(0.057693384490972686L);
    
    // Group 5: axis_l4 (8 nodes) - zero weights
    for (int i = 25; i <= 32; ++i) ws[i] = static_cast<Real>(0.0L);
    
    // Group 6: face_ll (24 nodes)
    for (int i = 33; i <= 56; ++i) ws[i] = static_cast<Real>(0.001365391469893991L);
    
    // Group 7: face_ls (48 nodes)
    for (int i = 57; i <= 104; ++i) ws[i] = static_cast<Real>(0.022543144647178933L);
    
    // Group 8: volume_lll (32 nodes)
    for (int i = 105; i <= 136; ++i) ws[i] = static_cast<Real>(0.001770878225839135L);
    
    // Group 9: volume_ssss (16 nodes)
    for (int i = 137; i <= 152; ++i) ws[i] = static_cast<Real>(0.015718757184571740L);
    
    return ws;
  }
};

// Degree-7 raw rule for 4D (shares nodes with degree-9)
template <> struct raw_rule_fam<7, 4, family_9_7> {
  static constexpr bool available = true;
  static constexpr std::size_t size = 153; // Same nodes as degree-9
  
  template <class Real>
  static constexpr std::array<std::array<Real,4>, size> nodes_with_zero() {
    return raw_rule_fam<9, 4, family_9_7>::template nodes_with_zero<Real>();
  }
  
  template <class Real>
  static constexpr std::array<Real, size> weights_with_zero() {
    std::array<Real, size> ws{};
    
    // Group 1: center
    ws[0] = static_cast<Real>(14.340132246424625251L);
    
    // Group 2: axis_l1 (8 nodes) -> axis_l2 in degree-7
    for (int i = 1; i <= 8; ++i) ws[i] = static_cast<Real>(0.123593980320432464L);
    
    // Group 3: axis_l2 (8 nodes) -> axis_l4 in degree-7  
    for (int i = 9; i <= 16; ++i) ws[i] = static_cast<Real>(0.010369774823699232L);
    
    // Group 4: axis_l3 (8 nodes) - zero weights in degree-7
    for (int i = 17; i <= 24; ++i) ws[i] = static_cast<Real>(0.0L);
    
    // Group 5: axis_l4 (8 nodes) -> axis_l5 in degree-7
    for (int i = 25; i <= 32; ++i) ws[i] = static_cast<Real>(0.009708933337374199L);
    
    // Group 6: face_ll (24 nodes) -> face_ss in degree-7
    for (int i = 33; i <= 56; ++i) ws[i] = static_cast<Real>(0.009708933337374199L);
    
    // Group 7: face_ls (first 8 nodes) -> face_ls in degree-7
    for (int i = 57; i <= 64; ++i) ws[i] = static_cast<Real>(0.022196457020333755L);
    
    // Remaining nodes have zero weight for degree-7
    for (int i = 65; i <= 152; ++i) ws[i] = static_cast<Real>(0.0L);
    
    return ws;
  }
};