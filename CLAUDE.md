# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## CRITICAL: Code Quality Standards

**You MUST maintain Boost library maintainer-level code quality. This means:**
- Follow exact Boost.Math conventions and patterns consistently
- Ensure all code aligns with what Boost library maintainers would expect
- Prioritize robustness, portability, and performance in that order
- Use constexpr and template metaprogramming appropriately
- Maintain backward compatibility and avoid breaking changes

## Repository Overview

This is the Boost.Math library repository - a comprehensive C++14 mathematical library that provides header-only implementations of advanced mathematical functions. The library is part of the Boost C++ Libraries collection.

## Current Work Context: Genz-Malik Cubature Integration

The repository is currently on the `cubature` branch with active development on the Genz-Malik 9/7 integration tables for numerical cubature.

### Design Principles & Architecture Decisions
When making any design or architecture decisions:
- **Always prefer what Boost library maintainers would want**
- Align all choices with Boost conventions and expectations
- Maintain consistency with existing Boost.Math patterns
- Prioritize clarity and maintainability over cleverness
- Ensure all template metaprogramming follows Boost style

### GM 9/7 Implementation Requirements
**AUTHORITATIVE SOURCE**: Use `shaiko/genz_malik_9_7_tables.hpp` as the single source of truth for GM 9/7 implementation:
- Refactor `rule_fam` to generate weights/nodes generically from tables
- Apply 2^-D scaling for [0,1]^D mapping consistently
- **EXCLUDE** axis_lp_null (zero-weight orbits) from node generation
- Add comprehensive validation in debug builds
- Maintain backward compatibility across D=2–15 for degrees 7 and 9
- Use table-driven approach with `gm7<D>` and `gm9<D>` groups

### Current State
- Integrating generated Genz-Malik 9/7 tables from `shaiko/genz_malik_9_7_tables.hpp`
- Converting dimension-specific specializations to table-driven implementations using:
  - `gm7<D>::groups` and `gm9<D>::groups` for weights (from DCUHRE-finalized tables)
  - Shared radii `gm_lambda0..3` for nodes
  - Scaling factor `s = 2^(-D)` for [0,1]^D mapping
- Debug validator `gm_validate_weight_sums<D>()` checks weight sums

### Active Issues
- Build failures in `include/boost/math/cubature/detail/gm_rules.hpp` due to:
  - Duplicate/stray code blocks (especially around 3D degree-9/7 sections)
  - Template specializations outside namespace scope
  - Unused constants causing -Werror warnings (6D deg-9 nodes)
- Include path temporarily adjusted to `../../../../../shaiko/...` (to be reverted once builds pass)

### Test Files for Validation
- `test/gm97_weight_sums_test.cpp` - Validates weight sums equal 1.0 for D=2-15
- `test/gm97_exactness_test.cpp` - Tests polynomial exactness for degree-7/9 rules
- `test/gm_rules_exactness_test.cpp` - General rule exactness testing

## Build Commands

```bash
# Configure build (from repository root)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

# Build all targets
cd build && make -j

# Build specific test targets
make gm97_weight_sums_test
make gm97_exactness_test  
make gm_rules_exactness_test
make cubature_smoke_test

# Run tests
ctest                    # Run all tests
./test/gm97_weight_sums_test   # Run specific test

# Clean build
make clean

# Build with Release optimizations
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

## Code Architecture

### Cubature Module Structure
```
include/boost/math/cubature/
├── adaptive.hpp           # Adaptive cubature algorithms
├── concepts.hpp           # C++20 concepts for integrands/domains
├── policies.hpp          # Integration with boost::math::policies
├── qmc.hpp              # Quasi-Monte Carlo integration
├── regions.hpp          # Domain types (hypercube, simplex)
├── simplex.hpp          # Simplex-specific integration
├── sparse_grid.hpp      # Smolyak sparse grid methods
├── transforms.hpp       # Domain/coordinate transforms
├── workspace.hpp        # Memory management utilities
└── detail/
    ├── adaptivity.hpp    # Adaptive subdivision logic
    ├── cc_rules.hpp     # Clenshaw-Curtis rules
    ├── genz_malik_9_7_tables.hpp  # Generated tables
    ├── gm_rules.hpp     # Genz-Malik rule implementations
    ├── gp_rules.hpp     # Gauss-Patterson rules
    ├── sobol_owen.hpp   # Sobol sequences
    └── vectorize.hpp    # Vector integrand support
```

### Key Implementation Details

1. **Table-Driven Architecture**: The Genz-Malik rules use constexpr tables for compile-time optimization
2. **Dimension Support**: Currently targeting D=2-15 with table-driven implementations
3. **Precision Agnostic**: Supports multiple precision types including `boost::multiprecision`
4. **Header-Only**: Following Boost.Math conventions, all implementations are header-only

## Development Patterns

### When modifying gm_rules.hpp:

#### Mandatory Requirements:
1. **Namespace hierarchy**: All code MUST be within `namespace boost::math::cubature::detail::gm`
2. **Table-driven implementation**: Use `gm7<D>` and `gm9<D>` groups from `shaiko/genz_malik_9_7_tables.hpp`
3. **Shared radii**: Node generation MUST use `gm_lambda0`, `gm_lambda1`, `gm_lambda2`, `gm_lambda3`
4. **Scaling**: Apply `s = 2^(-D)` scaling for [0,1]^D consistently
5. **Exclusions**: NEVER include axis_lp_null (zero-weight) orbits in node enumeration

#### Debug Validation:
```cpp
#ifdef BOOST_MATH_DEBUG_CUBATURE
template<std::size_t D>
inline void gm_validate_weight_sums() {
    // Validate raw sum over [-1,1]^D equals 2^D
    // Validate scaled sum over [0,1]^D equals 1.0
    // Assert on failure with detailed error message
}
#endif
```

#### Node Enumeration Order (MUST match table group order):
1. Center (G1)
2. Axis λ₁ (G2)
3. Axis λ₂ (G3) 
4. Axis λ₃ (G4, deg-9 only)
5. Pair λ₁,λ₁ (G6)
6. Pair λ₁,λ₂ (G7, deg-9 only)
7. Triple λ₁,λ₁,λ₁ (G8, deg-9 only, D≥3)
8. Full diagonal λ₀ (G9)

### Testing Strategy:

#### Required Test Coverage:
1. **Weight sum validation**: All dimensions D=2-15 must sum to exactly 1.0
2. **Polynomial exactness**: Degree-7 rules exact for degree≤7, degree-9 for degree≤9
3. **Compilation**: Zero warnings with `-Wall -Wextra -Werror -pedantic`
4. **Precision types**: Test with float, double, long double, cpp_dec_float_50
5. **Cross-validation**: Compare against DCUHRE reference values

#### Test Execution Order:
```bash
# After any rule modification:
make clean
make gm97_weight_sums_test && ./test/gm97_weight_sums_test
make gm97_exactness_test && ./test/gm97_exactness_test
make gm_rules_exactness_test && ./test/gm_rules_exactness_test

# Full validation:
ctest --output-on-failure
```

## Common Issues and Solutions

### Issue: "expected unqualified-id" errors
**Cause**: Stray code blocks outside class/namespace scope
**Fix**: Check for duplicate code sections and ensure proper scoping

### Issue: "not in a namespace enclosing 'gm'"
**Cause**: Template specializations declared outside namespace
**Fix**: Ensure all specializations are inside namespace gm

### Issue: Unused constant warnings
**Cause**: Legacy literal constants not removed after table conversion
**Fix**: Remove unused l0, l1, l2, l3 constants in dimension-specific implementations

## Important Files to Reference

### Core Documentation (MUST READ)
- `shaiko/TECHNICAL_BLUEPRINT.md` - Engineering blueprint with milestones and acceptance criteria
- `shaiko/ARCHITECTURE.md` - Detailed architecture, scope, and API design
- `shaiko/API.md` - Public API specifications with usage examples
- `shaiko/ALGORITHMS.md` - Mathematical algorithms and implementation details
- `shaiko/TESTING_VALIDATION.md` - Test requirements and validation procedures
- `shaiko/REFERENCES.md` - Academic references and citations

### Genz-Malik Specific
- **`shaiko/genz_malik_9_7_tables.hpp`** - AUTHORITATIVE SOURCE for GM 9/7 tables
- `shaiko/9-7_Updates/GM_9-7_2D-5D.md` - Weight/node specifications for D=2-5
- `shaiko/9-7_Updates/GM_9-7_CSV_README.md` - CSV format documentation
- `shaiko/9-7_Updates/GM_9-7_VALIDATION_INDEX.md` - Validation data index
- `shaiko/genz_malik_*d_*.md` - Dimension-specific rule documentation

### Implementation Files
- `include/boost/math/cubature/detail/gm_rules.hpp` - Main rule implementation
- `include/boost/math/cubature/adaptive.hpp` - Adaptive cubature driver
- `include/boost/math/cubature/workspace.hpp` - Memory management
- `include/boost/math/cubature/regions.hpp` - Domain types

## Compiler Requirements

- C++14 compliant compiler minimum
- Tested compilers:
  - GCC 5+
  - Clang 5+
  - Visual Studio 2015+
- Standalone mode available with `BOOST_MATH_STANDALONE` macro

## Boost.Math Code Style & Best Practices

### Core Principles
- **Header-only design**: All implementations must be in headers
- **Policy-based design**: Integrate with `boost::math::policies` consistently
- **Compile-time computation**: Use constexpr and template metaprogramming extensively
- **Precision agnostic**: Support float, double, long double, and multiprecision types
- **Deterministic behavior**: Results must be reproducible regardless of thread count

### Template & Namespace Conventions
```cpp
// Correct Boost.Math namespace hierarchy
namespace boost { namespace math { namespace cubature { namespace detail {
namespace gm {  // Sub-namespace for Genz-Malik specifics

// Use trailing return types for complex template functions
template <class Real, std::size_t D>
constexpr auto compute_weights() -> std::array<Real, node_count<D>()> {
    // Implementation
}

}}}}}} // Close all namespaces properly
```

### Error Handling & Policies
- Use Boost.Math error policies, not exceptions directly
- Support user-defined error handling via policy templates
- Include comprehensive static_asserts for compile-time validation

### Performance Considerations
- Minimize runtime allocations (prefer stack/compile-time arrays)
- Use expression templates where appropriate
- Use cache-friendly data layouts when possible
- Cache frequently computed values at compile-time

### Testing Requirements
- Every rule must pass exactness tests for its stated polynomial degree
- Weight sums must equal 1.0 (or 2^D for [-1,1]^D) to machine precision
- Include both debug and release mode validation
- Test with multiple precision types (float, double, cpp_dec_float)

### Documentation Standards
- Use Doxygen-style comments for public APIs
- Include mathematical references in comments (cite papers/algorithms)
- Provide usage examples in header comments
- Document complexity guarantees and precision requirements

### Specific to Cubature Module
- All rules must be affine-invariant
- Support both scalar and vector-valued integrands
- Maintain separation between algorithm and data (tables)
- Use compile-time dispatch for dimension-specific optimizations
- Validate against reference implementations (DCUHRE, Cuba, cubature)