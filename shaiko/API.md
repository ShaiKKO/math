# Boost.Math::cubature — API (draft)

## Headers

```cpp
#include <boost/math/cubature/adaptive.hpp>
#include <boost/math/cubature/sparse_grid.hpp>
#include <boost/math/cubature/qmc.hpp>
#include <boost/math/cubature/simplex.hpp>
#include <boost/math/cubature/transforms.hpp>
```

## Concepts (informal)

```cpp
template <class F, class Real>
concept ScalarIntegrand = requires(F f, std::span<const Real> x) {
  { f(x) } -> std::convertible_to<Real>;
};

template <class F, class Real>
concept VectorIntegrand = requires(F f, std::span<const Real> x, Real* out, std::size_t m) {
  { f(x, out, m) } -> std::same_as<void>;
};
```

## Domain types

```cpp
template <class Real> struct hypercube { std::vector<Real> lower, upper; };
template <class Real> struct simplex  { std::vector<std::vector<Real>> vertices; /* d+1 vertices, each of length d */ };
```

## Result type

```cpp
enum class status_code { success, maxeval_reached, cancelled };
template <class Real>
struct result {
  Real value;
  Real error;
  std::uint64_t evaluations;
  status_code status;
};
```

## Adaptive

```cpp
template <class Real, class F, class Policy = boost::math::policies::policy<>>
  requires (ScalarIntegrand<F, Real> || VectorIntegrand<F, Real>)
auto integrate_adaptive(const F& f, const hypercube<Real>& box,
                        Real abs_tol, Real rel_tol,
                        std::uint64_t max_eval = 0,
                        Policy const& pol = Policy{}) -> result<Real>;
```

## Sparse‑grid (Smolyak)

```cpp
template <class Real, class F, class Policy = boost::math::policies::policy<>>
auto integrate_sparse_grid(const F& f, const hypercube<Real>& box,
                           unsigned level,
                           Policy const& pol = Policy{}) -> result<Real>;
```

## QMC / RQMC

```cpp
template <class Real, class F, class Policy = boost::math::policies::policy<>, class Sobol = boost::random::sobol_engine<>>
auto integrate_qmc(const F& f, const hypercube<Real>& box,
                   std::size_t n_points_per_rep,
                   unsigned replicates = 1, bool tent_transform = true,
                   Policy const& pol = Policy{}) -> result<Real>;
```

Note: For integrate_qmc, result<Real>::error reports the replicate-based standard error (std error) across the k replicates.


## Simplex

```cpp
template <class Real, class F, class Policy = boost::math::policies::policy<>>
auto integrate_simplex(const F& f, const simplex<Real>& tri_or_simplex,
                       Real abs_tol, Real rel_tol,
                       Policy const& pol = Policy{}) -> result<Real>;
```

## Examples

### 1) Scalar integrand (adaptive)

```cpp
using Real = double;
auto f = [](std::span<const Real> x){ return std::exp(-(x[0]*x[0] + x[1]*x[1])); };
boost::math::cubature::hypercube<Real> box{{-5,-5},{5,5}};
auto res = integrate_adaptive<Real>(f, box, 1e-10, 1e-10);
```

### 2) Vector integrand (QMC)

```cpp
auto g = [](std::span<const double> x, double* out, std::size_t m){
  out[0] = std::sin(x[0]) * std::sin(x[1]);
  if(m>1) out[1] = std::cos(x[0] + x[1]);
};
hypercube<double> box{{0,0},{M_PI,M_PI}};
auto r = integrate_qmc<double>(g, box, 1<<18, /*rep*/ 8, /*tent*/ true);
```

### 3) Simplex (triangle via Duffy)

```cpp
auto tri = simplex<double>{/* three vertices */};
auto F = [](std::span<const double> x){ return x[0] + x[1]; };
auto r = integrate_simplex(F, tri, 1e-8, 1e-8);
```

**Notes**: Policies and multiprecision types are supported by replacing `double` with `boost::multiprecision::cpp_dec_float_50`. citeturn6search1
