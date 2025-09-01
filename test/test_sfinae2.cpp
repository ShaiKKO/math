#include <iostream>
#include <vector>
#include <type_traits>

// Helper to detect valid expressions
template<typename... Ts> struct make_void { typedef void type; };
template<typename... Ts> using void_t = typename make_void<Ts...>::type;

// SFINAE helper for detecting integrand signature  
template <typename Real, typename Func, typename Point, typename = void>
struct integrand_traits {
  static constexpr int signature = 0;
};

// Specialization for f(vector) signature
template <typename Real, typename Func, typename Point>
struct integrand_traits<Real, Func, Point,
    void_t<decltype(Real(std::declval<const Func&>()(std::declval<const Point&>())))>> {
  static constexpr int signature = 1;
};

// Specialization for f(pointer, size) signature
template <typename Real, typename Func, typename Point>
struct integrand_traits<Real, Func, Point,
    void_t<decltype(Real(std::declval<const Func&>()(
        std::declval<const Real*>(), 
        std::declval<std::size_t>())))>> {
  static constexpr int signature = 2;
};

// Specialization for f(pointer) signature
template <typename Real, typename Func, typename Point>
struct integrand_traits<Real, Func, Point,
    void_t<decltype(Real(std::declval<const Func&>()(
        std::declval<const Real*>())))>> {
  static constexpr int signature = 3;
};

int main() {
  using Real = double;
  using Point = std::vector<Real>;
  
  auto f_vec = [](const std::vector<Real>& x) -> Real { return x[0] + x[1]; };
  auto f_ptr = [](const Real* x, std::size_t /*n*/) -> Real { return x[0] + x[1]; };
  auto f_legacy = [](const Real* x) -> Real { return x[0] + x[1]; };
  
  std::cout << "f_vec signature: " << integrand_traits<Real, decltype(f_vec), Point>::signature << std::endl;
  std::cout << "f_ptr signature: " << integrand_traits<Real, decltype(f_ptr), Point>::signature << std::endl;
  std::cout << "f_legacy signature: " << integrand_traits<Real, decltype(f_legacy), Point>::signature << std::endl;
  
  return 0;
}