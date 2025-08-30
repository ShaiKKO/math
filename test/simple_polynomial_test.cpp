// Test simple polynomial integration 
// To verify vector integration is working

#include <boost/math/cubature/adaptive.hpp>
#include <iostream>
#include <cmath>

using namespace boost::math::cubature;

int main() {
  // Very simple test: constant function [1, 1, 1]
  auto f_const = [](const double* x, double* out, std::size_t m) {
    for (std::size_t i = 0; i < m; ++i) {
      out[i] = 1.0;
    }
  };
  
  hypercube<double> box(2);
  box.lower = {0, 0};
  box.upper = {1, 1};
  
  std::cout << "Test 1: Constant function [1,1,1] over unit square" << std::endl;
  auto results = integrate_adaptive_vector<double>(
    f_const, box, 3, 1e-10, 1e-10, 100);
  
  for (std::size_t i = 0; i < 3; ++i) {
    std::cout << "Component " << i << ": value=" << results[i].value 
              << " (expected 1.0), error=" << results[i].error << std::endl;
  }
  
  // Test with scalar version for comparison
  std::cout << "\nComparison with scalar version:" << std::endl;
  auto f_scalar = [](const double* x, std::size_t) -> double {
    return 1.0;
  };
  
  auto scalar_result = integrate_adaptive<double>(f_scalar, box, 1e-10, 1e-10);
  std::cout << "Scalar result: value=" << scalar_result.value 
            << ", error=" << scalar_result.error 
            << ", status=" << static_cast<int>(scalar_result.status) << std::endl;
  
  return 0;
}