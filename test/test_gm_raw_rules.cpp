// Test that raw rules are working correctly

#include <boost/math/cubature/detail/gm_rules.hpp>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace boost::math::cubature::detail::gm;

int main() {
  std::cout << "Testing GM 9/7 raw rules for 2D\n" << std::endl;
  
  // Get nodes and weights
  auto nodes_9 = raw_rule_fam<9, 2, family_9_7>::nodes_with_zero<double>();
  auto weights_9 = raw_rule_fam<9, 2, family_9_7>::weights_with_zero<double>();
  auto weights_7 = raw_rule_fam<7, 2, family_9_7>::weights_with_zero<double>();
  
  std::cout << "Number of nodes: " << nodes_9.size() << std::endl;
  std::cout << "Degree-9 weights size: " << weights_9.size() << std::endl;
  std::cout << "Degree-7 weights size: " << weights_7.size() << std::endl;
  
  // Test integration of constant function (should be 1.0 for unit square)
  double sum_9 = 0, sum_7 = 0;
  for (size_t i = 0; i < nodes_9.size(); ++i) {
    sum_9 += weights_9[i];
    sum_7 += weights_7[i];
  }
  
  std::cout << "\nWeight sums (should be 1.0):" << std::endl;
  std::cout << "Degree-9: " << sum_9 << std::endl;
  std::cout << "Degree-7: " << sum_7 << std::endl;
  
  // Test integration of x^2
  double integral_9 = 0, integral_7 = 0;
  for (size_t i = 0; i < nodes_9.size(); ++i) {
    double x = nodes_9[i][0];
    integral_9 += x * x * weights_9[i];
    integral_7 += x * x * weights_7[i];
  }
  
  std::cout << "\nIntegral of x^2 (should be 1/3 = 0.333...):" << std::endl;
  std::cout << "Degree-9: " << integral_9 << std::endl;
  std::cout << "Degree-7: " << integral_7 << std::endl;
  
  // Show first few nodes and weights
  std::cout << "\nFirst 10 nodes and weights:" << std::endl;
  std::cout << std::setw(5) << "i" 
            << std::setw(15) << "x" 
            << std::setw(15) << "y" 
            << std::setw(15) << "w9" 
            << std::setw(15) << "w7" << std::endl;
  
  for (size_t i = 0; i < std::min(size_t(10), nodes_9.size()); ++i) {
    std::cout << std::setw(5) << i
              << std::setw(15) << nodes_9[i][0]
              << std::setw(15) << nodes_9[i][1]
              << std::setw(15) << weights_9[i]
              << std::setw(15) << weights_7[i] << std::endl;
  }
  
  return 0;
}