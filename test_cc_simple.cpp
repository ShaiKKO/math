
#include <iostream>
#include <vector>
#include "./include/boost/math/cubature/detail/cc_rules.hpp"

using namespace boost::math::cubature::detail;

int main() {
  std::cout << "Testing CC rules directly..." << std::endl;
  
  auto nodes = clenshaw_curtis<double>::get_nodes(3);
  auto weights = clenshaw_curtis<double>::get_weights(3);
  
  std::cout << "Level 3: " << nodes.size() << " nodes" << std::endl;
  for (size_t i = 0; i < nodes.size(); ++i) {
    std::cout << "  Node " << i << ": " << nodes[i] << ", weight: " << weights[i] << std::endl;
  }
  
  return 0;
}

