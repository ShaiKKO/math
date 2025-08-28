
#include <iostream>
#include <vector>

int main() {
  // Test multi-index generation for 2D, level 3
  int dimension = 2;
  int level = 3;
  
  std::cout << "Multi-indices for d=" << dimension << ", l=" << level << std::endl;
  std::cout << "Range: " << (level - dimension + 1) << " <= |i| <= " << level << std::endl;
  
  // These should be included:
  // |i| = 2: (0,2), (1,1), (2,0)  
  // |i| = 3: (0,3), (1,2), (2,1), (3,0)
  
  for (int i0 = 0; i0 <= level; ++i0) {
    for (int i1 = 0; i1 <= level; ++i1) {
      int sum = i0 + i1;
      if (sum >= level - dimension + 1 && sum <= level) {
        int diff = level - sum;
        
        // Binomial coefficient C(d-1, diff) = C(1, diff)
        int binom = (diff == 0 || diff == 1) ? 1 : 0;
        
        // Sign
        int sign = (diff % 2 == 0) ? 1 : -1;
        
        int coeff = sign * binom;
        
        if (coeff != 0) {
          std::cout << "  [" << i0 << "," << i1 << "] sum=" << sum 
                    << " diff=" << diff << " coeff=" << coeff << std::endl;
        }
      }
    }
  }
  
  return 0;
}

