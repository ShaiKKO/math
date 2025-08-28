#include <boost/math/cubature/detail/reliability_metrics.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>

using namespace boost::math::cubature::detail;

template <typename Real>
void test_convergence_history() {
  std::cout << "Testing convergence_history..." << std::endl;
  
  convergence_history<Real> history;
  
  // Simulate a converging sequence
  history.record(Real(1.0), Real(0.1), 100, 1);
  history.record(Real(0.99), Real(0.05), 200, 2);
  history.record(Real(0.985), Real(0.025), 400, 4);
  history.record(Real(0.982), Real(0.012), 800, 8);
  history.record(Real(0.981), Real(0.006), 1600, 16);
  
  // Compute convergence rate (should be negative for convergence)
  Real rate = history.compute_convergence_rate();
  std::cout << "  Convergence rate: " << rate << std::endl;
  assert(rate < 0 && "Convergence rate should be negative");
  
  // Check monotonicity
  bool monotone = history.is_monotone();
  std::cout << "  Is monotone: " << (monotone ? "yes" : "no") << std::endl;
  assert(monotone && "Error should decrease monotonically");
  
  // Compute error ratio
  Real ratio = history.compute_error_ratio();
  std::cout << "  Error ratio (final/initial): " << ratio << std::endl;
  assert(ratio < 1 && "Error should decrease");
  
  // Estimate condition
  Real condition = history.estimate_condition();
  std::cout << "  Condition estimate: " << condition << std::endl;
  assert(condition >= 1 && condition <= 10 && "Condition should be reasonable");
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_chi2_reliability() {
  std::cout << "Testing chi2_reliability..." << std::endl;
  
  // Test case 1: Uniform errors (good reliability)
  {
    std::vector<Real> uniform_errors = {
      Real(0.01), Real(0.012), Real(0.009), Real(0.011), 
      Real(0.010), Real(0.008), Real(0.013), Real(0.011)
    };
    Real total_error = Real(0.084);
    
    Real chi2_prob = chi2_reliability<Real>::compute(uniform_errors, total_error);
    std::cout << "  Uniform errors chi2 probability: " << chi2_prob << std::endl;
    assert(chi2_prob > 0.3 && "Uniform errors should have good reliability");
  }
  
  // Test case 2: One outlier (poor reliability)
  {
    std::vector<Real> outlier_errors = {
      Real(0.001), Real(0.001), Real(0.001), Real(0.001),
      Real(0.1), Real(0.001), Real(0.001), Real(0.001)  // One large outlier
    };
    Real total_error = Real(0.107);
    
    Real chi2_prob = chi2_reliability<Real>::compute(outlier_errors, total_error);
    std::cout << "  Outlier errors chi2 probability: " << chi2_prob << std::endl;
    assert(chi2_prob < 0.5 && "Outlier errors should have poor reliability");
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_reliability_calculator() {
  std::cout << "Testing reliability_calculator..." << std::endl;
  
  // Test good convergence scenario
  {
    Real chi2_prob = Real(0.8);          // Good chi-squared
    Real convergence_rate = Real(-1.5);  // Good convergence
    Real error_ratio = Real(0.01);       // Significant improvement
    bool is_monotone = true;             // Monotonic
    Real condition_est = Real(2.0);      // Well-conditioned
    
    Real reliability = reliability_calculator<Real>::compute_reliability_factor(
        chi2_prob, convergence_rate, error_ratio, is_monotone, condition_est);
    
    std::cout << "  Good convergence reliability: " << reliability << std::endl;
    assert(reliability > 0.7 && "Good convergence should have high reliability");
  }
  
  // Test poor convergence scenario
  {
    Real chi2_prob = Real(0.1);          // Poor chi-squared
    Real convergence_rate = Real(-0.2);  // Slow convergence
    Real error_ratio = Real(0.8);        // Little improvement
    bool is_monotone = false;            // Non-monotonic
    Real condition_est = Real(8.0);      // Ill-conditioned
    
    Real reliability = reliability_calculator<Real>::compute_reliability_factor(
        chi2_prob, convergence_rate, error_ratio, is_monotone, condition_est);
    
    std::cout << "  Poor convergence reliability: " << reliability << std::endl;
    assert(reliability < 0.4 && "Poor convergence should have low reliability");
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_full_metrics_computation() {
  std::cout << "Testing full metrics computation..." << std::endl;
  
  // Setup convergence history
  convergence_history<Real> history;
  history.record(Real(1.0), Real(0.1), 100, 1);
  history.record(Real(0.99), Real(0.05), 250, 3);
  history.record(Real(0.985), Real(0.025), 500, 7);
  history.record(Real(0.982), Real(0.012), 1000, 15);
  history.record(Real(0.981), Real(0.006), 2000, 31);
  
  // Setup region errors
  std::vector<Real> region_errors;
  for (int i = 0; i < 31; ++i) {
    // Generate somewhat uniform errors with slight variation
    Real base_error = Real(0.006) / 31;
    Real variation = Real(1.0) + Real(0.2) * std::sin(Real(i));
    region_errors.push_back(base_error * variation);
  }
  
  Real final_error = Real(0.006);
  std::size_t max_depth = 5;
  
  // Compute full metrics
  auto metrics = reliability_calculator<Real>::compute_metrics(
      history, region_errors, final_error, max_depth);
  
  // Display results
  std::cout << "  Reliability Metrics:" << std::endl;
  std::cout << "    Chi2 probability: " << metrics.chi2_probability << std::endl;
  std::cout << "    Convergence rate: " << metrics.convergence_rate << std::endl;
  std::cout << "    Error ratio: " << metrics.error_ratio << std::endl;
  std::cout << "    Reliability factor: " << metrics.reliability_factor << std::endl;
  std::cout << "    Refinement depth: " << metrics.refinement_depth << std::endl;
  std::cout << "    Regions processed: " << metrics.regions_processed << std::endl;
  std::cout << "    Monotone convergence: " << (metrics.monotone_convergence ? "yes" : "no") << std::endl;
  std::cout << "    Condition estimate: " << metrics.condition_estimate << std::endl;
  std::cout << "    Is reliable: " << (metrics.is_reliable() ? "yes" : "no") << std::endl;
  std::cout << "    Convergence quality: " << metrics.convergence_quality() << std::endl;
  
  // Verify results are reasonable
  assert(metrics.chi2_probability >= 0 && metrics.chi2_probability <= 1);
  assert(metrics.reliability_factor >= 0 && metrics.reliability_factor <= 1);
  assert(metrics.refinement_depth == max_depth);
  assert(metrics.regions_processed == 31);
  
  std::cout << "  PASSED" << std::endl;
}

int main() {
  std::cout << "\n===== Testing Reliability Metrics =====\n" << std::endl;
  
  using Real = double;
  
  test_convergence_history<Real>();
  test_chi2_reliability<Real>();
  test_reliability_calculator<Real>();
  test_full_metrics_computation<Real>();
  
  std::cout << "\n===== All Tests PASSED =====\n" << std::endl;
  
  return 0;
}