#include <boost/math/cubature/detail/vector_adapter.hpp>
#include <boost/math/cubature/detail/genz_malik_evaluator.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>

using namespace boost::math::cubature::detail;

// Test vector integrand: f(x) = [sin(x[0]), cos(x[1]), exp(-||x||^2)]
template <typename Real>
class test_vector_integrand {
public:
  void operator()(const Real* x, Real* out, std::size_t m) const {
    assert(m >= 3);
    out[0] = std::sin(x[0]);
    out[1] = std::cos(x[1]);
    
    // Compute ||x||^2 for third component
    Real norm_sq = x[0]*x[0] + x[1]*x[1];
    out[2] = std::exp(-norm_sq);
  }
};

template <typename Real>
void test_vector_result_aggregator() {
  std::cout << "Testing vector_result_aggregator..." << std::endl;
  
  const std::size_t num_components = 3;
  
  // Test L-infinity norm (default)
  {
    vector_result_aggregator<Real> agg(num_components, error_norm::l_infinity);
    
    std::vector<Real> values1 = {Real(1.0), Real(2.0), Real(3.0)};
    std::vector<Real> errors1 = {Real(0.01), Real(0.02), Real(0.03)};
    
    agg.add_results(values1, errors1);
    
    Real total_error = agg.get_total_error();
    std::cout << "  L-infinity error: " << total_error << std::endl;
    assert(std::abs(total_error - Real(0.03)) < 1e-10);
    
    // Test subtraction
    agg.subtract_results(values1, errors1);
    total_error = agg.get_total_error();
    assert(std::abs(total_error) < 1e-10);
  }
  
  // Test L2 norm
  {
    vector_result_aggregator<Real> agg(num_components, error_norm::l2);
    
    std::vector<Real> values = {Real(1.0), Real(1.0), Real(1.0)};
    std::vector<Real> errors = {Real(0.01), Real(0.01), Real(0.01)};
    
    agg.add_results(values, errors);
    
    Real total_error = agg.get_total_error();
    Real expected = std::sqrt(Real(3)) * Real(0.01);
    std::cout << "  L2 error: " << total_error << " (expected: " << expected << ")" << std::endl;
    assert(std::abs(total_error - expected) < 1e-10);
  }
  
  // Test L1 norm
  {
    vector_result_aggregator<Real> agg(num_components, error_norm::l1);
    
    std::vector<Real> values = {Real(1.0), Real(1.0), Real(1.0)};
    std::vector<Real> errors = {Real(0.01), Real(0.02), Real(0.03)};
    
    agg.add_results(values, errors);
    
    Real total_error = agg.get_total_error();
    Real expected = Real(0.06);
    std::cout << "  L1 error: " << total_error << " (expected: " << expected << ")" << std::endl;
    assert(std::abs(total_error - expected) < 1e-10);
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_vector_workspace() {
  std::cout << "Testing vector_workspace..." << std::endl;
  
  const std::size_t num_components = 5;
  const std::size_t max_nodes = 100;
  
  vector_workspace<Real> workspace(num_components, max_nodes);
  
  // Test buffer access
  Real* buffer0 = workspace.get_values_buffer(0);
  Real* buffer1 = workspace.get_values_buffer(1);
  Real* scratch = workspace.get_scratch_buffer();
  
  // Buffers should not overlap
  assert(buffer1 - buffer0 == static_cast<std::ptrdiff_t>(num_components));
  assert(scratch != buffer0 && scratch != buffer1);
  
  // Test writing to buffers
  for (std::size_t i = 0; i < num_components; ++i) {
    buffer0[i] = static_cast<Real>(i);
    scratch[i] = static_cast<Real>(i * 2);
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_scalar_vector_adapters() {
  std::cout << "Testing scalar/vector adapters..." << std::endl;
  
  // Test scalar to vector adapter
  {
    auto scalar_f = [](const Real* x) { return x[0] * x[0]; };
    scalar_to_vector_adapter<Real, decltype(scalar_f)> adapter(scalar_f);
    
    Real x[] = {Real(2.0)};
    Real out[1];
    adapter(x, out, 1);
    
    assert(std::abs(out[0] - Real(4.0)) < 1e-10);
    std::cout << "  Scalar to vector: PASSED" << std::endl;
  }
  
  // Test vector component adapter
  {
    test_vector_integrand<Real> vector_f;
    vector_component_adapter<Real, test_vector_integrand<Real>> adapter(vector_f, 1, 3);
    
    Real x[] = {Real(0), Real(0)};
    Real result = adapter(x);
    
    assert(std::abs(result - Real(1.0)) < 1e-10);  // cos(0) = 1
    std::cout << "  Vector component: PASSED" << std::endl;
  }
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_vectorized_evaluator() {
  std::cout << "Testing vectorized_evaluator..." << std::endl;
  
  const std::size_t Dim = 2;
  const std::size_t num_components = 3;
  
  // Create test points
  std::vector<std::array<Real, Dim>> points = {
    {Real(0), Real(0)},
    {Real(1), Real(0)},
    {Real(0), Real(1)},
    {Real(1), Real(1)}
  };
  
  // Test integrand
  test_vector_integrand<Real> f;
  
  // Evaluate batch
  std::vector<std::vector<Real>> values_per_component;
  vectorized_evaluator<Real>::template evaluate_batch<test_vector_integrand<Real>, Dim>(
      f, points, num_components, values_per_component);
  
  // Verify dimensions
  assert(values_per_component.size() == num_components);
  assert(values_per_component[0].size() == points.size());
  
  // Verify some values
  // At point (0,0): [sin(0), cos(0), exp(0)] = [0, 1, 1]
  assert(std::abs(values_per_component[0][0] - Real(0)) < 1e-10);
  assert(std::abs(values_per_component[1][0] - Real(1)) < 1e-10);
  assert(std::abs(values_per_component[2][0] - Real(1)) < 1e-10);
  
  // Test weighted sum
  std::vector<Real> values = {Real(1), Real(2), Real(3), Real(4)};
  std::vector<Real> weights = {Real(0.25), Real(0.25), Real(0.25), Real(0.25)};
  Real sum = vectorized_evaluator<Real>::compute_weighted_sum(values, weights);
  assert(std::abs(sum - Real(2.5)) < 1e-10);
  
  std::cout << "  PASSED" << std::endl;
}

template <typename Real>
void test_vector_region() {
  std::cout << "Testing vector_region..." << std::endl;
  
  const std::size_t dim = 2;
  const std::size_t num_components = 3;
  
  vector_region<Real> reg(dim, num_components);
  
  // Set component errors
  reg.component_errors = {Real(0.01), Real(0.02), Real(0.03)};
  
  // Test L-infinity norm
  reg.compute_aggregated_error(error_norm::l_infinity);
  assert(std::abs(reg.aggregated_error - Real(0.03)) < 1e-10);
  assert(std::abs(reg.error - Real(0.03)) < 1e-10);  // Should be stored in base class
  
  // Test L2 norm
  reg.compute_aggregated_error(error_norm::l2);
  Real expected_l2 = std::sqrt(Real(0.01*0.01 + 0.02*0.02 + 0.03*0.03));
  assert(std::abs(reg.aggregated_error - expected_l2) < 1e-10);
  
  // Test L1 norm
  reg.compute_aggregated_error(error_norm::l1);
  assert(std::abs(reg.aggregated_error - Real(0.06)) < 1e-10);
  
  std::cout << "  PASSED" << std::endl;
}

int main() {
  std::cout << "\n===== Testing Vector Adapters =====\n" << std::endl;
  
  using Real = double;
  
  test_vector_result_aggregator<Real>();
  test_vector_workspace<Real>();
  test_scalar_vector_adapters<Real>();
  test_vectorized_evaluator<Real>();
  test_vector_region<Real>();
  
  std::cout << "\n===== All Tests PASSED =====\n" << std::endl;
  
  return 0;
}