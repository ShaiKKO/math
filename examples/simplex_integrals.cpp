#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <boost/math/cubature/cubature.hpp>
#include <boost/math/cubature/simplex.hpp>

using namespace boost::math::cubature;

// Example 1: Triangle integration
void example_triangle() {
    std::cout << "=== Example 1: Triangle Integration ===" << std::endl;
    
    using Real = double;
    
    // Triangle with vertices (0,0), (2,0), (1,1)
    std::array<std::array<Real, 2>, 3> vertices = {{
        {{0, 0}},
        {{2, 0}},
        {{1, 1}}
    }};
    
    // Test different functions
    
    // 1. Constant function
    auto constant = [](const Real* x) -> Real {
        return Real(1);
    };
    
    auto result1 = integrate_triangle<Real>(constant, vertices, 1e-10, 1e-10);
    std::cout << "Constant function (area): " << result1.value << std::endl;
    std::cout << "  Expected: 1.0" << std::endl;
    
    // 2. Linear function f(x,y) = x + y
    auto linear = [](const Real* x) -> Real {
        return x[0] + x[1];
    };
    
    auto result2 = integrate_triangle<Real>(linear, vertices, 1e-10, 1e-10);
    std::cout << "Linear function ∫(x+y): " << result2.value << std::endl;
    
    // 3. Quadratic function f(x,y) = x² + y²
    auto quadratic = [](const Real* x) -> Real {
        return x[0]*x[0] + x[1]*x[1];
    };
    
    auto result3 = integrate_triangle<Real>(quadratic, vertices, 1e-10, 1e-10);
    std::cout << "Quadratic function ∫(x²+y²): " << result3.value << std::endl;
    
    std::cout << std::endl;
}

// Example 2: Tetrahedron integration
void example_tetrahedron() {
    std::cout << "=== Example 2: Tetrahedron Integration ===" << std::endl;
    
    using Real = double;
    
    // Standard tetrahedron with vertices at origin and unit vectors
    std::array<std::array<Real, 3>, 4> vertices = {{
        {{0, 0, 0}},
        {{1, 0, 0}},
        {{0, 1, 0}},
        {{0, 0, 1}}
    }};
    
    // 1. Volume calculation
    auto constant = [](const Real* x) -> Real {
        return Real(1);
    };
    
    auto result1 = integrate_tetrahedron<Real>(constant, vertices, 1e-10, 1e-10);
    std::cout << "Volume of standard tetrahedron: " << result1.value << std::endl;
    std::cout << "  Expected: " << Real(1)/Real(6) << std::endl;
    
    // 2. Center of mass calculation
    auto x_coord = [](const Real* x) -> Real {
        return x[0];
    };
    
    auto result2 = integrate_tetrahedron<Real>(x_coord, vertices, 1e-10, 1e-10);
    Real centroid_x = result2.value / (Real(1)/Real(6)); // Divide by volume
    std::cout << "Centroid x-coordinate: " << centroid_x << std::endl;
    std::cout << "  Expected: 0.25" << std::endl;
    
    // 3. Moment of inertia
    auto moment = [](const Real* x) -> Real {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    };
    
    auto result3 = integrate_tetrahedron<Real>(moment, vertices, 1e-10, 1e-10);
    std::cout << "Second moment ∫(x²+y²+z²): " << result3.value << std::endl;
    
    std::cout << std::endl;
}

// Example 3: General simplex in higher dimensions
void example_general_simplex() {
    std::cout << "=== Example 3: 4D Simplex Integration ===" << std::endl;
    
    using Real = double;
    
    // 4D standard simplex
    simplex<Real> simp;
    simp.vertices = {
        {0, 0, 0, 0},
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    
    // Volume of 4D simplex
    auto constant = [](const Real* x) -> Real {
        return Real(1);
    };
    
    auto result = integrate_simplex<Real>(constant, simp, 1e-10, 1e-10);
    Real expected_volume = Real(1) / Real(24); // Volume = 1/4!
    
    std::cout << "Volume of 4D standard simplex: " << result.value << std::endl;
    std::cout << "  Expected: " << expected_volume << std::endl;
    std::cout << "  Error: " << std::abs(result.value - expected_volume) << std::endl;
    std::cout << "  Evaluations: " << result.evaluations << std::endl;
    
    std::cout << std::endl;
}

// Example 4: Rotated and translated simplex
void example_transformed_simplex() {
    std::cout << "=== Example 4: Transformed Triangle ===" << std::endl;
    
    using Real = double;
    
    // Triangle rotated and translated
    std::array<std::array<Real, 2>, 3> vertices = {{
        {{1, 1}},
        {{3, 2}},
        {{2, 4}}
    }};
    
    // Integrate a Gaussian-like function
    auto gaussian = [](const Real* x) -> Real {
        Real cx = Real(2), cy = Real(2.5); // Center
        Real dx = x[0] - cx, dy = x[1] - cy;
        return std::exp(-(dx*dx + dy*dy));
    };
    
    auto result = integrate_triangle<Real>(gaussian, vertices, 1e-8, 1e-8);
    
    std::cout << "Integral of Gaussian over triangle: " << result.value << std::endl;
    std::cout << "  Error estimate: " << result.error << std::endl;
    std::cout << "  Evaluations: " << result.evaluations << std::endl;
    
    std::cout << std::endl;
}

// Example 5: Physical application - Mass of triangular plate
void example_physical_application() {
    std::cout << "=== Example 5: Mass of Non-uniform Triangular Plate ===" << std::endl;
    
    using Real = double;
    
    // Triangular plate with vertices
    std::array<std::array<Real, 2>, 3> vertices = {{
        {{0, 0}},
        {{4, 0}},
        {{2, 3}}
    }};
    
    // Density function that varies with position
    // ρ(x,y) = 1 + 0.1*x + 0.2*y (kg/m²)
    auto density = [](const Real* x) -> Real {
        return Real(1) + Real(0.1) * x[0] + Real(0.2) * x[1];
    };
    
    auto mass_result = integrate_triangle<Real>(density, vertices, 1e-10, 1e-10);
    
    std::cout << "Total mass of plate: " << mass_result.value << " kg" << std::endl;
    
    // Calculate center of mass
    auto x_moment = [](const Real* x) -> Real {
        Real rho = Real(1) + Real(0.1) * x[0] + Real(0.2) * x[1];
        return x[0] * rho;
    };
    
    auto y_moment = [](const Real* x) -> Real {
        Real rho = Real(1) + Real(0.1) * x[0] + Real(0.2) * x[1];
        return x[1] * rho;
    };
    
    auto mx = integrate_triangle<Real>(x_moment, vertices, 1e-10, 1e-10);
    auto my = integrate_triangle<Real>(y_moment, vertices, 1e-10, 1e-10);
    
    Real com_x = mx.value / mass_result.value;
    Real com_y = my.value / mass_result.value;
    
    std::cout << "Center of mass: (" << com_x << ", " << com_y << ")" << std::endl;
    
    std::cout << std::endl;
}

// Example 6: Comparing with adaptive integration on box
void example_comparison_with_box() {
    std::cout << "=== Example 6: Simplex vs Box Integration ===" << std::endl;
    
    using Real = double;
    
    // Right triangle with vertices (0,0), (1,0), (0,1)
    std::array<std::array<Real, 2>, 3> vertices = {{
        {{0, 0}},
        {{1, 0}},
        {{0, 1}}
    }};
    
    // Function to integrate
    auto f = [](const Real* x) -> Real {
        return std::sin(x[0]) * std::cos(x[1]);
    };
    
    // Method 1: Direct simplex integration
    auto simplex_result = integrate_triangle<Real>(f, vertices, 1e-8, 1e-8);
    std::cout << "Simplex integration:" << std::endl;
    std::cout << "  Value: " << simplex_result.value << std::endl;
    std::cout << "  Evaluations: " << simplex_result.evaluations << std::endl;
    
    // Method 2: Box integration with boundary check
    auto f_box = [&f](const Real* x) -> Real {
        // Check if point is inside triangle
        if (x[0] >= 0 && x[1] >= 0 && x[0] + x[1] <= 1) {
            return f(x);
        }
        return Real(0);
    };
    
    hypercube<Real> box(2);
    box.lower = {0, 0};
    box.upper = {1, 1};
    
    auto box_result = integrate_adaptive<Real>(f_box, box, 1e-8, 1e-8);
    std::cout << "Box integration with boundary:" << std::endl;
    std::cout << "  Value: " << box_result.value << std::endl;
    std::cout << "  Evaluations: " << box_result.evaluations << std::endl;
    
    std::cout << "Difference: " << std::abs(simplex_result.value - box_result.value) << std::endl;
    std::cout << "\nNote: Simplex integration is more efficient for simplex domains!" << std::endl;
    
    std::cout << std::endl;
}

int main() {
    std::cout << std::fixed << std::setprecision(10);
    
    example_triangle();
    example_tetrahedron();
    example_general_simplex();
    example_transformed_simplex();
    example_physical_application();
    example_comparison_with_box();
    
    std::cout << "All simplex examples completed!" << std::endl;
    
    return 0;
}