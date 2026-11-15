// Cross-language check against the Python reference solver.
//
// tests/data holds the output of pinn-vs-fem-1d-bar's `run_fem` command for
// the `hetero` case (see tests/data/README.md and scripts/generate_golden.py).
// That case has no simple closed form and its FEM solution is not exact, so
// agreement to round-off means both codes assemble and solve the same system.

#include "cli/cases.hpp"
#include "fem1d/bar_solver.hpp"
#include "fem1d/io.hpp"
#include "test_helpers.hpp"

#include <catch2/catch_test_macros.hpp>
#include <string>

using fem1d::BarSolver;
using fem1d::BoundaryCondition;
using fem1d::Dirichlet;
using fem1d::Mesh;
using fem1d::Neumann;
using fem1d::Robin;
using fem1d::test::max_abs_diff;

namespace {

fem1d::Solution golden(const std::string& name) {
    return fem1d::read_csv(std::string(FEM1D_TEST_DATA_DIR) + "/" + name + ".csv");
}

void check_matches_golden(const std::string& name, std::size_t num_elements,
                          const BoundaryCondition& left, const BoundaryCondition& right) {
    const auto reference = golden(name);
    const auto mesh = Mesh::uniform(1.0, num_elements);
    const auto result =
        BarSolver{}.solve(mesh, fem1d::cli::reference_case("hetero").bar, left, right);

    REQUIRE(reference.x.size() == num_elements + 1);
    CHECK(max_abs_diff(result.x, reference.x) < 1e-15);
    CHECK(max_abs_diff(result.u, reference.u) < 1e-13);
}

} // namespace

SCENARIO("The C++ solver reproduces the Python reference solver", "[golden]") {
    GIVEN("the heterogeneous bar: two-layer E, tapered A, sine load") {
        WHEN("it is clamped at x = 0 and free at x = L, with 40 elements") {
            THEN("the displacements match the Python output to round-off") {
                check_matches_golden("hetero_neumann", 40, Dirichlet(0.0), Neumann(0.0));
            }
        }

        WHEN("it rests on a spring 10 u(L) + EA u'(L) = 0, with 40 elements") {
            THEN("the displacements match the Python output to round-off") {
                check_matches_golden("hetero_robin", 40, Dirichlet(0.0), Robin(10.0, 1.0, 0.0));
            }
        }

        WHEN("it is pulled at x = 0 and held at u(L) = 0.2, with 200 elements") {
            THEN("the displacements match the Python output to round-off") {
                check_matches_golden("hetero_mixed_fine", 200, Neumann(0.5), Dirichlet(0.2));
            }
        }
    }
}
