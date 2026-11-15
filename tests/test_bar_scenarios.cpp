// Behaviour of the whole solver, written as BDD scenarios. Each closed-form
// case is ported from tests/test_fem_analytic.py in the Python reference
// repository (E = A = 1, L = 1, 40 elements). Linear elements reproduce these
// solutions exactly at the nodes, so the tolerance is round-off.

#include "fem1d/bar_solver.hpp"
#include "test_helpers.hpp"

#include <catch2/catch_test_macros.hpp>
#include <memory>
#include <stdexcept>
#include <utility>

using fem1d::BarProperties;
using fem1d::BarSolver;
using fem1d::BoundaryCondition;
using fem1d::constant;
using fem1d::DenseLuSolver;
using fem1d::Dirichlet;
using fem1d::Mesh;
using fem1d::Neumann;
using fem1d::Robin;
using fem1d::ScalarField;
using fem1d::SingularMatrixError;
using fem1d::Solution;
using fem1d::ThomasSolver;
using fem1d::TridiagonalLuSolver;
using fem1d::test::max_nodal_error;

namespace {

constexpr double tolerance = 1e-12;

Solution solve_unit_bar(const ScalarField& body_force, const BoundaryCondition& left,
                        const BoundaryCondition& right) {
    BarProperties bar;
    bar.body_force = body_force;
    return BarSolver{}.solve(Mesh::uniform(1.0, 40), bar, left, right);
}

} // namespace

SCENARIO("A bar clamped at x = 0 and pulled at the tip", "[scenario]") {
    GIVEN("u(0) = 0, EA u'(1) = P = 0.6 and no body force") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(0.0), Dirichlet(0.0), Neumann(0.6));

            THEN("the displacement grows linearly: u = P x") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return 0.6 * x; }) < tolerance);
            }
        }
    }
}

SCENARIO("A bar clamped at x = 0 under its own uniform weight", "[scenario]") {
    GIVEN("-u'' = 1, u(0) = 0 and a free end u'(1) = 0") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(1.0), Dirichlet(0.0), Neumann(0.0));

            THEN("u = x - x^2 / 2") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return x - x * x / 2.0; }) <
                      tolerance);
            }
        }
    }
}

SCENARIO("A bar fixed at both ends under a uniform load", "[scenario]") {
    GIVEN("-u'' = 1 and u(0) = u(1) = 0") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(1.0), Dirichlet(0.0), Dirichlet(0.0));

            THEN("u = x (1 - x) / 2") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return x * (1.0 - x) / 2.0; }) <
                      tolerance);
            }
        }
    }
}

SCENARIO("The force is applied at the left end instead", "[scenario]") {
    GIVEN("EA u'(0) = P = 0.6 and u(1) = 0") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(0.0), Neumann(0.6), Dirichlet(0.0));

            THEN("u = P (x - 1)") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return 0.6 * (x - 1.0); }) <
                      tolerance);
            }
        }
    }
}

SCENARIO("A bar supported by a spring", "[scenario]") {
    GIVEN("u(0) = 0 and a spring at the tip: 2 u(1) + EA u'(1) = 1") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(0.0), Dirichlet(0.0), Robin(2.0, 1.0, 1.0));

            THEN("u = x / 3") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return x / 3.0; }) < tolerance);
            }
        }
    }

    GIVEN("the Robin condition at the left end: 2 u(0) + EA u'(0) = 1, u(1) = 0") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(0.0), Robin(2.0, 1.0, 1.0), Dirichlet(0.0));

            THEN("u = 1 - x") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return 1.0 - x; }) < tolerance);
            }
        }
    }

    GIVEN("beta different from one: 2 u(1) + 2 EA u'(1) = 1, u(0) = 0") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(0.0), Dirichlet(0.0), Robin(2.0, 2.0, 1.0));

            THEN("u = x / 4") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return x / 4.0; }) < tolerance);
            }
        }
    }

    GIVEN("beta = 0, so the condition becomes 2 u(1) = 1, and u(0) = 0") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(0.0), Dirichlet(0.0), Robin(2.0, 0.0, 1.0));

            THEN("u = x / 2") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return x / 2.0; }) < tolerance);
            }
        }
    }
}

SCENARIO("Prescribed displacements that are not zero", "[scenario]") {
    GIVEN("u(0) = 1 and u(1) = 0") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(0.0), Dirichlet(1.0), Dirichlet(0.0));

            THEN("u = 1 - x") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return 1.0 - x; }) < tolerance);
            }
        }
    }

    GIVEN("u(0) = 1 and a tip force EA u'(1) = 0.6") {
        WHEN("the bar is solved") {
            const auto s = solve_unit_bar(constant(0.0), Dirichlet(1.0), Neumann(0.6));

            THEN("u = 1 + P x") {
                CHECK(max_nodal_error(s.x, s.u, [](double x) { return 1.0 + 0.6 * x; }) <
                      tolerance);
            }
        }
    }
}

SCENARIO("A bar that is free at both ends", "[scenario]") {
    GIVEN("forces at both ends and no support") {
        WHEN("the bar is solved") {
            THEN("the solver reports a singular system instead of returning garbage") {
                CHECK_THROWS_AS(solve_unit_bar(constant(0.0), Neumann(-1.0), Neumann(1.0)),
                                SingularMatrixError);
            }
        }
    }
}

SCENARIO("The linear solver is a replaceable strategy", "[scenario]") {
    GIVEN("a heterogeneous bar with a spring support") {
        BarProperties bar;
        bar.youngs_modulus = [](double x) { return x < 0.5 ? 1.0 : 3.0; };
        bar.area = [](double x) { return 1.0 + 0.5 * x; };
        bar.body_force = [](double x) { return 2.0 * x; };
        const auto mesh = Mesh::uniform(1.0, 40);

        WHEN("it is solved with each of the three linear solvers") {
            const auto solve_with = [&](std::unique_ptr<fem1d::LinearSolver> solver) {
                return BarSolver{std::move(solver)}.solve(mesh, bar, Dirichlet(0.0),
                                                          Robin(10.0, 1.0, 0.0));
            };
            const auto thomas = solve_with(std::make_unique<ThomasSolver>());
            const auto tridiagonal_lu = solve_with(std::make_unique<TridiagonalLuSolver>());
            const auto dense = solve_with(std::make_unique<DenseLuSolver>());

            THEN("all give the same displacements") {
                CHECK(fem1d::test::max_abs_diff(thomas.u, dense.u) < tolerance);
                CHECK(fem1d::test::max_abs_diff(tridiagonal_lu.u, dense.u) < tolerance);
                CHECK(thomas.x == mesh.nodes());
            }
        }
    }
}

TEST_CASE("A bar solver needs a linear solver", "[scenario]") {
    CHECK_THROWS_AS(BarSolver{nullptr}, std::invalid_argument);
}
