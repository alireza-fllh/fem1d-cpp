// Mesh refinement study. Unlike the polynomial cases in test_bar_scenarios,
// a sine load is not reproduced exactly: the midpoint rule for the load
// introduces an O(h^2) error, which is the expected rate for linear elements.

#include "fem1d/bar_solver.hpp"
#include "test_helpers.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>
#include <vector>

using Catch::Matchers::WithinAbs;
using fem1d::BarProperties;
using fem1d::BarSolver;
using fem1d::Dirichlet;
using fem1d::Mesh;
using fem1d::Neumann;

SCENARIO("Linear elements converge at second order", "[convergence]") {
    const double pi = std::acos(-1.0);

    GIVEN("-u'' = 2 sin(2 pi x), u(0) = 0, u'(1) = 0") {
        BarProperties bar;
        bar.body_force = [pi](double x) { return 2.0 * std::sin(2.0 * pi * x); };
        const auto exact = [pi](double x) {
            return std::sin(2.0 * pi * x) / (2.0 * pi * pi) - x / pi;
        };

        WHEN("the mesh is refined from 10 to 160 elements") {
            std::vector<double> errors;
            for (std::size_t n : {10U, 20U, 40U, 80U, 160U}) {
                const auto s =
                    BarSolver{}.solve(Mesh::uniform(1.0, n), bar, Dirichlet(0.0), Neumann(0.0));
                errors.push_back(fem1d::test::max_nodal_error(s.x, s.u, exact));
            }

            THEN("each halving of h divides the nodal error by about four") {
                for (std::size_t i = 0; i + 1 < errors.size(); ++i) {
                    const double rate = std::log2(errors[i] / errors[i + 1]);
                    INFO("refinement " << i << ": error " << errors[i] << " -> " << errors[i + 1]);
                    CHECK_THAT(rate, WithinAbs(2.0, 0.1));
                }
            }
        }
    }
}
