// How each boundary condition changes the assembled system, checked on the
// matrix and right-hand side directly (no solve).

#include "fem1d/assembler.hpp"
#include "fem1d/boundary_conditions.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <stdexcept>
#include <vector>

using Catch::Matchers::WithinAbs;
using fem1d::Assembler;
using fem1d::BarProperties;
using fem1d::constant;
using fem1d::Dirichlet;
using fem1d::End;
using fem1d::LinearSystem;
using fem1d::Mesh;
using fem1d::Neumann;
using fem1d::Robin;

namespace {

// Three elements of length 1/3 with EA = 1 and f A = 3:
// K = 3 * tridiag(-1, 2, -1) with corner entries 3, F = (0.5, 1, 1, 0.5).
LinearSystem three_element_system() {
    BarProperties bar;
    bar.body_force = constant(3.0);
    return Assembler::assemble(Mesh::uniform(1.0, 3), bar);
}

} // namespace

TEST_CASE("Dirichlet replaces the end equation and keeps the matrix symmetric", "[bc]") {
    auto system = three_element_system();
    Dirichlet(2.0).apply(system, End::Left);

    const auto& k = system.matrix;
    CHECK(k(0, 0) == 1.0);
    CHECK(k(0, 1) == 0.0);
    CHECK(k(1, 0) == 0.0);
    CHECK(system.rhs[0] == 2.0);
    // The known value moves to the neighbour's right-hand side: 1 - (-3)(2) = 7.
    CHECK_THAT(system.rhs[1], WithinAbs(7.0, 1e-13));
    CHECK(k.lower() == k.upper());
}

TEST_CASE("Dirichlet at the right end acts on the last node", "[bc]") {
    auto system = three_element_system();
    Dirichlet(-1.0).apply(system, End::Right);

    const auto& k = system.matrix;
    CHECK(k(3, 3) == 1.0);
    CHECK(k(3, 2) == 0.0);
    CHECK(k(2, 3) == 0.0);
    CHECK(system.rhs[3] == -1.0);
    // 1 - (-3)(-1) = -2
    CHECK_THAT(system.rhs[2], WithinAbs(-2.0, 1e-13));
}

TEST_CASE("Neumann adds the end force with the sign of the outward normal", "[bc]") {
    auto system = three_element_system();
    const auto k_before = system.matrix;

    Neumann(0.6).apply(system, End::Right);
    Neumann(0.4).apply(system, End::Left);

    CHECK_THAT(system.rhs[3], WithinAbs(0.5 + 0.6, 1e-15));
    CHECK_THAT(system.rhs[0], WithinAbs(0.5 - 0.4, 1e-15));
    CHECK(system.matrix.diagonal() == k_before.diagonal());
}

TEST_CASE("Robin adds alpha/beta to the end stiffness and g/beta to the load", "[bc]") {
    auto system = three_element_system();
    Robin(2.0, 4.0, 1.0).apply(system, End::Right);

    CHECK_THAT(system.matrix(3, 3), WithinAbs(3.0 + 0.5, 1e-14));
    CHECK_THAT(system.rhs[3], WithinAbs(0.5 + 0.25, 1e-15));
}

TEST_CASE("Robin at the left end enters with a negative sign", "[bc]") {
    auto system = three_element_system();
    Robin(2.0, 4.0, 1.0).apply(system, End::Left);

    CHECK_THAT(system.matrix(0, 0), WithinAbs(3.0 - 0.5, 1e-14));
    CHECK_THAT(system.rhs[0], WithinAbs(0.5 - 0.25, 1e-15));
}

TEST_CASE("Robin with beta = 0 behaves like Dirichlet u = g / alpha", "[bc]") {
    auto robin_system = three_element_system();
    auto dirichlet_system = three_element_system();

    Robin(2.0, 0.0, 1.0).apply(robin_system, End::Right);
    Dirichlet(0.5).apply(dirichlet_system, End::Right);

    CHECK(robin_system.matrix.diagonal() == dirichlet_system.matrix.diagonal());
    CHECK(robin_system.matrix.upper() == dirichlet_system.matrix.upper());
    CHECK(robin_system.rhs == dirichlet_system.rhs);
}

TEST_CASE("Robin with alpha = beta = 0 is not a condition", "[bc]") {
    CHECK_THROWS_AS(Robin(0.0, 0.0, 1.0), std::invalid_argument);
}

TEST_CASE("The order in which the two ends are applied does not matter", "[bc]") {
    const Dirichlet dirichlet(0.3);
    const Neumann neumann(0.7);
    const Robin robin(2.0, 1.0, 0.5);
    const Robin robin_as_dirichlet(4.0, 0.0, 0.6);
    const std::vector<const fem1d::BoundaryCondition*> conditions{&dirichlet, &neumann, &robin,
                                                                  &robin_as_dirichlet};

    // On one element both ends touch the same off-diagonal entries.
    for (const std::size_t elements : {std::size_t{1}, std::size_t{3}}) {
        BarProperties bar;
        bar.body_force = constant(3.0);
        const auto assembled = Assembler::assemble(Mesh::uniform(1.0, elements), bar);

        for (const auto* left : conditions) {
            for (const auto* right : conditions) {
                auto left_first = assembled;
                left->apply(left_first, End::Left);
                right->apply(left_first, End::Right);

                auto right_first = assembled;
                right->apply(right_first, End::Right);
                left->apply(right_first, End::Left);

                // Same entries; F may differ by round-off, since floating-point
                // additions in a different order are not exactly equal.
                CHECK(left_first.matrix.diagonal() == right_first.matrix.diagonal());
                CHECK(left_first.matrix.lower() == right_first.matrix.lower());
                CHECK(left_first.matrix.upper() == right_first.matrix.upper());
                for (std::size_t i = 0; i < assembled.rhs.size(); ++i) {
                    CHECK_THAT(left_first.rhs[i], WithinAbs(right_first.rhs[i], 1e-14));
                }
            }
        }
    }
}
