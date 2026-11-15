#include "fem1d/assembler.hpp"
#include "test_helpers.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <numeric>
#include <vector>

using Catch::Matchers::WithinAbs;
using fem1d::Assembler;
using fem1d::BarProperties;
using fem1d::constant;
using fem1d::Mesh;

TEST_CASE("Two unit elements give the classic 1D stiffness matrix and load", "[assembler]") {
    BarProperties bar;
    bar.body_force = constant(1.0);

    // h = 0.5, EA / h = 2, f A h / 2 = 0.25
    const auto system = Assembler::assemble(Mesh::uniform(1.0, 2), bar);
    const auto& k = system.matrix;

    REQUIRE(k.size() == 3);
    CHECK(k.diagonal() == std::vector<double>{2.0, 4.0, 2.0});
    CHECK(k.upper() == std::vector<double>{-2.0, -2.0});
    CHECK(k.lower() == std::vector<double>{-2.0, -2.0});
    CHECK(system.rhs == std::vector<double>{0.25, 0.5, 0.25});
}

TEST_CASE("The assembled stiffness is symmetric and singular before boundary conditions",
          "[assembler]") {
    BarProperties bar;
    bar.youngs_modulus = [](double x) { return x < 0.5 ? 1.0 : 3.0; };
    bar.area = [](double x) { return 1.0 + 0.5 * x; };

    const auto system = Assembler::assemble(Mesh::uniform(1.0, 7), bar);
    const auto& k = system.matrix;

    CHECK(k.lower() == k.upper());

    // A rigid-body translation u = const must produce no internal force.
    const std::vector<double> ones(k.size(), 1.0);
    const std::vector<double> zeros(k.size(), 0.0);
    CHECK(fem1d::test::max_abs_diff(k.multiply(ones), zeros) < 1e-12);
}

TEST_CASE("The total assembled load equals f A L for constant data", "[assembler]") {
    BarProperties bar;
    bar.area = constant(2.0);
    bar.body_force = constant(1.5);

    const auto system = Assembler::assemble(Mesh::uniform(3.0, 11), bar);
    const double total = std::accumulate(system.rhs.begin(), system.rhs.end(), 0.0);

    CHECK_THAT(total, WithinAbs(1.5 * 2.0 * 3.0, 1e-12));
}

TEST_CASE("Assembly works on a non-uniform mesh", "[assembler]") {
    // Element lengths 0.25 and 0.75 -> stiffnesses 4 and 4/3.
    const auto system = Assembler::assemble(Mesh({0.0, 0.25, 1.0}), BarProperties{});
    const auto& k = system.matrix;

    CHECK_THAT(k(0, 0), WithinAbs(4.0, 1e-14));
    CHECK_THAT(k(1, 1), WithinAbs(4.0 + 4.0 / 3.0, 1e-14));
    CHECK_THAT(k(2, 2), WithinAbs(4.0 / 3.0, 1e-14));
    CHECK_THAT(k(1, 2), WithinAbs(-4.0 / 3.0, 1e-14));
}
