#include "fem1d/element.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <stdexcept>

using Catch::Matchers::WithinAbs;
using fem1d::BarProperties;
using fem1d::constant;
using fem1d::LinearBarElement;

TEST_CASE("Element stiffness is (EA/h) [[1, -1], [-1, 1]] with EA taken at the midpoint",
          "[element]") {
    BarProperties bar;
    bar.youngs_modulus = [](double x) { return 2.0 + x; };
    bar.area = constant(3.0);

    // h = 0.5, midpoint 0.25: EA = 2.25 * 3 = 6.75, EA / h = 13.5
    const LinearBarElement element(0.0, 0.5);
    const auto m = element.compute(bar);

    CHECK_THAT(m.stiffness[0][0], WithinAbs(13.5, 1e-13));
    CHECK_THAT(m.stiffness[0][1], WithinAbs(-13.5, 1e-13));
    CHECK_THAT(m.stiffness[1][0], WithinAbs(-13.5, 1e-13));
    CHECK_THAT(m.stiffness[1][1], WithinAbs(13.5, 1e-13));
}

TEST_CASE("Element load splits f A h equally between its two nodes", "[element]") {
    BarProperties bar;
    bar.area = constant(3.0);
    bar.body_force = [](double x) { return x; };

    // midpoint 0.25: f A = 0.75, total load f A h = 0.375
    const LinearBarElement element(0.0, 0.5);
    const auto m = element.compute(bar);

    CHECK_THAT(m.load[0], WithinAbs(0.1875, 1e-15));
    CHECK_THAT(m.load[1], WithinAbs(0.1875, 1e-15));
}

TEST_CASE("A rigid-body translation produces no element force", "[element]") {
    BarProperties bar;
    bar.youngs_modulus = [](double x) { return 1.0 + 10.0 * x * x; };

    const auto k = LinearBarElement(0.3, 0.45).compute(bar).stiffness;

    CHECK_THAT(k[0][0] + k[0][1], WithinAbs(0.0, 1e-12));
    CHECK_THAT(k[1][0] + k[1][1], WithinAbs(0.0, 1e-12));
}

TEST_CASE("An element needs a positive length", "[element]") {
    CHECK_THROWS_AS(LinearBarElement(1.0, 1.0), std::invalid_argument);
    CHECK_THROWS_AS(LinearBarElement(1.0, 0.5), std::invalid_argument);
}
