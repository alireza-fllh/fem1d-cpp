#include "fem1d/mesh.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <stdexcept>
#include <vector>

using Catch::Matchers::WithinAbs;
using fem1d::Mesh;

TEST_CASE("A uniform mesh has equally spaced nodes on [0, L]", "[mesh]") {
    const auto mesh = Mesh::uniform(2.0, 4);

    REQUIRE(mesh.num_elements() == 4);
    REQUIRE(mesh.num_nodes() == 5);

    const std::vector<double> expected{0.0, 0.5, 1.0, 1.5, 2.0};
    for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK_THAT(mesh.nodes()[i], WithinAbs(expected[i], 1e-15));
    }
    CHECK_THAT(mesh.element_length(1), WithinAbs(0.5, 1e-15));
    CHECK_THAT(mesh.element_midpoint(1), WithinAbs(0.75, 1e-15));
}

TEST_CASE("The last node of a uniform mesh is exactly L", "[mesh]") {
    const auto mesh = Mesh::uniform(1.0, 3);
    REQUIRE(mesh.nodes().front() == 0.0);
    REQUIRE(mesh.nodes().back() == 1.0);
}

TEST_CASE("A mesh can have non-uniform nodes", "[mesh]") {
    const Mesh mesh({0.0, 0.1, 0.4, 1.0});

    REQUIRE(mesh.num_elements() == 3);
    CHECK_THAT(mesh.element_length(0), WithinAbs(0.1, 1e-15));
    CHECK_THAT(mesh.element_length(2), WithinAbs(0.6, 1e-15));
    CHECK_THAT(mesh.element_midpoint(2), WithinAbs(0.7, 1e-15));
}

TEST_CASE("A mesh rejects invalid input", "[mesh]") {
    CHECK_THROWS_AS(Mesh::uniform(1.0, 0), std::invalid_argument);
    CHECK_THROWS_AS(Mesh::uniform(0.0, 4), std::invalid_argument);
    CHECK_THROWS_AS(Mesh::uniform(-1.0, 4), std::invalid_argument);
    CHECK_THROWS_AS(Mesh(std::vector<double>{0.0}), std::invalid_argument);
    CHECK_THROWS_AS(Mesh({0.0, 0.5, 0.5, 1.0}), std::invalid_argument);
    CHECK_THROWS_AS(Mesh({0.0, 1.0, 0.5}), std::invalid_argument);
}

TEST_CASE("Asking for an element that does not exist throws", "[mesh]") {
    const auto mesh = Mesh::uniform(1.0, 2);
    CHECK_THROWS_AS(mesh.element_length(2), std::out_of_range);
    CHECK_THROWS_AS(mesh.element_midpoint(2), std::out_of_range);
}
