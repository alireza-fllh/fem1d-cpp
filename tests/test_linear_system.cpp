#include "fem1d/linear_system.hpp"

#include <catch2/catch_test_macros.hpp>
#include <stdexcept>
#include <vector>

using fem1d::TridiagonalMatrix;

TEST_CASE("A new tridiagonal matrix is zero", "[linear_system]") {
    const TridiagonalMatrix a(3);

    REQUIRE(a.size() == 3);
    REQUIRE(a.diagonal().size() == 3);
    REQUIRE(a.lower().size() == 2);
    REQUIRE(a.upper().size() == 2);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            CHECK(a(i, j) == 0.0);
        }
    }
}

TEST_CASE("Entries inside the band can be written and read back", "[linear_system]") {
    TridiagonalMatrix a(3);
    a.at(0, 0) = 1.0;
    a.at(0, 1) = 2.0;
    a.at(1, 0) = 3.0;
    a.at(2, 2) += 4.0;

    CHECK(a(0, 0) == 1.0);
    CHECK(a(0, 1) == 2.0);
    CHECK(a(1, 0) == 3.0);
    CHECK(a(2, 2) == 4.0);
    CHECK(a.upper()[0] == 2.0);
    CHECK(a.lower()[0] == 3.0);
}

TEST_CASE("Entries outside the band read as zero but cannot be written", "[linear_system]") {
    TridiagonalMatrix a(4);

    CHECK(a(0, 2) == 0.0);
    CHECK(a(3, 0) == 0.0);
    CHECK_THROWS_AS(a.at(0, 2), std::out_of_range);
    CHECK_THROWS_AS(a.at(3, 1), std::out_of_range);
}

TEST_CASE("Indices beyond the matrix size throw", "[linear_system]") {
    TridiagonalMatrix a(2);

    CHECK_THROWS_AS(a(2, 2), std::out_of_range);
    CHECK_THROWS_AS(a.at(2, 1), std::out_of_range);
    CHECK_THROWS_AS(a.at(1, 2), std::out_of_range);
}

TEST_CASE("A tridiagonal matrix needs at least one row", "[linear_system]") {
    CHECK_THROWS_AS(TridiagonalMatrix(0), std::invalid_argument);
}

TEST_CASE("Matrix-vector product with tridiag(-1, 2, -1)", "[linear_system]") {
    TridiagonalMatrix a(4);
    for (std::size_t i = 0; i < 4; ++i) {
        a.at(i, i) = 2.0;
        if (i + 1 < 4) {
            a.at(i, i + 1) = -1.0;
            a.at(i + 1, i) = -1.0;
        }
    }

    CHECK(a.multiply({1.0, 1.0, 1.0, 1.0}) == std::vector<double>{1.0, 0.0, 0.0, 1.0});
    CHECK(a.multiply({1.0, 2.0, 3.0, 4.0}) == std::vector<double>{0.0, 0.0, 0.0, 5.0});
    CHECK_THROWS_AS(a.multiply({1.0, 2.0}), std::invalid_argument);
}
