#include "fem1d/linear_solver.hpp"
#include "test_helpers.hpp"

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <random>
#include <stdexcept>
#include <vector>

using fem1d::DenseLuSolver;
using fem1d::SingularMatrixError;
using fem1d::ThomasSolver;
using fem1d::TridiagonalMatrix;
using fem1d::test::max_abs_diff;

namespace {

TridiagonalMatrix second_difference(std::size_t n) {
    TridiagonalMatrix a(n);
    for (std::size_t i = 0; i < n; ++i) {
        a.at(i, i) = 2.0;
        if (i + 1 < n) {
            a.at(i, i + 1) = -1.0;
            a.at(i + 1, i) = -1.0;
        }
    }
    return a;
}

TridiagonalMatrix random_diagonally_dominant(std::size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> off(-1.0, 1.0);
    TridiagonalMatrix a(n);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        a.at(i, i + 1) = off(rng);
        a.at(i + 1, i) = off(rng);
    }
    for (std::size_t i = 0; i < n; ++i) {
        a.at(i, i) = 3.0 + off(rng);
    }
    return a;
}

} // namespace

TEMPLATE_TEST_CASE("Solvers recover a known solution of tridiag(-1, 2, -1)", "[solver]",
                   ThomasSolver, DenseLuSolver) {
    // A x = b with x = (1, 2, 3, 4, 5)
    const std::vector<double> b{0.0, 0.0, 0.0, 0.0, 6.0};
    const std::vector<double> expected{1.0, 2.0, 3.0, 4.0, 5.0};

    const auto x = TestType{}.solve(second_difference(5), b);

    CHECK(max_abs_diff(x, expected) < 1e-13);
}

TEMPLATE_TEST_CASE("Solvers reject a singular matrix", "[solver]", ThomasSolver, DenseLuSolver) {
    // [[1, -1], [-1, 1]]: the free-free bar of one element.
    TridiagonalMatrix a(2);
    a.at(0, 0) = 1.0;
    a.at(0, 1) = -1.0;
    a.at(1, 0) = -1.0;
    a.at(1, 1) = 1.0;

    CHECK_THROWS_AS(TestType{}.solve(a, {1.0, -1.0}), SingularMatrixError);
}

TEMPLATE_TEST_CASE("Solvers reject a right-hand side of the wrong size", "[solver]", ThomasSolver,
                   DenseLuSolver) {
    CHECK_THROWS_AS(TestType{}.solve(second_difference(3), {1.0, 2.0}), std::invalid_argument);
}

TEST_CASE("Thomas and dense LU agree on a random diagonally dominant system", "[solver]") {
    const std::size_t n = 50;
    const auto a = random_diagonally_dominant(n, 42U);

    std::mt19937 rng(7U);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<double> b(n);
    for (auto& value : b) {
        value = dist(rng);
    }

    const auto x_thomas = ThomasSolver{}.solve(a, b);
    const auto x_dense = DenseLuSolver{}.solve(a, b);

    CHECK(max_abs_diff(x_thomas, x_dense) < 1e-12);
    CHECK(max_abs_diff(a.multiply(x_thomas), b) < 1e-12);
}

TEST_CASE("Only the pivoting solver handles a zero leading diagonal entry", "[solver]") {
    // [[0, 1], [1, 1]] x = (1, 2) has the solution x = (1, 1), but Thomas
    // does not pivot and stops at the zero in the top-left corner.
    TridiagonalMatrix a(2);
    a.at(0, 1) = 1.0;
    a.at(1, 0) = 1.0;
    a.at(1, 1) = 1.0;
    const std::vector<double> b{1.0, 2.0};

    CHECK_THROWS_AS(ThomasSolver{}.solve(a, b), SingularMatrixError);
    CHECK(max_abs_diff(DenseLuSolver{}.solve(a, b), {1.0, 1.0}) < 1e-14);
}

TEMPLATE_TEST_CASE("Solvers handle a 1 x 1 system", "[solver]", ThomasSolver, DenseLuSolver) {
    TridiagonalMatrix a(1);
    a.at(0, 0) = 4.0;

    CHECK(TestType{}.solve(a, {2.0}) == std::vector<double>{0.5});
}
