#pragma once

#include "fem1d/linear_system.hpp"

#include <stdexcept>
#include <vector>

namespace fem1d {

/// Raised when a solver meets a zero (or numerically zero) pivot.
class SingularMatrixError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Strategy for solving A x = b with a tridiagonal A.
class LinearSolver {
public:
    virtual ~LinearSolver() = default;

    /// \throws std::invalid_argument if `b.size() != a.size()`.
    /// \throws SingularMatrixError if a pivot is numerically zero.
    [[nodiscard]] virtual std::vector<double> solve(const TridiagonalMatrix& a,
                                                    const std::vector<double>& b) const = 0;

protected:
    LinearSolver() = default;
    LinearSolver(const LinearSolver&) = default;
    LinearSolver& operator=(const LinearSolver&) = default;
    LinearSolver(LinearSolver&&) = default;
    LinearSolver& operator=(LinearSolver&&) = default;
};

/// Thomas algorithm: Gaussian elimination specialised to tridiagonal
/// matrices, without pivoting. O(n) time and memory.
///
/// Without pivoting it is only stable when the matrix is diagonally dominant,
/// so it checks that first. Stiffness matrices of the bar with physical
/// supports (Dirichlet, Neumann, springs with positive stiffness) satisfy it.
class ThomasSolver final : public LinearSolver {
public:
    /// \throws std::domain_error if some row has |a_ii| < sum_{j != i} |a_ij|.
    [[nodiscard]] std::vector<double> solve(const TridiagonalMatrix& a,
                                            const std::vector<double>& b) const override;
};

/// Gaussian elimination with partial pivoting, specialised to tridiagonal
/// matrices (the algorithm of LAPACK's dgtsv). Row exchanges add one extra
/// upper diagonal, so it stays O(n) and works for any non-singular matrix.
class TridiagonalLuSolver final : public LinearSolver {
public:
    [[nodiscard]] std::vector<double> solve(const TridiagonalMatrix& a,
                                            const std::vector<double>& b) const override;
};

/// LU factorisation with partial pivoting on a dense copy of the matrix.
/// O(n^3) time and O(n^2) memory, the same as the Python reference solver's
/// `numpy.linalg.solve`. Kept as a reference and as a benchmark baseline.
class DenseLuSolver final : public LinearSolver {
public:
    [[nodiscard]] std::vector<double> solve(const TridiagonalMatrix& a,
                                            const std::vector<double>& b) const override;
};

} // namespace fem1d
