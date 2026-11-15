#include "fem1d/linear_solver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace fem1d {

namespace {

void check_sizes(const TridiagonalMatrix& a, const std::vector<double>& b) {
    if (b.size() != a.size()) {
        throw std::invalid_argument("LinearSolver: right-hand side size does not match the matrix");
    }
}

/// Pivots smaller than this are treated as zero. Elimination on an n x n
/// matrix accumulates round-off of about n * eps * max|a_ij|, so a pivot
/// below that size carries no reliable information.
double pivot_tolerance(const TridiagonalMatrix& a) {
    double scale = 0.0;
    for (const auto* band : {&a.lower(), &a.diagonal(), &a.upper()}) {
        for (double value : *band) {
            scale = std::max(scale, std::abs(value));
        }
    }
    return static_cast<double>(a.size()) * std::numeric_limits<double>::epsilon() * scale;
}

void check_pivot(double pivot, double tolerance, std::size_t row) {
    // Written as !(|p| > tol) so that a NaN pivot is rejected too.
    if (!(std::abs(pivot) > tolerance)) {
        throw SingularMatrixError("LinearSolver: matrix is singular (zero pivot in row " +
                                  std::to_string(row) + ")");
    }
}

/// Weak diagonal dominance, with a little slack for round-off in the sums.
void check_diagonally_dominant(const TridiagonalMatrix& a) {
    const auto& lower = a.lower();
    const auto& diagonal = a.diagonal();
    const auto& upper = a.upper();
    const std::size_t n = a.size();
    constexpr double slack = 1.0 - 1e-12;

    for (std::size_t i = 0; i < n; ++i) {
        double off_diagonal = 0.0;
        if (i > 0) {
            off_diagonal += std::abs(lower[i - 1]);
        }
        if (i + 1 < n) {
            off_diagonal += std::abs(upper[i]);
        }
        if (std::abs(diagonal[i]) < slack * off_diagonal) {
            throw std::domain_error("ThomasSolver: row " + std::to_string(i) +
                                    " is not diagonally dominant; use TridiagonalLuSolver");
        }
    }
}

} // namespace

std::vector<double> ThomasSolver::solve(const TridiagonalMatrix& a,
                                        const std::vector<double>& b) const {
    check_sizes(a, b);
    check_diagonally_dominant(a);
    const std::size_t n = a.size();
    const auto& lower = a.lower();
    const auto& diagonal = a.diagonal();
    const auto& upper = a.upper();
    const double tolerance = pivot_tolerance(a);

    // Forward elimination: scale each row so the pivot becomes 1. c holds the
    // modified super-diagonal; x first holds the modified right-hand side.
    std::vector<double> c(n - 1);
    std::vector<double> x(n);

    double pivot = diagonal[0];
    check_pivot(pivot, tolerance, 0);
    if (n > 1) {
        c[0] = upper[0] / pivot;
    }
    x[0] = b[0] / pivot;

    for (std::size_t i = 1; i < n; ++i) {
        pivot = diagonal[i] - lower[i - 1] * c[i - 1];
        check_pivot(pivot, tolerance, i);
        if (i + 1 < n) {
            c[i] = upper[i] / pivot;
        }
        x[i] = (b[i] - lower[i - 1] * x[i - 1]) / pivot;
    }

    // Back substitution.
    for (std::size_t i = n - 1; i-- > 0;) {
        x[i] -= c[i] * x[i + 1];
    }
    return x;
}

std::vector<double> TridiagonalLuSolver::solve(const TridiagonalMatrix& a,
                                               const std::vector<double>& b) const {
    check_sizes(a, b);
    const std::size_t n = a.size();
    const double tolerance = pivot_tolerance(a);

    // Working copies. A row exchange at step i moves the old a(i+1, i+2) into
    // row i, so the upper triangle gains a second super-diagonal `upper2`.
    std::vector<double> lower = a.lower();
    std::vector<double> diagonal = a.diagonal();
    std::vector<double> upper = a.upper();
    std::vector<double> upper2(n > 2 ? n - 2 : 0, 0.0);
    std::vector<double> x = b;

    for (std::size_t i = 0; i + 1 < n; ++i) {
        if (std::abs(diagonal[i]) >= std::abs(lower[i])) {
            // No exchange: eliminate a(i+1, i) with row i.
            check_pivot(diagonal[i], tolerance, i);
            const double factor = lower[i] / diagonal[i];
            diagonal[i + 1] -= factor * upper[i];
            x[i + 1] -= factor * x[i];
        } else {
            // Exchange rows i and i+1, then eliminate with the new row i.
            //   row i   : [d_i,  u_i,     0      ]      [l_i, d_i+1, u_i+1]
            //   row i+1 : [l_i,  d_i+1,   u_i+1  ]  ->  [0,   u_i - f d_i+1, -f u_i+1]
            const double factor = diagonal[i] / lower[i];
            diagonal[i] = lower[i];
            const double old_diagonal = diagonal[i + 1];
            diagonal[i + 1] = upper[i] - factor * old_diagonal;
            if (i + 2 < n) {
                upper2[i] = upper[i + 1];
                upper[i + 1] = -factor * upper2[i];
            }
            upper[i] = old_diagonal;
            const double old_rhs = x[i];
            x[i] = x[i + 1];
            x[i + 1] = old_rhs - factor * x[i + 1];
        }
    }
    check_pivot(diagonal[n - 1], tolerance, n - 1);

    // Back substitution with the upper triangle (three diagonals).
    for (std::size_t i = n; i-- > 0;) {
        double sum = x[i];
        if (i + 1 < n) {
            sum -= upper[i] * x[i + 1];
        }
        if (i + 2 < n) {
            sum -= upper2[i] * x[i + 2];
        }
        x[i] = sum / diagonal[i];
    }
    return x;
}

std::vector<double> DenseLuSolver::solve(const TridiagonalMatrix& a,
                                         const std::vector<double>& b) const {
    check_sizes(a, b);
    const std::size_t n = a.size();
    const double tolerance = pivot_tolerance(a);

    // Dense row-major copy: m[i * n + j] = a(i, j).
    std::vector<double> m(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = (i > 0 ? i - 1 : 0); j < std::min(n, i + 2); ++j) {
            m[i * n + j] = a(i, j);
        }
    }
    std::vector<double> x = b;

    // Gaussian elimination with partial pivoting, i.e. P A = L U applied to b.
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(m[i * n + k]) > std::abs(m[p * n + k])) {
                p = i;
            }
        }
        check_pivot(m[p * n + k], tolerance, k);
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(m[k * n + j], m[p * n + j]);
            }
            std::swap(x[k], x[p]);
        }

        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = m[i * n + k] / m[k * n + k];
            for (std::size_t j = k + 1; j < n; ++j) {
                m[i * n + j] -= factor * m[k * n + j];
            }
            x[i] -= factor * x[k];
        }
    }

    // Back substitution with the upper triangle.
    for (std::size_t i = n; i-- > 0;) {
        double sum = x[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= m[i * n + j] * x[j];
        }
        x[i] = sum / m[i * n + i];
    }
    return x;
}

} // namespace fem1d
