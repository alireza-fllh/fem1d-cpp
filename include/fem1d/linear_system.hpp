#pragma once

#include <cstddef>
#include <vector>

namespace fem1d {

/// Square tridiagonal matrix stored as three diagonals.
///
/// Linear elements in 1D only couple neighbouring nodes, so the global
/// stiffness matrix is tridiagonal and needs O(n) storage instead of O(n^2).
class TridiagonalMatrix {
public:
    /// A zero matrix of size n x n.
    /// \throws std::invalid_argument if `n == 0`.
    explicit TridiagonalMatrix(std::size_t n);

    [[nodiscard]] std::size_t size() const noexcept;

    /// Entry (row, col); zero outside the band.
    /// \throws std::out_of_range if row or col is not below size().
    [[nodiscard]] double operator()(std::size_t row, std::size_t col) const;

    /// Writable entry (row, col) inside the band.
    /// \throws std::out_of_range if |row - col| > 1 or an index is not below size().
    double& at(std::size_t row, std::size_t col);

    /// Main diagonal, size n.
    [[nodiscard]] const std::vector<double>& diagonal() const noexcept;
    /// Entries (i + 1, i), size n - 1.
    [[nodiscard]] const std::vector<double>& lower() const noexcept;
    /// Entries (i, i + 1), size n - 1.
    [[nodiscard]] const std::vector<double>& upper() const noexcept;

    /// Matrix-vector product.
    /// \throws std::invalid_argument if `x.size() != size()`.
    [[nodiscard]] std::vector<double> multiply(const std::vector<double>& x) const;

private:
    std::vector<double> lower_;
    std::vector<double> diagonal_;
    std::vector<double> upper_;
};

/// The global system K u = F.
struct LinearSystem {
    TridiagonalMatrix matrix;
    std::vector<double> rhs;
};

} // namespace fem1d
