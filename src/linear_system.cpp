#include "fem1d/linear_system.hpp"

#include <stdexcept>

namespace fem1d {

TridiagonalMatrix::TridiagonalMatrix(std::size_t n) {
    if (n == 0) {
        throw std::invalid_argument("TridiagonalMatrix: size must be at least 1");
    }
    lower_.assign(n - 1, 0.0);
    diagonal_.assign(n, 0.0);
    upper_.assign(n - 1, 0.0);
}

std::size_t TridiagonalMatrix::size() const noexcept {
    return diagonal_.size();
}

double TridiagonalMatrix::operator()(std::size_t row, std::size_t col) const {
    if (row >= size() || col >= size()) {
        throw std::out_of_range("TridiagonalMatrix: index out of range");
    }
    if (row == col) {
        return diagonal_[row];
    }
    if (col == row + 1) {
        return upper_[row];
    }
    if (row == col + 1) {
        return lower_[col];
    }
    return 0.0;
}

double& TridiagonalMatrix::at(std::size_t row, std::size_t col) {
    if (row >= size() || col >= size()) {
        throw std::out_of_range("TridiagonalMatrix: index out of range");
    }
    if (row == col) {
        return diagonal_[row];
    }
    if (col == row + 1) {
        return upper_[row];
    }
    if (row == col + 1) {
        return lower_[col];
    }
    throw std::out_of_range("TridiagonalMatrix: entry is outside the band");
}

const std::vector<double>& TridiagonalMatrix::diagonal() const noexcept {
    return diagonal_;
}

const std::vector<double>& TridiagonalMatrix::lower() const noexcept {
    return lower_;
}

const std::vector<double>& TridiagonalMatrix::upper() const noexcept {
    return upper_;
}

std::vector<double> TridiagonalMatrix::multiply(const std::vector<double>& x) const {
    const std::size_t n = size();
    if (x.size() != n) {
        throw std::invalid_argument("TridiagonalMatrix::multiply: size mismatch");
    }
    std::vector<double> y(n);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = diagonal_[i] * x[i];
        if (i > 0) {
            sum += lower_[i - 1] * x[i - 1];
        }
        if (i + 1 < n) {
            sum += upper_[i] * x[i + 1];
        }
        y[i] = sum;
    }
    return y;
}

} // namespace fem1d
