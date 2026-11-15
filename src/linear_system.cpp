#include "fem1d/linear_system.hpp"

#include "stub.hpp"

namespace fem1d {

TridiagonalMatrix::TridiagonalMatrix(std::size_t /*n*/) {
    detail::not_implemented();
}

std::size_t TridiagonalMatrix::size() const noexcept {
    return diagonal_.size();
}

double TridiagonalMatrix::operator()(std::size_t /*row*/, std::size_t /*col*/) const {
    detail::not_implemented();
}

double& TridiagonalMatrix::at(std::size_t /*row*/, std::size_t /*col*/) {
    detail::not_implemented();
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

std::vector<double> TridiagonalMatrix::multiply(const std::vector<double>& /*x*/) const {
    detail::not_implemented();
}

} // namespace fem1d
