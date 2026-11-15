#include "fem1d/linear_solver.hpp"

#include "stub.hpp"

namespace fem1d {

std::vector<double> ThomasSolver::solve(const TridiagonalMatrix& /*a*/,
                                        const std::vector<double>& /*b*/) const {
    detail::not_implemented();
}

std::vector<double> DenseLuSolver::solve(const TridiagonalMatrix& /*a*/,
                                         const std::vector<double>& /*b*/) const {
    detail::not_implemented();
}

} // namespace fem1d
