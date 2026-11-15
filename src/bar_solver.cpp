#include "fem1d/bar_solver.hpp"

#include "fem1d/assembler.hpp"
#include "stub.hpp"

#include <utility>

namespace fem1d {

BarSolver::BarSolver(std::unique_ptr<LinearSolver> linear_solver)
    : linear_solver_(std::move(linear_solver)) {
    detail::not_implemented();
}

Solution BarSolver::solve(const Mesh& /*mesh*/, const BarProperties& /*bar*/,
                          const BoundaryCondition& /*left*/,
                          const BoundaryCondition& /*right*/) const {
    detail::not_implemented();
}

} // namespace fem1d
