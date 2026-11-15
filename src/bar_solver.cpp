#include "fem1d/bar_solver.hpp"

#include "fem1d/assembler.hpp"

#include <stdexcept>
#include <utility>

namespace fem1d {

BarSolver::BarSolver(std::unique_ptr<LinearSolver> linear_solver)
    : linear_solver_(std::move(linear_solver)) {
    if (!linear_solver_) {
        throw std::invalid_argument("BarSolver: linear solver must not be null");
    }
}

Solution BarSolver::solve(const Mesh& mesh, const BarProperties& bar, const BoundaryCondition& left,
                          const BoundaryCondition& right) const {
    auto system = Assembler::assemble(mesh, bar);

    // The order does not matter: natural conditions only add to the diagonal
    // and right-hand side of their own end node, while a Dirichlet condition
    // reads and clears only the off-diagonal entries next to its node.
    left.apply(system, End::Left);
    right.apply(system, End::Right);

    return Solution{mesh.nodes(), linear_solver_->solve(system.matrix, system.rhs)};
}

} // namespace fem1d
