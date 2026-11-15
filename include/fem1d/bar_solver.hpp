#pragma once

#include "fem1d/bar_properties.hpp"
#include "fem1d/boundary_conditions.hpp"
#include "fem1d/linear_solver.hpp"
#include "fem1d/mesh.hpp"

#include <memory>
#include <vector>

namespace fem1d {

/// Nodal coordinates and nodal displacements.
struct Solution {
    std::vector<double> x;
    std::vector<double> u;
};

/// Solves -(E A u')' = f A on a mesh: assembles, applies both boundary
/// conditions and calls the linear solver it was given.
class BarSolver {
public:
    /// \throws std::invalid_argument if `linear_solver` is null.
    explicit BarSolver(
        std::unique_ptr<LinearSolver> linear_solver = std::make_unique<ThomasSolver>());

    [[nodiscard]] Solution solve(const Mesh& mesh, const BarProperties& bar,
                                 const BoundaryCondition& left,
                                 const BoundaryCondition& right) const;

private:
    std::unique_ptr<LinearSolver> linear_solver_;
};

} // namespace fem1d
