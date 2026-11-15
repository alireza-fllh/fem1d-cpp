#pragma once

#include "fem1d/linear_system.hpp"

namespace fem1d {

/// The two ends of the bar: x = 0 (Left) and x = L (Right).
enum class End { Left, Right };

/// A condition imposed on the assembled system at one end of the bar.
///
/// Integrating by parts leaves the boundary term n EA u'(x_b) v(x_b), with
/// outward normal n = -1 at x = 0 and n = +1 at x = L. Natural conditions
/// (Neumann, Robin) enter through that term; the essential condition
/// (Dirichlet) replaces the equation of the end node.
///
/// New kinds of condition are added by deriving from this class, without
/// changing the assembler or the solver.
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    /// Imposes the condition on `system` at the node of `end`.
    virtual void apply(LinearSystem& system, End end) const = 0;

protected:
    BoundaryCondition() = default;
    BoundaryCondition(const BoundaryCondition&) = default;
    BoundaryCondition& operator=(const BoundaryCondition&) = default;
    BoundaryCondition(BoundaryCondition&&) = default;
    BoundaryCondition& operator=(BoundaryCondition&&) = default;
};

/// Prescribed displacement: u = value.
///
/// The known column is moved to the right-hand side before the row and
/// column are cleared, so the matrix stays symmetric and tridiagonal.
class Dirichlet final : public BoundaryCondition {
public:
    explicit Dirichlet(double displacement);
    void apply(LinearSystem& system, End end) const override;

private:
    double displacement_;
};

/// Prescribed axial force: EA u' = P.
class Neumann final : public BoundaryCondition {
public:
    explicit Neumann(double axial_force);
    void apply(LinearSystem& system, End end) const override;

private:
    double axial_force_;
};

/// Mixed condition: alpha u + beta EA u' = g.
///
/// With beta = 0 it reduces to the Dirichlet condition u = g / alpha.
class Robin final : public BoundaryCondition {
public:
    /// \throws std::invalid_argument if alpha and beta are both zero.
    Robin(double alpha, double beta, double g);
    void apply(LinearSystem& system, End end) const override;

private:
    double alpha_;
    double beta_;
    double g_;
};

} // namespace fem1d
