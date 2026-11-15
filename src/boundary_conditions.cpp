#include "fem1d/boundary_conditions.hpp"

#include <stdexcept>

namespace fem1d {

namespace {

struct EndNode {
    std::size_t index;
    double outward_normal;
};

EndNode end_node(const LinearSystem& system, End end) {
    if (end == End::Left) {
        return {0, -1.0};
    }
    return {system.matrix.size() - 1, 1.0};
}

} // namespace

Dirichlet::Dirichlet(double displacement) : displacement_(displacement) {}

void Dirichlet::apply(LinearSystem& system, End end) const {
    auto& k = system.matrix;
    auto& f = system.rhs;
    const std::size_t i = end_node(system, end).index;

    // Move the known column to the right-hand side, then clear the row and
    // the column. Only the neighbours of node i are coupled to it.
    if (i > 0) {
        f[i - 1] -= k(i - 1, i) * displacement_;
        k.at(i - 1, i) = 0.0;
        k.at(i, i - 1) = 0.0;
    }
    if (i + 1 < k.size()) {
        f[i + 1] -= k(i + 1, i) * displacement_;
        k.at(i + 1, i) = 0.0;
        k.at(i, i + 1) = 0.0;
    }
    k.at(i, i) = 1.0;
    f[i] = displacement_;
}

Neumann::Neumann(double axial_force) : axial_force_(axial_force) {}

void Neumann::apply(LinearSystem& system, End end) const {
    const auto [i, n] = end_node(system, end);
    system.rhs[i] += n * axial_force_;
}

Robin::Robin(double alpha, double beta, double g) : alpha_(alpha), beta_(beta), g_(g) {
    if (alpha_ == 0.0 && beta_ == 0.0) {
        throw std::invalid_argument("Robin: alpha and beta cannot both be zero");
    }
}

void Robin::apply(LinearSystem& system, End end) const {
    if (beta_ == 0.0) {
        Dirichlet(g_ / alpha_).apply(system, end);
        return;
    }
    // EA u' = (g - alpha u) / beta enters through the boundary term n EA u' v.
    const auto [i, n] = end_node(system, end);
    system.matrix.at(i, i) += n * alpha_ / beta_;
    system.rhs[i] += n * g_ / beta_;
}

} // namespace fem1d
