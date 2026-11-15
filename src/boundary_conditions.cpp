#include "fem1d/boundary_conditions.hpp"

#include "stub.hpp"

namespace fem1d {

Dirichlet::Dirichlet(double displacement) : displacement_(displacement) {}

void Dirichlet::apply(LinearSystem& /*system*/, End /*end*/) const {
    detail::not_implemented();
}

Neumann::Neumann(double axial_force) : axial_force_(axial_force) {}

void Neumann::apply(LinearSystem& /*system*/, End /*end*/) const {
    detail::not_implemented();
}

Robin::Robin(double alpha, double beta, double g) : alpha_(alpha), beta_(beta), g_(g) {
    detail::not_implemented();
}

void Robin::apply(LinearSystem& /*system*/, End /*end*/) const {
    detail::not_implemented();
}

} // namespace fem1d
