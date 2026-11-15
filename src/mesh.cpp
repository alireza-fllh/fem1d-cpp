#include "fem1d/mesh.hpp"

#include "stub.hpp"

namespace fem1d {

Mesh::Mesh(std::vector<double> /*nodes*/) {
    detail::not_implemented();
}

Mesh Mesh::uniform(double /*length*/, std::size_t /*num_elements*/) {
    detail::not_implemented();
}

std::size_t Mesh::num_nodes() const noexcept {
    return nodes_.size();
}

std::size_t Mesh::num_elements() const noexcept {
    return nodes_.size() - 1;
}

const std::vector<double>& Mesh::nodes() const noexcept {
    return nodes_;
}

double Mesh::element_length(std::size_t /*element*/) const {
    detail::not_implemented();
}

double Mesh::element_midpoint(std::size_t /*element*/) const {
    detail::not_implemented();
}

} // namespace fem1d
