#include "fem1d/mesh.hpp"

#include <stdexcept>
#include <utility>

namespace fem1d {

namespace {

void check_element(const Mesh& mesh, std::size_t element) {
    if (element >= mesh.num_elements()) {
        throw std::out_of_range("Mesh: element index out of range");
    }
}

} // namespace

Mesh::Mesh(std::vector<double> nodes) : nodes_(std::move(nodes)) {
    if (nodes_.size() < 2) {
        throw std::invalid_argument("Mesh: at least two nodes are required");
    }
    for (std::size_t i = 0; i + 1 < nodes_.size(); ++i) {
        // Written as !(a < b) so that NaN coordinates are rejected too.
        if (!(nodes_[i] < nodes_[i + 1])) {
            throw std::invalid_argument("Mesh: nodes must be strictly increasing");
        }
    }
}

Mesh Mesh::uniform(double length, std::size_t num_elements) {
    if (!(length > 0.0)) {
        throw std::invalid_argument("Mesh::uniform: length must be positive");
    }
    if (num_elements == 0) {
        throw std::invalid_argument("Mesh::uniform: at least one element is required");
    }

    // Same construction as numpy.linspace: x_i = i * h, with the last node set to L exactly.
    const double h = length / static_cast<double>(num_elements);
    std::vector<double> nodes(num_elements + 1);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        nodes[i] = static_cast<double>(i) * h;
    }
    nodes.back() = length;
    return Mesh(std::move(nodes));
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

double Mesh::element_length(std::size_t element) const {
    check_element(*this, element);
    return nodes_[element + 1] - nodes_[element];
}

double Mesh::element_midpoint(std::size_t element) const {
    check_element(*this, element);
    return 0.5 * (nodes_[element] + nodes_[element + 1]);
}

} // namespace fem1d
