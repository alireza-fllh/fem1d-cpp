#pragma once

#include <cstddef>
#include <vector>

namespace fem1d {

/// One-dimensional mesh of two-node line elements.
///
/// Element `e` spans nodes `e` and `e + 1`.
class Mesh {
public:
    /// Builds a mesh from nodal coordinates.
    /// \throws std::invalid_argument if there are fewer than two nodes or the
    ///         coordinates are not strictly increasing.
    explicit Mesh(std::vector<double> nodes);

    /// Uniform mesh of `num_elements` elements on [0, length], with nodes
    /// placed like `numpy.linspace(0, length, num_elements + 1)`.
    /// \throws std::invalid_argument if `length <= 0` or `num_elements == 0`.
    static Mesh uniform(double length, std::size_t num_elements);

    [[nodiscard]] std::size_t num_nodes() const noexcept;
    [[nodiscard]] std::size_t num_elements() const noexcept;
    [[nodiscard]] const std::vector<double>& nodes() const noexcept;

    /// \throws std::out_of_range if `element >= num_elements()`.
    [[nodiscard]] double element_length(std::size_t element) const;

    /// \throws std::out_of_range if `element >= num_elements()`.
    [[nodiscard]] double element_midpoint(std::size_t element) const;

private:
    std::vector<double> nodes_;
};

} // namespace fem1d
