#pragma once

#include "fem1d/bar_properties.hpp"

#include <array>

namespace fem1d {

/// Local stiffness matrix and load vector of a two-node element.
struct ElementMatrices {
    std::array<std::array<double, 2>, 2> stiffness{};
    std::array<double, 2> load{};
};

/// Two-node bar element with linear shape functions.
///
/// E A and f A are evaluated once, at the element midpoint (one-point
/// quadrature), exactly as in the Python reference solver:
///
///     k = (E A / h) [[1, -1], [-1, 1]],    f_e = (f A h / 2) [1, 1].
class LinearBarElement {
public:
    /// \throws std::invalid_argument unless `x_left < x_right`.
    LinearBarElement(double x_left, double x_right);

    [[nodiscard]] double length() const noexcept;
    [[nodiscard]] double midpoint() const noexcept;

    [[nodiscard]] ElementMatrices compute(const BarProperties& bar) const;

private:
    double x_left_;
    double x_right_;
};

} // namespace fem1d
