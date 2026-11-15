#pragma once

#include <functional>

namespace fem1d {

/// A coefficient that varies along the bar, evaluated at a position x.
using ScalarField = std::function<double(double)>;

/// A field with the same value everywhere.
inline ScalarField constant(double value) {
    return [value](double /*x*/) { return value; };
}

/// Material, geometry and loading of the bar in -(E A u')' = f A.
///
/// `body_force` is a force per unit volume, so the load per unit length is
/// f(x) A(x). This matches the Python reference solver (fem.py), which
/// assembles f A into the load vector.
struct BarProperties {
    ScalarField youngs_modulus = constant(1.0);
    ScalarField area = constant(1.0);
    ScalarField body_force = constant(0.0);
};

} // namespace fem1d
