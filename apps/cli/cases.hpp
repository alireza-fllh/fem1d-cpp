#pragma once

#include "fem1d/bar_properties.hpp"

#include <string>
#include <vector>

namespace fem1d::cli {

/// A predefined problem of the Python reference solver
/// (`physics.get_case_config` in pinn-vs-fem-1d-bar).
struct ReferenceCase {
    BarProperties bar;
    /// Axial force at x = L used when the right end gets no explicit condition.
    double tip_load = 0.0;
};

/// Names accepted by `reference_case`, in the order the Python CLI lists them.
[[nodiscard]] const std::vector<std::string>& reference_case_names();

/// - `body_force`: E = A = 1, f = 1, tip load 0.
/// - `tip_load`:   E = A = 1, f = 0, tip load 1.
/// - `hetero`:     E = 1 for x < 0.5 and 3 otherwise, A = 1 + x/2,
///                 f = 2 sin(2 pi x), tip load 0.
/// \throws std::invalid_argument for any other name.
[[nodiscard]] ReferenceCase reference_case(const std::string& name);

} // namespace fem1d::cli
