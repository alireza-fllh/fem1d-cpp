#pragma once

#include "fem1d/bar_properties.hpp"
#include "fem1d/linear_system.hpp"
#include "fem1d/mesh.hpp"

namespace fem1d {

/// Builds the global system by summing element contributions into the
/// global degrees of freedom. Boundary conditions are not applied here.
class Assembler {
public:
    [[nodiscard]] static LinearSystem assemble(const Mesh& mesh, const BarProperties& bar);
};

} // namespace fem1d
