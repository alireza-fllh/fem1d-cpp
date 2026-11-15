#include "fem1d/assembler.hpp"

#include "fem1d/element.hpp"

#include <vector>

namespace fem1d {

LinearSystem Assembler::assemble(const Mesh& mesh, const BarProperties& bar) {
    const std::size_t n = mesh.num_nodes();
    LinearSystem system{TridiagonalMatrix(n), std::vector<double>(n, 0.0)};
    const auto& x = mesh.nodes();

    for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
        const auto m = LinearBarElement(x[e], x[e + 1]).compute(bar);

        // Element e connects global nodes e and e + 1.
        for (std::size_t a = 0; a < 2; ++a) {
            for (std::size_t b = 0; b < 2; ++b) {
                system.matrix.at(e + a, e + b) += m.stiffness[a][b];
            }
            system.rhs[e + a] += m.load[a];
        }
    }
    return system;
}

} // namespace fem1d
