#include "fem1d/assembler.hpp"

#include "stub.hpp"

namespace fem1d {

LinearSystem Assembler::assemble(const Mesh& /*mesh*/, const BarProperties& /*bar*/) {
    detail::not_implemented();
}

} // namespace fem1d
