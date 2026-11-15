#pragma once

#include <stdexcept>

namespace fem1d::detail {

[[noreturn]] inline void not_implemented() {
    throw std::logic_error("not implemented");
}

} // namespace fem1d::detail
