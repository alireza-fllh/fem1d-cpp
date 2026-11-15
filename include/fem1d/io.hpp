#pragma once

#include "fem1d/bar_solver.hpp"

#include <iosfwd>
#include <string>

namespace fem1d {

/// Writes a solution as CSV with the header `x,u`, the same layout as the
/// Python reference solver. Values are written with enough digits to read
/// back the exact same doubles.
void write_csv(std::ostream& out, const Solution& solution);

/// \throws std::runtime_error if the file cannot be written.
void write_csv(const std::string& path, const Solution& solution);

/// Reads a CSV with the header `x,u` and one node per line.
/// \throws std::runtime_error on a missing header or a malformed line.
[[nodiscard]] Solution read_csv(std::istream& in);

/// \throws std::runtime_error if the file cannot be opened or is malformed.
[[nodiscard]] Solution read_csv(const std::string& path);

} // namespace fem1d
