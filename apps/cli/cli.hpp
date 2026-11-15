#pragma once

// Command-line front end with the same flags as the Python reference solver's
// `python -m src.training.run_fem`, so one command line runs either code.
// It lives in a library, apart from main(), so the tests can drive it.

#include "fem1d/boundary_conditions.hpp"
#include "fem1d/linear_solver.hpp"

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem1d::cli {

/// Bad command line: unknown flag, missing value, malformed number, ...
class UsageError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

enum class BcKind { Dirichlet, Neumann, Robin };

/// Parameters of one end; which ones are used depends on `kind`.
struct EndOptions {
    BcKind kind = BcKind::Dirichlet;
    double u = 0.0;     // Dirichlet: u = u
    double P = 0.0;     // Neumann:   EA u' = P
    double alpha = 0.0; // Robin:     alpha u + beta EA u' = g
    double beta = 1.0;
    double g = 0.0;
};

struct Options {
    std::string case_name = "body_force";
    std::size_t elements = 40;
    EndOptions left;
    /// Unset means "use the case default": Neumann with the case's tip load.
    std::optional<BcKind> right_kind;
    EndOptions right;
    /// Unset means the case's tip load.
    std::optional<double> right_force;
    std::string out_dir = "data/outputs";
    std::string solver = "lu";
    bool help = false;
};

/// Parses the arguments after the program name. Accepts `--flag value` and
/// `--flag=value`.
/// \throws UsageError on any malformed argument.
[[nodiscard]] Options parse_args(const std::vector<std::string>& args);

/// The conditions at both ends after the case defaults are filled in.
struct ResolvedEnds {
    EndOptions left;
    EndOptions right;
};
[[nodiscard]] ResolvedEnds resolve_ends(const Options& options, double tip_load);

[[nodiscard]] std::unique_ptr<BoundaryCondition> make_boundary_condition(const EndOptions& end);

/// "thomas", "lu" (pivoting tridiagonal LU) or "dense".
/// \throws UsageError for any other name.
[[nodiscard]] std::unique_ptr<LinearSolver> make_linear_solver(const std::string& name);

/// Human-readable form, e.g. "robin(alpha=10, beta=1, g=0)".
[[nodiscard]] std::string describe(const EndOptions& end);

[[nodiscard]] std::string usage();

/// Runs the whole program and writes `<out>/fem_<case>.csv`.
/// \returns 0 on success, 1 if the problem cannot be solved, 2 on a usage
///          error (as argparse does).
int run(const std::vector<std::string>& args, std::ostream& out, std::ostream& err);

} // namespace fem1d::cli
