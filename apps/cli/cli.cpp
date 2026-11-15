#include "cli/cli.hpp"

#include "cli/cases.hpp"
#include "fem1d/bar_solver.hpp"
#include "fem1d/io.hpp"
#include "fem1d/mesh.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <functional>
#include <map>
#include <ostream>
#include <sstream>

namespace fem1d::cli {

namespace {

const std::vector<std::string> solver_names{"thomas", "lu", "dense"};

std::string join(const std::vector<std::string>& items, const char* separator) {
    std::string result;
    for (const auto& item : items) {
        result += (result.empty() ? "" : separator) + item;
    }
    return result;
}

void check_choice(const std::string& flag, const std::string& value,
                  const std::vector<std::string>& choices) {
    if (std::find(choices.begin(), choices.end(), value) == choices.end()) {
        throw UsageError(flag + ": invalid choice '" + value + "' (choose from " +
                         join(choices, ", ") + ")");
    }
}

double parse_double(const std::string& flag, const std::string& text) {
    try {
        std::size_t used = 0;
        const double value = std::stod(text, &used);
        if (used == text.size()) {
            return value;
        }
    } catch (const std::logic_error&) {
        // Reported below.
    }
    throw UsageError(flag + ": expected a number, got '" + text + "'");
}

std::size_t parse_positive_int(const std::string& flag, const std::string& text) {
    // std::stoul accepts "-3" and wraps it around, so check the digits first.
    const bool digits_only = !text.empty() && std::all_of(text.begin(), text.end(), [](char c) {
        return std::isdigit(static_cast<unsigned char>(c)) != 0;
    });
    if (digits_only) {
        try {
            const auto value = std::stoul(text);
            if (value > 0) {
                return value;
            }
        } catch (const std::out_of_range&) {
            // Reported below.
        }
    }
    throw UsageError(flag + ": expected a positive integer, got '" + text + "'");
}

BcKind parse_kind(const std::string& flag, const std::string& text) {
    check_choice(flag, text, {"dirichlet", "neumann", "robin"});
    if (text == "dirichlet") {
        return BcKind::Dirichlet;
    }
    return text == "neumann" ? BcKind::Neumann : BcKind::Robin;
}

using Setter = std::function<void(Options&, const std::string& flag, const std::string& value)>;

/// Every flag that takes a value, with what it sets.
const std::map<std::string, Setter>& setters() {
    // Shorthand for flags that set a double.
    const auto number = [](EndOptions Options::*end, double EndOptions::*field) -> Setter {
        return [end, field](Options& o, const std::string& flag, const std::string& value) {
            (o.*end).*field = parse_double(flag, value);
        };
    };

    static const std::map<std::string, Setter> table{
        {"--case",
         [](Options& o, const std::string& flag, const std::string& value) {
             check_choice(flag, value, reference_case_names());
             o.case_name = value;
         }},
        {"--elements",
         [](Options& o, const std::string& flag, const std::string& value) {
             o.elements = parse_positive_int(flag, value);
         }},
        {"--left-type", [](Options& o, const std::string& flag,
                           const std::string& value) { o.left.kind = parse_kind(flag, value); }},
        {"--right-type", [](Options& o, const std::string& flag,
                            const std::string& value) { o.right_kind = parse_kind(flag, value); }},
        {"--u0", number(&Options::left, &EndOptions::u)},
        {"--P0", number(&Options::left, &EndOptions::P)},
        {"--alpha0", number(&Options::left, &EndOptions::alpha)},
        {"--beta0", number(&Options::left, &EndOptions::beta)},
        {"--g0", number(&Options::left, &EndOptions::g)},
        {"--uL", number(&Options::right, &EndOptions::u)},
        {"--PL", [](Options& o, const std::string& flag,
                    const std::string& value) { o.right_force = parse_double(flag, value); }},
        {"--alphaL", number(&Options::right, &EndOptions::alpha)},
        {"--betaL", number(&Options::right, &EndOptions::beta)},
        {"--gL", number(&Options::right, &EndOptions::g)},
        {"--out",
         [](Options& o, const std::string&, const std::string& value) { o.out_dir = value; }},
        {"--solver",
         [](Options& o, const std::string& flag, const std::string& value) {
             check_choice(flag, value, solver_names);
             o.solver = value;
         }},
    };
    return table;
}

std::string format_number(double value) {
    std::ostringstream s;
    s << value;
    return s.str();
}

} // namespace

Options parse_args(const std::vector<std::string>& args) {
    Options options;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        if (arg == "--help" || arg == "-h") {
            options.help = true;
            continue;
        }

        std::string flag = arg;
        std::optional<std::string> value;
        if (const auto eq = arg.find('='); eq != std::string::npos) {
            flag = arg.substr(0, eq);
            value = arg.substr(eq + 1);
        }

        const auto setter = setters().find(flag);
        if (setter == setters().end()) {
            throw UsageError("unrecognised argument '" + arg + "'");
        }
        if (!value) {
            if (i + 1 == args.size()) {
                throw UsageError(flag + ": expected a value");
            }
            value = args[++i];
        }
        setter->second(options, flag, *value);
    }
    return options;
}

ResolvedEnds resolve_ends(const Options& options, double tip_load) {
    ResolvedEnds ends{options.left, options.right};
    ends.right.kind = options.right_kind.value_or(BcKind::Neumann);
    ends.right.P = options.right_force.value_or(tip_load);
    return ends;
}

std::unique_ptr<BoundaryCondition> make_boundary_condition(const EndOptions& end) {
    switch (end.kind) {
    case BcKind::Dirichlet:
        return std::make_unique<Dirichlet>(end.u);
    case BcKind::Neumann:
        return std::make_unique<Neumann>(end.P);
    case BcKind::Robin:
        return std::make_unique<Robin>(end.alpha, end.beta, end.g);
    }
    throw std::logic_error("make_boundary_condition: unknown kind");
}

std::unique_ptr<LinearSolver> make_linear_solver(const std::string& name) {
    if (name == "thomas") {
        return std::make_unique<ThomasSolver>();
    }
    if (name == "lu") {
        return std::make_unique<TridiagonalLuSolver>();
    }
    if (name == "dense") {
        return std::make_unique<DenseLuSolver>();
    }
    throw UsageError("--solver: invalid choice '" + name + "'");
}

std::string describe(const EndOptions& end) {
    switch (end.kind) {
    case BcKind::Dirichlet:
        return "dirichlet(u=" + format_number(end.u) + ")";
    case BcKind::Neumann:
        return "neumann(P=" + format_number(end.P) + ")";
    case BcKind::Robin:
        return "robin(alpha=" + format_number(end.alpha) + ", beta=" + format_number(end.beta) +
               ", g=" + format_number(end.g) + ")";
    }
    return "unknown";
}

std::string usage() {
    return "usage: fem1d [options]\n"
           "\n"
           "Solves -(E A u')' = f A on (0, 1) with linear finite elements and writes\n"
           "<out>/fem_<case>.csv. The flags are those of the Python reference solver's\n"
           "`python -m src.training.run_fem`, plus --solver.\n"
           "\n"
           "  --case NAME            " +
           join(reference_case_names(), " | ") +
           " (default body_force)\n"
           "  --elements N           number of elements (default 40)\n"
           "  --left-type KIND       dirichlet | neumann | robin (default dirichlet)\n"
           "  --right-type KIND      dirichlet | neumann | robin (default: neumann with\n"
           "                         the case's tip load)\n"
           "  --u0, --uL VALUE       Dirichlet:  u = VALUE (default 0)\n"
           "  --P0, --PL VALUE       Neumann:    EA u' = VALUE (default 0; --PL defaults\n"
           "                         to the case's tip load)\n"
           "  --alpha0, --alphaL     Robin:      alpha u + beta EA u' = g\n"
           "  --beta0, --betaL       (defaults alpha = 0, beta = 1, g = 0)\n"
           "  --g0, --gL\n"
           "  --solver NAME          " +
           join(solver_names, " | ") +
           " (default lu)\n"
           "  --out DIR              output directory (default data/outputs)\n"
           "  -h, --help             show this message\n";
}

int run(const std::vector<std::string>& args, std::ostream& out, std::ostream& err) {
    Options options;
    try {
        options = parse_args(args);
    } catch (const UsageError& e) {
        err << "fem1d: error: " << e.what() << "\n\n" << usage();
        return 2;
    }
    if (options.help) {
        out << usage();
        return 0;
    }

    try {
        const auto problem = reference_case(options.case_name);
        const auto ends = resolve_ends(options, problem.tip_load);
        out << "[FEM] BCs: left=" << describe(ends.left) << " | right=" << describe(ends.right)
            << '\n';

        const BarSolver solver(make_linear_solver(options.solver));
        const auto solution =
            solver.solve(Mesh::uniform(1.0, options.elements), problem.bar,
                         *make_boundary_condition(ends.left), *make_boundary_condition(ends.right));

        const std::filesystem::path dir(options.out_dir);
        std::filesystem::create_directories(dir);
        const auto path = dir / ("fem_" + options.case_name + ".csv");
        write_csv(path.string(), solution);
        out << "Saved FEM solution to " << path.string() << '\n';
        return 0;
    } catch (const std::exception& e) {
        err << "fem1d: error: " << e.what() << '\n';
        return 1;
    }
}

} // namespace fem1d::cli
