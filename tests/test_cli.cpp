#include "cli/cli.hpp"
#include "fem1d/io.hpp"
#include "test_helpers.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using fem1d::cli::BcKind;
using fem1d::cli::parse_args;
using fem1d::cli::run;
using fem1d::cli::UsageError;
using fem1d::test::max_abs_diff;

namespace {

/// An empty directory of its own for each test, since ctest may run them in parallel.
std::filesystem::path scratch_dir(const std::string& name) {
    const auto dir = std::filesystem::temp_directory_path() / ("fem1d_cli_test_" + name);
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    return dir;
}

struct RunResult {
    int status;
    std::string out;
    std::string err;
};

RunResult run_cli(const std::vector<std::string>& args) {
    std::ostringstream out;
    std::ostringstream err;
    const int status = run(args, out, err);
    return {status, out.str(), err.str()};
}

fem1d::Solution golden(const std::string& name) {
    return fem1d::read_csv(std::string(FEM1D_TEST_DATA_DIR) + "/" + name + ".csv");
}

} // namespace

TEST_CASE("No arguments give the defaults of the Python CLI", "[cli]") {
    const auto options = parse_args({});
    CHECK(options.case_name == "body_force");
    CHECK(options.elements == 40);
    CHECK(options.left.kind == BcKind::Dirichlet);
    CHECK(options.left.u == 0.0);
    CHECK_FALSE(options.right_kind.has_value());
    CHECK_FALSE(options.right_force.has_value());
    CHECK(options.out_dir == "data/outputs");
    CHECK(options.solver == "lu");
}

TEST_CASE("Flags take their value as the next argument or after '='", "[cli]") {
    const auto a = parse_args({"--elements", "80", "--right-type", "robin", "--alphaL", "10"});
    const auto b = parse_args({"--elements=80", "--right-type=robin", "--alphaL=10"});
    for (const auto& options : {a, b}) {
        CHECK(options.elements == 80);
        REQUIRE(options.right_kind.has_value());
        CHECK(*options.right_kind == BcKind::Robin);
        CHECK(options.right.alpha == 10.0);
    }
}

TEST_CASE("Every boundary parameter has its own flag", "[cli]") {
    const auto o =
        parse_args({"--u0", "1", "--P0", "2", "--alpha0", "3", "--beta0", "4", "--g0", "5",
                    "--uL", "6", "--PL", "7", "--alphaL", "8", "--betaL", "9", "--gL", "10"});
    CHECK(o.left.u == 1.0);
    CHECK(o.left.P == 2.0);
    CHECK(o.left.alpha == 3.0);
    CHECK(o.left.beta == 4.0);
    CHECK(o.left.g == 5.0);
    CHECK(o.right.u == 6.0);
    REQUIRE(o.right_force.has_value());
    CHECK(*o.right_force == 7.0);
    CHECK(o.right.alpha == 8.0);
    CHECK(o.right.beta == 9.0);
    CHECK(o.right.g == 10.0);
}

TEST_CASE("Malformed command lines are rejected", "[cli]") {
    using Args = std::vector<std::string>;
    CHECK_THROWS_AS(parse_args(Args{"--no-such-flag", "1"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"elements"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--elements"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--elements", "0"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--elements", "-3"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--elements", "4.5"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--u0", "abc"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--u0", "1.0x"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--case", "cantilever"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--left-type", "spring"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--solver", "cholesky"}), UsageError);
    CHECK_THROWS_AS(parse_args(Args{"--help=yes"}), UsageError);
}

TEST_CASE("The right end defaults to the case's tip load", "[cli]") {
    SECTION("no right condition given") {
        const auto ends = fem1d::cli::resolve_ends(parse_args({}), 1.0);
        CHECK(ends.right.kind == BcKind::Neumann);
        CHECK(ends.right.P == 1.0);
    }
    SECTION("an explicit force wins") {
        const auto ends = fem1d::cli::resolve_ends(parse_args({"--PL", "0.25"}), 1.0);
        CHECK(ends.right.kind == BcKind::Neumann);
        CHECK(ends.right.P == 0.25);
    }
    SECTION("another kind of condition ignores the tip load") {
        const auto ends = fem1d::cli::resolve_ends(parse_args({"--right-type", "dirichlet"}), 1.0);
        CHECK(ends.right.kind == BcKind::Dirichlet);
    }
}

SCENARIO("The command line reproduces the Python runs behind the golden data", "[cli][golden]") {
    // The same arguments as in tests/data/README.md.
    struct Run {
        std::string golden;
        std::vector<std::string> args;
    };
    const std::vector<Run> runs{
        {"hetero_neumann",
         {"--case", "hetero", "--left-type", "dirichlet", "--u0", "0.0", "--right-type", "neumann",
          "--PL", "0.0"}},
        {"hetero_robin",
         {"--case", "hetero", "--left-type", "dirichlet", "--u0", "0.0", "--right-type", "robin",
          "--alphaL", "10.0", "--betaL", "1.0", "--gL", "0.0"}},
        {"hetero_mixed_fine",
         {"--case", "hetero", "--elements", "200", "--left-type", "neumann", "--P0", "0.5",
          "--right-type", "dirichlet", "--uL", "0.2"}},
    };

    for (const auto& r : runs) {
        GIVEN("the arguments of the Python run for " + r.golden) {
            const auto dir = scratch_dir(r.golden);
            auto args = r.args;
            args.insert(args.end(), {"--out", dir.string()});

            WHEN("fem1d is run with them") {
                const auto result = run_cli(args);

                THEN("it writes fem_hetero.csv with the Python solution") {
                    REQUIRE(result.status == 0);
                    CHECK(result.err.empty());
                    CHECK_THAT(result.out, ContainsSubstring("fem_hetero.csv"));

                    const auto written = fem1d::read_csv((dir / "fem_hetero.csv").string());
                    const auto reference = golden(r.golden);
                    REQUIRE(written.x.size() == reference.x.size());
                    CHECK(max_abs_diff(written.x, reference.x) < 1e-15);
                    CHECK(max_abs_diff(written.u, reference.u) < 1e-13);
                }
            }
        }
    }
}

SCENARIO("The tip_load case defaults to a unit force at x = L", "[cli]") {
    GIVEN("only --case tip_load") {
        const auto dir = scratch_dir("tip_load");
        WHEN("fem1d is run") {
            const auto result = run_cli({"--case", "tip_load", "--out", dir.string()});
            THEN("u = P x / (EA) = x") {
                REQUIRE(result.status == 0);
                const auto written = fem1d::read_csv((dir / "fem_tip_load.csv").string());
                CHECK(written.u.back() == Approx(1.0).epsilon(1e-12));
                CHECK(max_abs_diff(written.u, written.x) < 1e-12);
            }
        }
    }
}

TEST_CASE("Every linear solver can be chosen from the command line", "[cli]") {
    for (const std::string solver : {"thomas", "lu", "dense"}) {
        const auto dir = scratch_dir("solver_" + solver);
        const auto result =
            run_cli({"--case", "hetero", "--solver", solver, "--out", dir.string()});
        CHECK(result.status == 0);
    }
}

TEST_CASE("Errors are reported on stderr with an exit status", "[cli]") {
    SECTION("--help prints the usage and succeeds") {
        const auto result = run_cli({"--help"});
        CHECK(result.status == 0);
        CHECK_THAT(result.out, ContainsSubstring("--left-type"));
    }
    SECTION("a usage error exits with 2, as argparse does") {
        const auto result = run_cli({"--elements", "zero"});
        CHECK(result.status == 2);
        CHECK_THAT(result.err, ContainsSubstring("--elements"));
    }
    SECTION("a bar free at both ends cannot be solved") {
        const auto dir = scratch_dir("free_free");
        const auto result = run_cli({"--left-type", "neumann", "--out", dir.string()});
        CHECK(result.status == 1);
        CHECK_THAT(result.err, ContainsSubstring("singular"));
        CHECK_FALSE(std::filesystem::exists(dir / "fem_body_force.csv"));
    }
    SECTION("Thomas refuses a matrix it cannot solve stably") {
        const auto dir = scratch_dir("thomas_robin");
        const auto result = run_cli({"--left-type", "robin", "--alpha0", "10", "--solver", "thomas",
                                     "--out", dir.string()});
        CHECK(result.status == 1);
        CHECK_THAT(result.err, ContainsSubstring("diagonally dominant"));
    }
}
