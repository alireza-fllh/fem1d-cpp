// Wall-clock comparison of the three linear solvers on the bar's stiffness
// system, plus the cost of assembly.
//
// Each time is the fastest of several batches of repeated calls, which
// filters out most interference from the rest of the machine. Run a Release
// build: `fem1d_bench` (full) or `fem1d_bench --quick` (smoke test).

#include "cli/cases.hpp"
#include "fem1d/assembler.hpp"
#include "fem1d/boundary_conditions.hpp"
#include "fem1d/linear_solver.hpp"
#include "fem1d/mesh.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <optional>
#include <vector>

#ifndef FEM1D_BUILD_TYPE
#define FEM1D_BUILD_TYPE "unknown"
#endif

#ifdef __VERSION__
#define FEM1D_COMPILER __VERSION__
#else
#define FEM1D_COMPILER "unknown"
#endif

namespace {

using Clock = std::chrono::steady_clock;

/// Keeps the compiler from discarding results it can see are unused.
volatile double sink = 0.0;

/// Seconds per call of `f`: the fastest of `batches` batches, each of enough
/// calls to last at least `min_batch_seconds`.
template <class F>
double seconds_per_call(F&& f, double min_batch_seconds, int batches = 5) {
    std::size_t calls = 1;
    double best = std::numeric_limits<double>::infinity();
    for (int b = 0; b < batches;) {
        const auto start = Clock::now();
        for (std::size_t i = 0; i < calls; ++i) {
            f();
        }
        const double elapsed = std::chrono::duration<double>(Clock::now() - start).count();
        if (elapsed < min_batch_seconds) {
            calls *= 2; // Too short to time reliably; calibrate and retry.
            continue;
        }
        best = std::min(best, elapsed / static_cast<double>(calls));
        ++b;
    }
    return best;
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result = std::max(result, std::abs(a[i] - b[i]));
    }
    return result;
}

struct Row {
    std::size_t unknowns;
    double assemble;
    double thomas;
    double tridiagonal_lu;
    std::optional<double> dense_lu;
    double lu_vs_thomas;
    std::optional<double> dense_vs_lu;
};

Row measure(std::size_t elements, std::size_t max_dense_unknowns, double min_batch_seconds) {
    const auto mesh = fem1d::Mesh::uniform(1.0, elements);
    const auto bar = fem1d::cli::reference_case("hetero").bar;

    Row row{};
    row.unknowns = mesh.num_nodes();
    row.assemble = seconds_per_call(
        [&] { sink = sink + fem1d::Assembler::assemble(mesh, bar).rhs.back(); }, min_batch_seconds);

    // Clamped at x = 0, free at x = L: diagonally dominant, so Thomas applies.
    auto system = fem1d::Assembler::assemble(mesh, bar);
    fem1d::Dirichlet(0.0).apply(system, fem1d::End::Left);
    fem1d::Neumann(0.0).apply(system, fem1d::End::Right);

    const auto time_solver = [&](const fem1d::LinearSolver& solver) {
        return seconds_per_call(
            [&] { sink = sink + solver.solve(system.matrix, system.rhs).back(); },
            min_batch_seconds);
    };

    const fem1d::ThomasSolver thomas;
    const fem1d::TridiagonalLuSolver tridiagonal_lu;
    row.thomas = time_solver(thomas);
    row.tridiagonal_lu = time_solver(tridiagonal_lu);
    const auto u_lu = tridiagonal_lu.solve(system.matrix, system.rhs);
    row.lu_vs_thomas = max_abs_diff(u_lu, thomas.solve(system.matrix, system.rhs));

    if (row.unknowns <= max_dense_unknowns) {
        const fem1d::DenseLuSolver dense_lu;
        row.dense_lu = time_solver(dense_lu);
        row.dense_vs_lu = max_abs_diff(dense_lu.solve(system.matrix, system.rhs), u_lu);
    }
    return row;
}

void print_time(double seconds) {
    if (seconds < 1e-3) {
        std::printf(" %9.2f us |", seconds * 1e6);
    } else {
        std::printf(" %9.2f ms |", seconds * 1e3);
    }
}

} // namespace

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
    if (argc > 2 || (argc == 2 && !quick)) {
        std::fprintf(stderr, "usage: fem1d_bench [--quick]\n");
        return 2;
    }

    const std::vector<std::size_t> elements =
        quick ? std::vector<std::size_t>{99, 399}
              : std::vector<std::size_t>{99, 399, 1599, 9'999, 99'999, 999'999};
    const std::size_t max_dense_unknowns = 1600;
    const double min_batch_seconds = quick ? 1e-3 : 0.05;

    std::printf("Build type: %s, compiler: %s\n\n", FEM1D_BUILD_TYPE, FEM1D_COMPILER);
    std::printf("| unknowns | assemble     | Thomas       | tridiag. LU  | dense LU     |"
                " Thomas ns/unknown | LU ns/unknown | max|LU-Thomas| | max|dense-LU| |\n");
    std::printf("|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");

    for (const auto n : elements) {
        const auto row = measure(n, max_dense_unknowns, min_batch_seconds);
        const auto per_unknown = [&](double seconds) {
            return seconds * 1e9 / static_cast<double>(row.unknowns);
        };
        std::printf("| %8zu |", row.unknowns);
        print_time(row.assemble);
        print_time(row.thomas);
        print_time(row.tridiagonal_lu);
        if (row.dense_lu) {
            print_time(*row.dense_lu);
        } else {
            std::printf("            - |");
        }
        std::printf(" %17.2f | %13.2f | %14.1e |", per_unknown(row.thomas),
                    per_unknown(row.tridiagonal_lu), row.lu_vs_thomas);
        if (row.dense_vs_lu) {
            std::printf(" %13.1e |\n", *row.dense_vs_lu);
        } else {
            std::printf("             - |\n");
        }
    }
    return 0;
}
