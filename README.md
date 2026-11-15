# fem1d-cpp

[![CI](https://github.com/alireza-fllh/fem1d-cpp/actions/workflows/ci.yml/badge.svg)](https://github.com/alireza-fllh/fem1d-cpp/actions/workflows/ci.yml)

A C++17 finite element solver for the axial displacement of a 1D elastic bar,

$$-\frac{d}{dx}\left(E(x)A(x)\frac{du}{dx}\right) = f(x)A(x), \quad x \in (0, L),$$

with Dirichlet, Neumann and Robin conditions at either end.

It is a port of the Python FEM reference solver in
[pinn-vs-fem-1d-bar](https://github.com/alireza-fllh/pinn-vs-fem-1d-bar) and reproduces its
output to round-off.

- **Design:** each FEM step is its own type, and boundary conditions and linear solvers are
  replaceable strategies. [DESIGN.md](DESIGN.md) explains the decisions and relates them to SOLID
  and GRASP.
- **Linear algebra:** tridiagonal storage with O(n) solvers. Thomas is guarded by a
  diagonal-dominance check, and the default pivoting tridiagonal LU follows LAPACK's `dgtsv`. A
  dense LU is kept as a reference. The tests found a case where an unguarded Thomas solver gives
  a wrong answer without warning; [DESIGN.md](DESIGN.md#the-linear-solver-is-a-strategy-and-the-default-pivots)
  tells that story.
- **Tests:** Catch2 unit tests, BDD scenarios against closed-form solutions, a measured O(h²)
  convergence rate, golden data from the Python solver, and end-to-end tests of the command line.
- **Tooling:** CMake, warnings as errors, AddressSanitizer and UndefinedBehaviorSanitizer, and
  clang-format, all checked in GitHub Actions with GCC and Clang.

## Using the library

```cpp
#include "fem1d/fem1d.hpp"

#include <iostream>

int main() {
    using namespace fem1d;

    BarProperties bar;
    bar.youngs_modulus = [](double x) { return x < 0.5 ? 1.0 : 3.0; };
    bar.area = [](double x) { return 1.0 + 0.5 * x; };
    bar.body_force = constant(1.0);

    const auto mesh = Mesh::uniform(1.0, 100);
    const BarSolver solver; // or BarSolver(std::make_unique<ThomasSolver>())
    const auto solution = solver.solve(mesh, bar, Dirichlet(0.0), Robin(10.0, 1.0, 0.0));

    write_csv(std::cout, solution); // x,u per node
}
```

Link against the `fem1d::fem1d` CMake target. Robin conditions are written as
αu + βEAu' = g, as in the Python solver.

## Repository layout

| Path | Contents |
|---|---|
| `include/fem1d/`, `src/` | The library |
| `apps/` | The `fem1d` command line |
| `bench/` | The solver benchmark |
| `tests/` | Catch2 tests; `tests/data/` holds the golden data |
| `scripts/generate_golden.py` | Regenerates the golden data from the Python solver |
| `cmake/` | Warning and sanitizer flags |

## Build and test

Requires CMake 3.25+ and a C++17 compiler. Catch2 is fetched automatically.

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Options:

| Option | Default | Effect |
|---|---|---|
| `FEM1D_BUILD_TESTS` | `ON` | Build the Catch2 test suite |
| `FEM1D_ENABLE_SANITIZERS` | `OFF` | Build with AddressSanitizer and UndefinedBehaviorSanitizer |
| `FEM1D_WARNINGS_AS_ERRORS` | `ON` | Treat compiler warnings as errors |
| `FEM1D_BUILD_BENCHMARKS` | `OFF` | Build `fem1d_bench` |

## Command line

`build/apps/fem1d` takes the same flags as the Python solver's `python -m src.training.run_fem`
and writes the same `<out>/fem_<case>.csv`, so one command line runs either code:

```bash
build/apps/fem1d --case hetero --right-type robin --alphaL 10 --betaL 1 --gL 0 --out results
```

It adds `--solver thomas|lu|dense` to pick the linear solver. `fem1d --help` lists every flag.

## Cross-check against the Python solver

`tests/data/` holds output of the Python solver's own `run_fem` command for the heterogeneous
bar (two-layer E, tapered A, sine load) with three sets of boundary conditions.
`tests/test_golden.cpp` solves the same problems and requires agreement to 1e-13; the measured
difference is below 2e-14. [`tests/data/README.md`](tests/data/README.md) records the commands
and the reference commit. To regenerate the data, in the Python repo's environment:

```bash
python scripts/generate_golden.py --reference /path/to/pinn-vs-fem-1d-bar
```

## Benchmark

`fem1d_bench` times assembly and the three linear solvers on the heterogeneous bar, clamped at
x = 0 and free at x = L. Each time is the fastest of five batches of repeated calls. One run on an
Intel i7-12700H under WSL2, GCC 13.3, Release build, single thread:

| unknowns | assemble | Thomas | tridiagonal LU | dense LU |
|---:|---:|---:|---:|---:|
| 100 | 2.1 us | 1.2 us | 1.6 us | 64 us |
| 400 | 7.6 us | 4.8 us | 6.2 us | 4.4 ms |
| 1 600 | 27 us | 19 us | 25 us | 453 ms |
| 10 000 | 171 us | 117 us | 160 us | - |
| 100 000 | 1.8 ms | 1.2 ms | 1.6 ms | - |
| 1 000 000 | 20 ms | 13 ms | 19 ms | - |

- Both tridiagonal solvers cost a constant 12-19 ns per unknown, so they are O(n). The times
  include allocating the solution and the input checks (pivot tolerance, and diagonal dominance
  for Thomas).
- Pivoting makes tridiagonal LU about 1.3 times slower than Thomas. It is still the default,
  because Thomas cannot solve every system a Robin condition can produce.
- Dense LU, the method the Python solver uses through `numpy.linalg.solve`, grows as n^3. At
  1 600 unknowns it is about 18 000 times slower than tridiagonal LU.
- All three agree: tridiagonal and dense LU give the same answer bit for bit on this
  diagonally dominant matrix, and Thomas differs by round-off (4e-14 at 10^6 unknowns).

```bash
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release -DFEM1D_BUILD_BENCHMARKS=ON
cmake --build build-release
build-release/bench/fem1d_bench
```

CI runs `fem1d_bench --quick` only to check that it works; timings from shared runners mean
little.

## License

MIT
