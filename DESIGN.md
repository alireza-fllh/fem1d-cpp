# Design

This document explains how fem1d is organised and why. It covers the main decisions, how they
relate to SOLID and GRASP, how the code is tested, and what it deliberately leaves out.

## Scope

fem1d solves the axial displacement of a 1D elastic bar,

$$-\frac{d}{dx}\left(E(x)A(x)\frac{du}{dx}\right) = f(x)A(x), \quad x \in (0, L),$$

with linear two-node elements and a Dirichlet, Neumann or Robin condition at each end. It is a
port of the FEM solver in [pinn-vs-fem-1d-bar](https://github.com/alireza-fllh/pinn-vs-fem-1d-bar)
and reproduces its output to round-off.

The problem is small on purpose. The point is the engineering around it: a clean separation of
the FEM steps, linear algebra chosen for the structure of the matrix, tests at every level, and
a build that checks itself. It is single-threaded and makes no performance claims beyond the
measured [benchmark](README.md#benchmark).

## Structure

One solve is a pipeline. Each stage is a separate type with one job:

```text
Mesh ──► LinearBarElement ──► Assembler ──► LinearSystem ──► BoundaryCondition ──► LinearSolver ──► Solution
nodes     local k, f per       sums into      tridiagonal K,    modifies K, F at     solves K u = F
          element              global K, F    vector F          one end node
```

`BarSolver` runs the pipeline. It owns a `LinearSolver` and borrows everything else.

| Type | Header | Responsibility |
|---|---|---|
| `Mesh` | `mesh.hpp` | Node coordinates; validates that they are strictly increasing (NaN is rejected) |
| `BarProperties` | `bar_properties.hpp` | E(x), A(x), f(x) as `std::function<double(double)>` |
| `LinearBarElement` | `element.hpp` | Local stiffness and load of one element (midpoint quadrature) |
| `Assembler` | `assembler.hpp` | Sums element contributions into the global system |
| `TridiagonalMatrix`, `LinearSystem` | `linear_system.hpp` | Banded storage of K, plus the right-hand side F |
| `BoundaryCondition` and subclasses | `boundary_conditions.hpp` | Imposes one condition at one end |
| `LinearSolver` and subclasses | `linear_solver.hpp` | Solves the tridiagonal system |
| `BarSolver`, `Solution` | `bar_solver.hpp` | Runs the pipeline; holds nodes and displacements |
| `write_csv`, `read_csv` | `io.hpp` | CSV in the Python solver's format |

Outside the library:

- `apps/` holds the `fem1d` command line. Its logic is a small library (`fem1d_cli`) separate from
  `main()`, so the tests can drive it.
- `bench/` holds the solver benchmark.

## Decisions

### Store the matrix as three diagonals

Linear elements in 1D couple each node only to its neighbours, so K is tridiagonal. Storing it as
three vectors makes assembly and solving O(n), where a dense matrix costs O(n²) memory and O(n³)
time. The Python solver uses a dense `numpy` array and `numpy.linalg.solve`; the benchmark shows
what that costs at 1 600 unknowns (453 ms dense against 25 µs tridiagonal).

`TridiagonalMatrix` gives two kinds of access:

- `operator()(r, c) const` reads any entry and returns 0 outside the band, so code can treat the
  matrix as a mathematical object.
- `at(r, c)` gives a reference and throws outside the band, so a write that would break the
  structure fails loudly instead of being lost.

### Boundary conditions are objects

`BoundaryCondition` has a single virtual method, `apply(LinearSystem&, End)`. `Dirichlet`,
`Neumann` and `Robin` derive from it and are `final`.

- **Adding a kind of condition** (a point spring, a contact condition) means adding a class. The
  assembler and the solver do not change.
- **Each condition knows its own sign convention.** Integration by parts leaves the boundary
  term n·EA·u'·v, with outward normal n = -1 at x = 0 and +1 at x = L. Natural conditions enter
  through that term: Neumann adds n·P to F, and Robin adds n·α/β to K and n·g/β to F.
- **Dirichlet uses elimination.** It moves the known column to the right-hand side and then
  clears the row and the column, so K stays symmetric and tridiagonal. A penalty method would
  avoid touching F but adds a large number to K and spoils its conditioning. Robin with β = 0 is
  the Dirichlet condition u = g/α and is handled as one.
- **The order of the two conditions does not matter.** Natural conditions only add to the
  diagonal and F at their own node. Dirichlet reads and clears only the off-diagonal entries next
  to its node, which the other end never changes. `test_boundary_conditions.cpp` checks every
  pair of conditions, including on a one-element mesh where both ends share those entries. K
  comes out identical and F agrees to round-off, because floating-point additions done in a
  different order can differ in the last bit.

`End` is an enum and not a second parameter object, because a 1D bar has exactly two ends.

### The linear solver is a strategy, and the default pivots

`LinearSolver` has three implementations:

| Solver | Method | Cost | Use |
|---|---|---|---|
| `ThomasSolver` | Gaussian elimination without pivoting | O(n) | Fastest; diagonally dominant matrices only |
| `TridiagonalLuSolver` | Partial pivoting, the algorithm of LAPACK's `dgtsv` | O(n) | Default; any non-singular tridiagonal matrix |
| `DenseLuSolver` | LU with partial pivoting on a dense copy | O(n³) | Reference and benchmark baseline |

`BarSolver` receives its solver through the constructor as a `std::unique_ptr<LinearSolver>` and
owns it. The default is `TridiagonalLuSolver`.

The default is not Thomas, the textbook choice for tridiagonal systems, because of a bug that
the tests found:

1. The Python solver writes the Robin condition as αu + βEAu' = g in terms of u', not of the
   outward flux. At x = L, α > 0 is a spring. At x = 0 the same α is a spring of *negative*
   stiffness: it subtracts α/β from K(0, 0).
2. With E = A = 1, 40 elements and 2u(0) + EAu'(0) = 1, K(0, 0) drops to 38 next to
   K(0, 1) = -40. The matrix is regular but neither diagonally dominant nor positive definite.
3. The first version of this port solved it with Thomas, the obvious choice. Without
   pivoting the pivots shrank towards zero around row 20, and the answer was wrong by 2e-3 with
   no error raised. The Python solver gets it right only because `numpy.linalg.solve` pivots.

The fix has two parts:

- `ThomasSolver` checks weak diagonal dominance first and throws `std::domain_error` if the
  matrix does not have it. It never returns an inaccurate answer.
- `BarSolver` defaults to `TridiagonalLuSolver`, which is still O(n) and about 1.3 times slower
  than Thomas in the benchmark.

A regression test keeps this case in the suite.

Every solver treats a pivot below n·ε·max|aᵢⱼ| as zero and throws `SingularMatrixError`, for
example for a bar that is free at both ends. The check is written as `!(|p| > tol)` so that a
NaN pivot is rejected too.

### Material fields are `std::function`

E, A and f are `std::function<double(double)>`, so a caller can pass a lambda, a function or a
functor without templates leaking into the public API. The price is one indirect call per
evaluation, three per element. The benchmark puts assembly at about 20 ns per element, in the
same range as the solve, which is acceptable here. A code that assembled millions of elements
repeatedly would template the element on the field types instead.

### Match the Python solver exactly

`LinearBarElement` evaluates EA and fA once, at the element midpoint, as the Python solver does.
Two-point Gauss quadrature would be more accurate for varying coefficients, but one-point
quadrature keeps the two codes doing the same arithmetic, so the golden tests can demand
agreement to 1e-13 and not just to discretisation error. It also keeps the O(h²) convergence of
linear elements, which `test_convergence.cpp` measures.

### Errors

Errors are exceptions from the standard hierarchy, thrown where the problem is detected:

| Exception | When |
|---|---|
| `std::invalid_argument` | Bad input: non-increasing nodes, zero elements, α = β = 0, mismatched sizes, a null solver |
| `std::out_of_range` | An index outside the matrix or the mesh |
| `std::domain_error` | `ThomasSolver` given a matrix that is not diagonally dominant |
| `fem1d::SingularMatrixError` (a `std::runtime_error`) | A numerically zero pivot |
| `std::runtime_error` | Unreadable or malformed CSV (the message gives the line number) |

The command line catches them at the top. It exits with 2 for a usage error, as argparse does,
and with 1 for a problem that cannot be solved, and it writes no output file in either case.

### Ownership and value types

- `Mesh`, `TridiagonalMatrix`, `LinearSystem` and `Solution` are plain value types, copyable and
  movable.
- The only owning pointer is `BarSolver`'s `std::unique_ptr<LinearSolver>`. There is no
  `shared_ptr` and no raw `new`.
- The polymorphic bases have protected copy and move operations, so a derived object cannot be
  sliced through a base reference.

## SOLID

| Principle | Where it shows |
|---|---|
| Single responsibility | Each stage of the pipeline is one type with one reason to change: the element formula, the assembly loop, a boundary condition, a linear solver, the file format. |
| Open/closed | New boundary conditions and new linear solvers are new classes. `Assembler` and `BarSolver` do not change. |
| Liskov substitution | `test_linear_solver.cpp` runs the same `TEMPLATE_TEST_CASE`s on all three solvers, and a scenario solves the same bar with each. The one deliberate difference, `ThomasSolver` refusing non-dominant matrices, is part of the documented interface (`\throws`), not a surprise. |
| Interface segregation | The interfaces each have one method: `apply` and `solve`. |
| Dependency inversion | `BarSolver` depends on the `LinearSolver` abstraction and receives the concrete solver from outside. The command line chooses it from `--solver`. |

## GRASP

| Pattern | Where it shows |
|---|---|
| Information expert | `LinearBarElement` computes its own matrices because it knows its nodes. Each boundary condition knows its sign convention and how it changes K and F. |
| Creator | `Assembler` creates the `LinearSystem` it fills. The CLI's `make_boundary_condition` and `make_linear_solver` create objects from parsed options, keeping that mapping out of the library. |
| Controller | `BarSolver` coordinates a solve; the command line's `run()` coordinates a program run. Neither does numerical work itself. |
| Low coupling | Boundary conditions see only `LinearSystem`, not the mesh or the material. Solvers see only the matrix and the right-hand side. |
| High cohesion | `linear_solver.cpp` holds only linear algebra. Parsing, files and printing are in `io.cpp` and `apps/`. |
| Polymorphism | Variation in boundary conditions and solvers is handled by virtual dispatch, not by `switch` statements spread through the code. The CLI's one `switch` on `BcKind` turns parsed text into those objects. |
| Protected variations | The rest of the code depends on two stable interfaces, so changes to how a condition or a solver works do not ripple outward. |

Not every class uses these patterns. `Assembler::assemble` is a static function, because
assembly has no state and only one variant. An interface would add indirection with nothing to
vary.

## Testing

The tests are written with Catch2. The first commit in the history contains the API and failing
tests only; the implementation followed. They work at several levels:

| Level | Files | What it shows |
|---|---|---|
| Unit | `test_mesh`, `test_element`, `test_linear_system`, `test_assembler`, `test_boundary_conditions`, `test_linear_solver`, `test_io` | Each type against hand-computed values and its error cases |
| Behaviour (BDD) | `test_bar_scenarios` | `SCENARIO`/`GIVEN`/`WHEN`/`THEN` against closed-form solutions: every analytic case of the Python repo, the free-free bar, and swapping the solver |
| Convergence | `test_convergence` | The nodal error falls as h² (measured rate 2 ± 0.1) for a problem whose FEM solution is not exact |
| Cross-language | `test_golden` | Agreement with the Python solver's output to 1e-13 on a heterogeneous bar |
| End-to-end | `test_cli` | The `fem1d` command line, given the Python command lines, writes the same files |

Some details that matter:

- **Golden data come from a case that is hard to get right.** Several analytic cases (a tip load
  on a uniform bar, for instance) are reproduced exactly by linear elements, so matching them
  says little. The `hetero` case has a discontinuous E, a tapered A and a sine load, so its FEM
  solution has a real discretisation error. Agreement to 1e-13 means both codes build and solve
  the same system. As a negative check, done by hand, ignoring the taper (A = 1) makes the
  comparison fail by 1e-2 or more. `scripts/generate_golden.py` regenerates the data from the Python repo's own
  command line and records the reference commit.
- **Bugs become tests.** The Thomas case above is in the suite as a regression test. A separate
  test covers the row-exchange branch of `TridiagonalLuSolver` directly, because diagonally
  dominant matrices never enter it.
- **CI** builds with GCC and Clang in Release, and with GCC in Debug under AddressSanitizer and
  UndefinedBehaviorSanitizer, with `-Wall -Wextra -Wpedantic -Wconversion -Wshadow` and other
  warnings as errors. It runs every test and the benchmark's quick mode, and checks formatting
  with clang-format.

## Not included

These were left out on purpose:

- **Higher-order elements and Gauss quadrature.** They would come in as another element type and
  a quadrature rule passed to it. The Python solver has neither, so they would break the golden
  comparison.
- **More than one dimension.** In 2D or 3D, K is no longer banded, so `TridiagonalMatrix` would
  give way to a sparse format (CSR) and the solvers to sparse direct or iterative methods. The
  `LinearSolver` interface would have to change to take a general matrix.
- **Parallelism.** At these sizes the whole solve takes milliseconds, and the tridiagonal solve
  is sequential by nature. Assembly could be parallelised by colouring the elements, but there is
  nothing here to measure the gain against, so the code makes no such claim.
- **Python bindings.** pybind11 bindings would let the Python repo call this solver directly.
  The CSV interface already makes the two interchangeable.
