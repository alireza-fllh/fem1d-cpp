# fem1d-cpp

A C++17 finite element solver for the axial displacement of a 1D elastic bar,

$$-\frac{d}{dx}\left(E(x)A(x)\frac{du}{dx}\right) = f(x)A(x), \quad x \in (0, L),$$

with Dirichlet, Neumann and Robin conditions at either end.

It is a port of the Python FEM reference solver in
[pinn-vs-fem-1d-bar](https://github.com/alireza-fllh/pinn-vs-fem-1d-bar), and is tested against
closed-form solutions and against that solver's output.

> **Status:** work in progress. The public API and the test suite are in place; the
> implementation follows test-first.

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

## License

MIT
