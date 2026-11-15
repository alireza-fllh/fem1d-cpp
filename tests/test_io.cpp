#include "fem1d/io.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <sstream>
#include <stdexcept>

using fem1d::read_csv;
using fem1d::Solution;
using fem1d::write_csv;

TEST_CASE("Writing and reading a solution gives back the exact same doubles", "[io]") {
    const Solution original{{0.0, 0.1, 1.0 / 3.0, 1.0}, {0.0, -2.5e-17, std::sqrt(2.0), 1e300}};

    std::stringstream buffer;
    write_csv(buffer, original);
    const auto restored = read_csv(buffer);

    CHECK(restored.x == original.x);
    CHECK(restored.u == original.u);
}

TEST_CASE("The CSV layout matches the Python reference solver", "[io]") {
    std::stringstream buffer;
    write_csv(buffer, Solution{{0.0, 0.5}, {1.0, 2.0}});

    std::string header;
    std::getline(buffer, header);
    CHECK(header == "x,u");
}

TEST_CASE("Reading accepts numpy's output format and Windows line endings", "[io]") {
    std::istringstream in("x,u\r\n"
                          "0.000000000000000000e+00,0.000000000000000000e+00\r\n"
                          "2.500000000000000139e-02,6.284344774907086283e-03\r\n");
    const auto s = read_csv(in);

    REQUIRE(s.x.size() == 2);
    CHECK(s.x[1] == 0.025);
    CHECK(s.u[1] == 6.284344774907086283e-03);
}

TEST_CASE("Reading rejects malformed input", "[io]") {
    std::istringstream no_header("0.0,1.0\n");
    CHECK_THROWS_AS(read_csv(no_header), std::runtime_error);

    std::istringstream one_column("x,u\n0.0\n");
    CHECK_THROWS_AS(read_csv(one_column), std::runtime_error);

    std::istringstream not_a_number("x,u\n0.0,abc\n");
    CHECK_THROWS_AS(read_csv(not_a_number), std::runtime_error);

    std::istringstream trailing_text("x,u\n0.0,1.0xyz\n");
    CHECK_THROWS_AS(read_csv(trailing_text), std::runtime_error);
}

TEST_CASE("Reading a file that does not exist throws", "[io]") {
    CHECK_THROWS_AS(read_csv(std::string("does/not/exist.csv")), std::runtime_error);
}
