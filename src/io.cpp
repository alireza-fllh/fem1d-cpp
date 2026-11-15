#include "fem1d/io.hpp"

#include <fstream>
#include <iomanip>
#include <istream>
#include <limits>
#include <locale>
#include <ostream>
#include <stdexcept>

namespace fem1d {

namespace {

void strip_carriage_return(std::string& line) {
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
}

double parse_number(const std::string& text, std::size_t line_number) {
    const auto fail = [&] {
        return std::runtime_error("read_csv: line " + std::to_string(line_number) +
                                  ": not a number: '" + text + "'");
    };
    try {
        std::size_t used = 0;
        const double value = std::stod(text, &used);
        if (used != text.size()) {
            throw fail();
        }
        return value;
    } catch (const std::logic_error&) {
        // std::stod reports bad input as invalid_argument or out_of_range.
        throw fail();
    }
}

} // namespace

void write_csv(std::ostream& out, const Solution& solution) {
    if (solution.x.size() != solution.u.size()) {
        throw std::invalid_argument("write_csv: x and u have different sizes");
    }
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "x,u\n";
    for (std::size_t i = 0; i < solution.x.size(); ++i) {
        out << solution.x[i] << ',' << solution.u[i] << '\n';
    }
}

void write_csv(const std::string& path, const Solution& solution) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("write_csv: cannot open '" + path + "' for writing");
    }
    write_csv(out, solution);
    if (!out) {
        throw std::runtime_error("write_csv: error while writing '" + path + "'");
    }
}

Solution read_csv(std::istream& in) {
    std::string line;
    if (!std::getline(in, line)) {
        throw std::runtime_error("read_csv: empty input");
    }
    strip_carriage_return(line);
    if (line != "x,u") {
        throw std::runtime_error("read_csv: expected the header 'x,u', got '" + line + "'");
    }

    Solution solution;
    std::size_t line_number = 1;
    while (std::getline(in, line)) {
        ++line_number;
        strip_carriage_return(line);
        if (line.empty()) {
            continue;
        }
        const auto comma = line.find(',');
        if (comma == std::string::npos) {
            throw std::runtime_error("read_csv: line " + std::to_string(line_number) +
                                     ": expected two comma-separated values");
        }
        solution.x.push_back(parse_number(line.substr(0, comma), line_number));
        solution.u.push_back(parse_number(line.substr(comma + 1), line_number));
    }
    return solution;
}

Solution read_csv(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("read_csv: cannot open '" + path + "'");
    }
    return read_csv(in);
}

} // namespace fem1d
