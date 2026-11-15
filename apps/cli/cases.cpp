#include "cli/cases.hpp"

#include <cmath>
#include <stdexcept>

namespace fem1d::cli {

const std::vector<std::string>& reference_case_names() {
    static const std::vector<std::string> names{"body_force", "tip_load", "hetero"};
    return names;
}

ReferenceCase reference_case(const std::string& name) {
    ReferenceCase result;
    if (name == "body_force") {
        result.bar.body_force = constant(1.0);
    } else if (name == "tip_load") {
        result.tip_load = 1.0;
    } else if (name == "hetero") {
        const double pi = std::acos(-1.0);
        result.bar.youngs_modulus = [](double x) { return x < 0.5 ? 1.0 : 3.0; };
        result.bar.area = [](double x) { return 1.0 + 0.5 * x; };
        result.bar.body_force = [pi](double x) { return 2.0 * std::sin(2.0 * pi * x); };
    } else {
        throw std::invalid_argument("unknown case '" + name + "'");
    }
    return result;
}

} // namespace fem1d::cli
