#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace fem1d::test {

/// Largest |a_i - b_i|.
inline double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("max_abs_diff: size mismatch");
    }
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result = std::max(result, std::abs(a[i] - b[i]));
    }
    return result;
}

/// Largest |u_i - exact(x_i)| over the nodes.
inline double max_nodal_error(const std::vector<double>& x, const std::vector<double>& u,
                              const std::function<double(double)>& exact) {
    std::vector<double> reference(x.size());
    std::transform(x.begin(), x.end(), reference.begin(), exact);
    return max_abs_diff(u, reference);
}

} // namespace fem1d::test
