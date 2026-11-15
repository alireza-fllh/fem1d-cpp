#include "fem1d/element.hpp"

#include <stdexcept>

namespace fem1d {

LinearBarElement::LinearBarElement(double x_left, double x_right)
    : x_left_(x_left), x_right_(x_right) {
    if (!(x_left_ < x_right_)) {
        throw std::invalid_argument("LinearBarElement: x_left must be less than x_right");
    }
}

double LinearBarElement::length() const noexcept {
    return x_right_ - x_left_;
}

double LinearBarElement::midpoint() const noexcept {
    return 0.5 * (x_left_ + x_right_);
}

ElementMatrices LinearBarElement::compute(const BarProperties& bar) const {
    const double h = length();
    const double xm = midpoint();
    const double area = bar.area(xm);

    // Same order of operations as fem.py, so the results agree to round-off.
    const double k = bar.youngs_modulus(xm) * area / h;
    const double f = bar.body_force(xm) * area * h * 0.5;

    ElementMatrices m;
    m.stiffness = {{{k, -k}, {-k, k}}};
    m.load = {f, f};
    return m;
}

} // namespace fem1d
