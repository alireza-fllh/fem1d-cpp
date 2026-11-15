#include "fem1d/element.hpp"

#include "stub.hpp"

namespace fem1d {

LinearBarElement::LinearBarElement(double x_left, double x_right)
    : x_left_(x_left), x_right_(x_right) {
    detail::not_implemented();
}

double LinearBarElement::length() const noexcept {
    return x_right_ - x_left_;
}

double LinearBarElement::midpoint() const noexcept {
    return 0.5 * (x_left_ + x_right_);
}

ElementMatrices LinearBarElement::compute(const BarProperties& /*bar*/) const {
    detail::not_implemented();
}

} // namespace fem1d
