#pragma once

#include "linalg.hpp"
#include "matrix.hpp"

namespace num {

matrix<double> gradient_linear_model(const matrix<double>& X, const matrix<double>& y,
                                     const double lr = .1,
                                     const int iterations = 1000) {
    matrix<double> coefs(X.shape().cols, 1, 0);
    const double weighted_lr = lr / double(X.shape().rows);

    for (int i = 0; i < iterations; ++i) {
        const matrix<double> grad = (2 * X.transpose()).dot(X.dot(coefs) - y);
        coefs = coefs - (weighted_lr * grad);
    }

    return coefs;
}

matrix<double> exact_linear_model(const matrix<double>& X, const matrix<double> y) {
    const matrix<double> Xt = X.transpose();
    matrix<double> coefs = inverse(Xt.dot(X)).dot(Xt).dot(y);
    return coefs;
}

} // namespace num
