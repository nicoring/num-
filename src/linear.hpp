#pragma once

#include "linalg.hpp"
#include "matrix.hpp"

namespace num {

matrix<double> gradient_linear_model(const matrix<double>& X, const matrix<double>& y,
                                     const double lr = .1, const double l2_norm = 0.,
                                     const int iterations = 1000) {
    matrix<double> coefs(X.shape().cols, 1, 0);
    const double N = double(X.shape().rows);

    for (int i = 0; i < iterations; ++i) {
        matrix<double> grad = (2 * X.transpose()).dot(X.dot(coefs) - y) / N;
        if (l2_norm > 0) {
            grad = grad + 2 * l2_norm * coefs;
        }
        coefs = coefs - (lr * grad);
    }

    return coefs;
}

matrix<double> exact_linear_model(const matrix<double>& X, const matrix<double> y) {
    const matrix<double> Xt = X.transpose();
    matrix<double> coefs = inverse(Xt.dot(X)).dot(Xt).dot(y);
    return coefs;
}

matrix<double> exact_linear_model(const matrix<double>& X, const matrix<double> y,
                                  const double l2_norm) {
    const matrix<double> Xt = X.transpose();
    const matrix<double> I = eye<double>(X.shape().rows);
    matrix<double> coefs = inverse(Xt.dot(X) + I * l2_norm).dot(Xt).dot(y);
    return coefs;
}

} // namespace num
