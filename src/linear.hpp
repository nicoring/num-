#include <optional>

#include "matrix.hpp"

#ifndef LINEAR_H
#define LINEAR_H

namespace num {

matrix<double> train_linear_model(const matrix<double>& X,
                                  const matrix<double>& y, const double lr = .1,
                                  const int iterations = 1000) {
    matrix<double> coefs(X.shape().cols, 1, 0);
    const double weighted_lr = lr / double(X.shape().rows);

    for (int i = 0; i < iterations; ++i) {
        const matrix<double> grad = (2 * X.transpose()).dot(X.dot(coefs) - y);
        coefs = coefs - (weighted_lr * grad);
    }

    return coefs;
}
} // namespace num

#endif