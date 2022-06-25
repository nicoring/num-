#include <gtest/gtest.h>

#include "../src/linear.hpp"

TEST(LinearSGD, Fit) {
    num::matrix<double> X = num::randn(10000, 10);
    num::matrix<double> w = num::randn(10, 1);
    num::matrix<double> y = X.dot(w);
    num::matrix<double> w_est = num::gradient_linear_model(X, y);
    num::matrix<double> y_pred = X.dot(w_est);
    double w_err = (w - w_est).sum();
    double y_err = (y - y_pred).sum();
    EXPECT_NEAR(w_err, 0, 0.000001);
    EXPECT_NEAR(y_err, 0, 0.000001);
}

TEST(LinearExact, Fit) {
    num::matrix<double> X = num::randn(10000, 10);
    num::matrix<double> w = num::randn(10, 1);
    num::matrix<double> y = X.dot(w);
    num::matrix<double> w_est = num::exact_linear_model(X, y);
    num::matrix<double> y_pred = X.dot(w_est);
    double w_err = (w - w_est).sum();
    double y_err = (y - y_pred).sum();
    EXPECT_NEAR(w_err, 0, 0.000001);
    EXPECT_NEAR(y_err, 0, 0.000001);
}
