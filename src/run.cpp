#include <iostream>

#include "io.hpp"
#include "linalg.hpp"
#include "linear.hpp"
#include "matrix.hpp"

void test_matrix() {
    num::matrix<double> empty_matrix;
    num::matrix<double> matrix(2, 2, 0);
    num::matrix<double> id(2, 2, 0);
    id.set(0, 0, 1);
    id.set(1, 1, 1);
    matrix.set(0, 0, 1);
    matrix.set(0, 1, 2);
    matrix.set(1, 0, 3);
    matrix.set(1, 1, 4);
    matrix.print();
    matrix.transpose().print();
    id.print();
    matrix.dot(id).print();
    id.dot(matrix).print();

    (1 + id + 1).print();
    (id + id).print();
    (id - 10).print();
    (id / 2).print();
    (id * 100).print();
}

void test_sgd() {
    num::matrix<double> X = num::randn(10000, 10);
    num::matrix<double> w = num::randn(10, 1);
    num::matrix<double> y = X.dot(w);
    num::matrix<double> w_est = num::gradient_linear_model(X, y);

    num::matrix<double> y_pred = X.dot(w_est);
    double w_err = (w - w_est).sum();
    double y_err = (y - y_pred).sum();
    w_est.print();
    std::cout << "w_err: " << w_err << " y_err: " << y_err << std::endl;
}

void test_lstsq() {
    num::matrix<double> X = num::randn(10000, 10);
    num::matrix<double> w = num::randn(10, 1);
    num::matrix<double> y = X.dot(w);
    num::matrix<double> w_est = num::exact_linear_model(X, y);

    num::matrix<double> y_pred = X.dot(w_est);
    double w_err = (w - w_est).sum();
    double y_err = (y - y_pred).sum();
    std::cout << "w:" << std::endl;
    w.print();
    std::cout << "w_est:" << std::endl;
    w_est.print();
    std::cout << "w_err: " << w_err << " y_err: " << y_err << std::endl;
}

void test_inverse() {
    num::matrix<double> matrix = num::eye<double>(3);
    matrix.set(2, 2, 0.5);
    // num::matrix<double> result = inverse(matrix);
    // result.print();

    num::matrix<double> m2({{1, -2, 5}, {3, 1, 0}, {1, 0, 2}});
    num::matrix<double> r2 = inverse(m2);
    r2.print();
}

void test_io() {
    num::matrix<double> mat = num::read_csv<double>("../data/winequality-red.csv", ';');
    // mat.print();
    std::cout << mat.shape() << std::endl;
}

int main() {
    test_io();
    return 0;
}
