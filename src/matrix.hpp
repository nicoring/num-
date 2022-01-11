#include <cstddef>
#include <iostream>
#include <ostream>
#include <random>
#include <sstream>
#include <vector>

#ifndef MATRIX_H
#define MATRIX_H

namespace num {

struct shape_t {
    uint32_t rows;
    uint32_t cols;
};

std::ostream& operator<<(std::ostream& out, const shape_t& shape) {
    out << "(" << shape.rows << ", " << shape.cols << ")";
    return out;
}

std::string to_string(const shape_t& shape) {
    std::ostringstream ss;
    ss << shape;
    return std::move(ss).str();
}

template <typename T> class matrix {
  private:
    std::vector<T> data;
    const uint32_t rows;
    const uint32_t cols;

  public:
    matrix() : matrix(0, 0){};

    matrix(uint32_t rows, uint32_t cols)
        : rows(rows), cols(cols), data(rows * cols){};

    matrix(uint32_t rows, uint32_t cols, T fill_value)
        : rows(rows), cols(cols), data(rows * cols, fill_value){};

    shape_t shape() const { return shape_t{rows, cols}; };

    bool empty() const { return rows == 0 && cols == 0; };

    inline T get(uint32_t row, uint32_t col) const {
        return data[row * cols + col];
    };

    inline void set(uint32_t row, uint32_t col, T value) {
        data[row * cols + col] = value;
    };

    void print() const {
        for (uint32_t row = 0; row < rows; row++) {
            for (uint32_t col = 0; col < cols; col++) {
                std::cout << get(row, col) << " ";
            }
            std::cout << std::endl;
        }
    };

    matrix dot(matrix other) const {
        if (cols != other.rows) {
            std::string my_shape = to_string(shape());
            std::string other_shape = to_string(other.shape());
            throw std::invalid_argument(
                "inner dimensions don't match for dot product: " + my_shape +
                " " + other_shape);
        }
        matrix new_mat(rows, other.cols);
        for (uint32_t new_row = 0; new_row < rows; new_row++) {
            for (uint32_t new_col = 0; new_col < other.cols; new_col++) {

                T value = 0;
                for (uint32_t inner = 0; inner < cols; inner++) {
                    value += get(new_row, inner) * other.get(inner, new_col);
                }
                new_mat.set(new_row, new_col, value);
            }
        }
        return new_mat;
    };

    matrix transpose() const {
        matrix new_mat(cols, rows);
        for (uint32_t row = 0; row < rows; row++) {
            for (uint32_t col = 0; col < cols; col++) {
                new_mat.set(col, row, get(row, col));
            }
        }
        return new_mat;
    };

    matrix& operator=(const matrix& m) {
        if (this == &m)
            return *this;
        if (m.rows != rows || m.cols != cols) {
            throw std::invalid_argument(
                "cols or rows don't match for assignment");
        }
        data.assign(m.data.begin(), m.data.end());
        return *this;
    }

#define add_operators(op)                                                      \
    matrix operator op(T rhs) const {                                          \
        matrix new_mat(rows, cols);                                            \
        for (uint32_t row = 0; row < rows; row++) {                            \
            for (uint32_t col = 0; col < cols; col++) {                        \
                new_mat.set(row, col, get(row, col) op rhs);                   \
            }                                                                  \
        }                                                                      \
        return new_mat;                                                        \
    };                                                                         \
                                                                               \
    friend matrix operator op(T lhs, const matrix& rhs) {                      \
        matrix<T> new_mat(rhs.rows, rhs.cols);                                 \
        for (uint32_t row = 0; row < rhs.rows; row++) {                        \
            for (uint32_t col = 0; col < rhs.cols; col++) {                    \
                new_mat.set(row, col, lhs op rhs.get(row, col));               \
            }                                                                  \
        }                                                                      \
        return new_mat;                                                        \
    };                                                                         \
                                                                               \
    matrix operator op(const matrix& rhs) const {                              \
        if (rhs.rows != rows || rhs.cols != cols) {                            \
            throw std::invalid_argument(                                       \
                "cols or rows don't match for operation ");                    \
        }                                                                      \
        matrix new_mat(rows, cols);                                            \
        for (uint32_t row = 0; row < rows; row++) {                            \
            for (uint32_t col = 0; col < cols; col++) {                        \
                new_mat.set(row, col, get(row, col) op rhs.get(row, col));     \
            }                                                                  \
        }                                                                      \
        return new_mat;                                                        \
    };

    add_operators(+);
    add_operators(-);
    add_operators(*);
    add_operators(/);
#undef add_operators

    T sum() const {
        T res = 0;
        for (uint32_t row = 0; row < rows; row++) {
            for (uint32_t col = 0; col < cols; col++) {
                res += get(row, col);
            }
        }
        return res;
    }
};

template <typename T> matrix<T> zeros(uint32_t rows, uint32_t cols) {
    return matrix<T>(rows, cols, 0);
}

template <typename T> matrix<T> ones(uint32_t rows, uint32_t cols) {
    return matrix<T>(rows, cols, 1);
}

template <typename T> matrix<T> eye(uint32_t rows, uint32_t cols) {
    matrix<T> mat(rows, cols, 0);
    uint32_t min_count = std::min(rows, cols);
    for (uint32_t i = 0; i < min_count; i++) {
        mat.set(i, i, 1);
    }
    return mat;
}

template <typename T> matrix<T> eye(uint32_t size) {
    return eye<T>(size, size);
}

matrix<double> rand(uint32_t rows, uint32_t cols) {
    matrix<double> mat(rows, cols);

    // Will be used to obtain a seed for the random number engine
    std::random_device rd;
    // Standard mersenne_twister_engine seeded with rd()
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform(0.0, 1.0);

    for (uint32_t row = 0; row < rows; row++) {
        for (uint32_t col = 0; col < cols; col++) {
            mat.set(row, col, uniform(gen));
        }
    }

    return mat;
}

matrix<double> randn(uint32_t rows, uint32_t cols, double mean = 0,
                     double stddev = 1) {
    matrix<double> mat(rows, cols);

    // Will be used to obtain a seed for the random number engine
    std::random_device rd;
    // Standard mersenne_twister_engine seeded with rd()
    std::mt19937 gen(rd());
    std::normal_distribution<> normal(mean, stddev);

    for (uint32_t row = 0; row < rows; row++) {
        for (uint32_t col = 0; col < cols; col++) {
            mat.set(row, col, normal(gen));
        }
    }

    return mat;
}

} // namespace num

#endif
