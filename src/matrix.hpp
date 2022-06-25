#pragma once

#include <iostream>
#include <ostream>
#include <random>
#include <span>
#include <sstream>
#include <vector>

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

template <typename T> class matrix_row {
  private:
    std::span<T> row;

  public:
    matrix_row(std::span<T> row) : row(row){};
    uint32_t length() { return row.size(); };

    const T& operator[](uint32_t col) const { return row[col]; };
    T& operator[](uint32_t col) { return row[col]; };

    matrix_row& operator=(T value) {
        for (T& elem : row) {
            elem = value;
        }
        return *this;
    }

    matrix_row& operator+=(T value) {
        for (T& elem : row) {
            elem += value;
        }
        return *this;
    }

    matrix_row& operator-=(T value) {
        for (T& elem : row) {
            elem -= value;
        }
        return *this;
    }

    matrix_row& operator*=(T value) {
        for (T& elem : row) {
            elem *= value;
        }
        return *this;
    }

    matrix_row& operator/=(T value) {
        for (T& elem : row) {
            elem /= value;
        }
        return *this;
    }
};

template <typename T> class matrix {
  private:
    std::vector<T> data;
    const uint32_t rows;
    const uint32_t cols;

    template <typename LoopFn> void for_each(LoopFn fn) const {
        for (uint32_t row = 0; row < rows; ++row) {
            for (uint32_t col = 0; col < cols; ++col) {
                fn(row, col, get(row, col));
            }
        }
    }

    template <typename U, typename LoopFn> matrix<U> map(LoopFn fn) const {
        matrix<U> new_mat(rows, cols);
        for (uint32_t row = 0; row < rows; ++row) {
            for (uint32_t col = 0; col < cols; ++col) {
                new_mat.set(row, col, fn(row, col, get(row, col)));
            }
        }
        return new_mat;
    }

    void ensure_same_dim(const matrix& other) const {
        if (other.rows != rows || other.cols != cols) {
            throw std::invalid_argument("cols or rows don't match for operation");
        }
    }

  public:
    matrix() : matrix(0, 0){};

    matrix(uint32_t rows, uint32_t cols) : rows(rows), cols(cols), data(rows * cols){};

    matrix(uint32_t rows, uint32_t cols, T fill_value)
        : rows(rows), cols(cols), data(rows * cols, fill_value){};

    matrix(uint32_t rows, uint32_t cols, const std::vector<T>& values)
        : matrix(rows, cols) {
        if (values.size() != rows * cols) {
            throw std::invalid_argument(
                "rows and cols do not match with passed values");
        }
        data.assign(values.begin(), values.end());
    };

    matrix(const std::vector<std::vector<T>>& values)
        : matrix(values.size(), values[0].size()) {
        for (u_int32_t row = 0; row < rows; row++) {
            if (values[row].size() != cols) {
                throw std::invalid_argument("Not all rows of passed values "
                                            "have the same number of columns");
            }
            for (u_int32_t col = 0; col < cols; col++) {
                set(row, col, values[row][col]);
            }
        }
    };

    matrix(const matrix<T>& other) : matrix(other.rows, other.cols) {
        for_each([&](uint32_t row, uint32_t col, T value) {
            set(row, col, other.get(row, col));
        });
    }

    shape_t shape() const { return shape_t{rows, cols}; };

    bool empty() const { return rows == 0 || cols == 0; };

    inline T get(uint32_t row, uint32_t col) const { return data[row * cols + col]; };

    matrix_row<T> get(uint32_t row) {
        uint32_t offset = row * cols;
        return matrix_row(std::span(data.data() + offset, cols));
    };

    const matrix_row<T> operator[](uint32_t row) const { return get(row); };
    matrix_row<T> operator[](uint32_t row) { return get(row); };

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
                "inner dimensions don't match for dot product: " + my_shape + " " +
                other_shape);
        }
        matrix new_mat(rows, other.cols);

// #pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
        for (uint32_t new_row = 0; new_row < rows; new_row++) {
// #pragma GCC ivdep
#pragma clang loop vectorize(enable) interleave(enable)
            for (uint32_t new_col = 0; new_col < other.cols; new_col++) {

// T value = 0;
// #pragma GCC ivdep
#pragma clang loop vectorize_width(4) interleave_count(2)
                for (uint32_t inner = 0; inner < cols; inner++) {
                    // value += get(new_row, inner) * other.get(inner, new_col);
                    const T to_add = get(new_row, inner) * other.get(inner, new_col);
                    const T current = new_mat.get(new_row, new_col);
                    new_mat.set(new_row, new_col, current + to_add);
                }
                // new_mat.set(new_row, new_col, value);
            }
        }
        return new_mat;
    };

    matrix transpose() const {
        matrix<T> new_mat{cols, rows};
        for_each(
            [&](uint32_t row, uint32_t col, T value) { new_mat.set(col, row, value); });
        return new_mat;
    };

    matrix& operator=(const matrix& m) {
        if (this == &m)
            return *this;
        if (m.rows != rows || m.cols != cols) {
            throw std::invalid_argument("cols or rows don't match for assignment");
        }
        data.assign(m.data.begin(), m.data.end());
        return *this;
    }

    // + operators

    matrix operator+(T rhs) const {
        return map<T>([&](uint32_t row, uint32_t col, T value) { return value + rhs; });
    }

    friend matrix operator+(T lhs, const matrix& rhs) {
        return rhs.map<T>(
            [&](uint32_t row, uint32_t col, T value) { return lhs + value; });
    }

    matrix operator+(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<T>([&](uint32_t row, uint32_t col, T value) {
            return value + rhs.get(row, col);
        });
    };

    // - operators

    matrix operator-(T rhs) const {
        return map<T>([&](uint32_t row, uint32_t col, T value) { return value - rhs; });
    }

    friend matrix operator-(T lhs, const matrix& rhs) {
        return rhs.map<T>(
            [&](uint32_t row, uint32_t col, T value) { return lhs - value; });
    }

    matrix operator-(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<T>([&](uint32_t row, uint32_t col, T value) {
            return value - rhs.get(row, col);
        });
    };

    // * operators

    matrix operator*(T rhs) const {
        return map<T>([&](uint32_t row, uint32_t col, T value) { return value * rhs; });
    }

    friend matrix operator*(T lhs, const matrix& rhs) {
        return rhs.map<T>(
            [&](uint32_t row, uint32_t col, T value) { return lhs * value; });
    }

    matrix operator*(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<T>([&](uint32_t row, uint32_t col, T value) {
            return value * rhs.get(row, col);
        });
    };

    // / operators

    matrix operator/(T rhs) const {
        return map<T>([&](uint32_t row, uint32_t col, T value) { return value / rhs; });
    }

    friend matrix operator/(T lhs, const matrix& rhs) {
        return rhs.map<T>(
            [&](uint32_t row, uint32_t col, T value) { return lhs / value; });
    }

    matrix operator/(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<T>([&](uint32_t row, uint32_t col, T value) {
            return value / rhs.get(row, col);
        });
    };

    // == operators

    matrix<bool> operator==(T rhs) const {
        return map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return value == rhs; });
    }

    friend matrix<bool> operator==(T lhs, const matrix& rhs) {
        return rhs.map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return lhs == value; });
    }

    matrix<bool> operator==(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<bool>([&](uint32_t row, uint32_t col, T value) {
            return value == rhs.get(row, col);
        });
    };

    // != operators

    matrix<bool> operator!=(T rhs) const {
        return map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return value != rhs; });
    }

    friend matrix<bool> operator!=(T lhs, const matrix& rhs) {
        return rhs.map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return lhs != value; });
    }

    matrix<bool> operator!=(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<bool>([&](uint32_t row, uint32_t col, T value) {
            return value != rhs.get(row, col);
        });
    };

    // <= operators

    matrix<bool> operator<=(T rhs) const {
        return map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return value <= rhs; });
    }

    friend matrix<bool> operator<=(T lhs, const matrix& rhs) {
        return rhs.map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return lhs <= value; });
    }

    matrix<bool> operator<=(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<bool>([&](uint32_t row, uint32_t col, T value) {
            return value <= rhs.get(row, col);
        });
    };

    // >= operators

    matrix<bool> operator>=(T rhs) const {
        return map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return value >= rhs; });
    }

    friend matrix<bool> operator>=(T lhs, const matrix& rhs) {
        return rhs.map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return lhs >= value; });
    }

    matrix<bool> operator>=(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<bool>([&](uint32_t row, uint32_t col, T value) {
            return value >= rhs.get(row, col);
        });
    };

    // < operators

    matrix<bool> operator<(T rhs) const {
        return map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return value < rhs; });
    }

    friend matrix<bool> operator<(T lhs, const matrix& rhs) {
        return rhs.map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return lhs < value; });
    }

    matrix<bool> operator<(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<bool>([&](uint32_t row, uint32_t col, T value) {
            return value < rhs.get(row, col);
        });
    };

    // > operators

    matrix<bool> operator>(T rhs) const {
        return map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return value > rhs; });
    }

    friend matrix<bool> operator>(T lhs, const matrix& rhs) {
        return rhs.map<bool>(
            [&](uint32_t row, uint32_t col, T value) { return lhs > value; });
    }

    matrix<bool> operator>(const matrix& rhs) const {
        ensure_same_dim(rhs);
        return map<bool>([&](uint32_t row, uint32_t col, T value) {
            return value > rhs.get(row, col);
        });
    };

    // aggregation methods

    T sum() const {
        T res = 0;
        for_each([&](uint32_t row, uint32_t col, T value) { res += value; });
        return res;
    }

    bool all() const {
        bool res = true;
        for_each([&](uint32_t row, uint32_t col, T value) { res &= value; });
        return res;
    }

    bool any() const {
        bool res = false;
        for_each([&](uint32_t row, uint32_t col, T value) { res |= value; });
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

template <typename T> matrix<T> eye(uint32_t size) { return eye<T>(size, size); }

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

matrix<double> randn(uint32_t rows, uint32_t cols, double mean = 0, double stddev = 1) {
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
