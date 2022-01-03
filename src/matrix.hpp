#include <vector>
#include <cstddef>
#include <iostream>

struct shape_t {
    size_t rows, cols;
};

template<class T>
class mat {
    std::vector<T> data;
    const size_t rows;
    const size_t cols;
    
    public:
    mat(): rows(0), cols(0) {};
    ~mat() {};
    mat(size_t rows, size_t cols): rows(rows), cols(cols) {
        data = std::vector<T>(rows * cols);
    }
    mat(size_t rows, size_t cols, T fill_value): rows(rows), cols(cols) {
        data = std::vector<T>(rows * cols, fill_value);
    }

    shape_t shape() { return shape_t{rows, cols}; };
    bool empty() {return rows == 0 and cols == 0; };

    inline T get(size_t row, size_t col) const {
        return data[rows * row + col];
    }

    inline void set(size_t row, size_t col, T value) {
        data[rows * row + col] = value;
    }

    void print() const {
        for (size_t row = 0; row < rows; row++) {
            for (size_t col = 0; col < cols; col++) {
                std::cout << get(row, col) << " ";
            }
            std::cout << std::endl;
        }
    }

    mat dot(mat other) const {
        if (not (cols == other.rows)) {
            throw std::invalid_argument("inner dimensions don't match for dot product");
        }
        mat new_mat = mat(rows, other.cols);
        for (size_t new_row = 0; new_row < rows; new_row++) {
            for (size_t new_col = 0; new_col < other.cols; new_col++) {
                T value = 0;
                for (size_t inner = 0; inner < cols; inner++) {
                    value += get(new_row, inner) * other.get(inner, new_col);
                }
                new_mat.set(new_row, new_col, value);
            }
        }
        return new_mat;
    };

    mat transpose() const {
        mat new_mat = mat(cols, rows);
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                new_mat.set(col, row, get(row, col));
            }
        }
        return new_mat;
    };

    mat operator+(T to_add) const {
        mat new_mat = mat(cols, rows);
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                new_mat.set(row, col, get(row, col) + to_add);
            }
        }
        return new_mat;
    }

    mat operator+(mat other) const {
        mat new_mat = mat(cols, rows);
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                new_mat.set(row, col, get(row, col) + other.get(row, col));
            }
        }
        return new_mat;
    }

    friend mat operator+(T to_add, mat m) { return m + to_add; }
};