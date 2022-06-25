#pragma once

#include "matrix.hpp"

namespace num {

template <typename T>
void eliminate(matrix_row<T> source_1, matrix_row<T> source_2, matrix_row<T> target_1,
               matrix_row<T> target_2, uint32_t col, T target) {
    T factor = (source_2[col] - target) / source_1[col];
    const uint32_t length = source_1.length();
    for (uint32_t i = 0; i < length; ++i) {
        source_2[i] = source_2[i] - factor * source_1[i];
        target_2[i] = target_2[i] - factor * target_1[i];
    }
}

template <typename T> void swap_rows(matrix_row<T> row1, matrix_row<T> row2) {
    T tmp;
    for (uint32_t i = 0; i < row1.length(); ++i) {
        tmp = row1[i];
        row1[i] = row2[i];
        row2[i] = tmp;
    }
}

template <typename T> void gaussian_elemination(matrix<T>& source, matrix<T>& target) {
    const uint32_t size = source.shape().rows;

    // build triangluar matrix
    for (uint32_t i = 0; i < size; ++i) {
        if (source[i][i] == 0) {
            bool found = false;
            for (uint32_t j = i + 1; j < size; ++j) {
                if (source[i][j] != 0) {
                    swap_rows(source[i], source[j]);
                    swap_rows(target[i], target[j]);
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error("Matrix is not invertible");
            }
        }
        for (uint32_t j = i + 1; j < size; ++j) {
            eliminate<T>(source[i], source[j], target[i], target[j], i, 0);
        }
    }

    for (uint32_t i = size - 1; i > 0; --i) {
        for (int32_t j = i - 1; j >= 0; --j) {
            eliminate<T>(source[i], source[j], target[i], target[j], i, 0);
        }
    }

    for (uint32_t i = 0; i < size; ++i) {
        eliminate<T>(source[i], source[i], target[i], target[i], i, 1);
    };
}

template <typename T> matrix<T> inverse(const matrix<T> mat) {
    shape_t shape = mat.shape();
    if (shape.rows != shape.cols) {
        throw std::invalid_argument("Matrix must be square");
    }
    matrix<T> mat_copy(mat);
    matrix<T> res = eye<T>(shape.rows);
    gaussian_elemination<T>(mat_copy, res);
    return res;
}
} // namespace num
