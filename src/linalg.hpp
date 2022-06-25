#pragma once

#include "matrix.hpp"

template <typename T>
void eliminate(num::matrix_row<T> source_1, num::matrix_row<T> source_2,
               num::matrix_row<T> target_1, num::matrix_row<T> target_2, uint32_t col,
               T target) {
    T factor = (source_2.get(col) - target) / source_1.get(col);
    const uint32_t length = source_1.length();
    for (uint32_t i = 0; i < length; ++i) {
        source_2.set(i, source_2.get(i) - factor * source_1.get(i));
        target_2.set(i, target_2.get(i) - factor * target_1.get(i));
    }
}

template <typename T> void swap_rows(num::matrix<T>& m, uint32_t row1, uint32_t row2) {
    T tmp;
    for (uint32_t i = 0; i < m.shape().cols; ++i) {
        tmp = m.get(row1, i);
        m.set(row1, i, m.get(row2, i));
        m.set(row2, i, tmp);
    }
}

template <typename T>
void gaussian_elemination(num::matrix<T>& source, num::matrix<T>& target) {
    const uint32_t size = source.shape().rows;

    // build triangluar matrix
    for (uint32_t i = 0; i < size; ++i) {
        if (source.get(i, i) == 0) {
            bool found = false;
            for (uint32_t j = i + 1; j < size; ++j) {
                if (source.get(i, j) != 0) {
                    swap_rows(source, i, j);
                    swap_rows(target, i, j);
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error("Matrix is not invertible");
            }
        }
        for (uint32_t j = i + 1; j < size; ++j) {
            eliminate<T>(source.get(i), source.get(j), target.get(i), target.get(j), i,
                         0);
        }
    }

    for (uint32_t i = size - 1; i > 0; --i) {
        for (int32_t j = i - 1; j >= 0; --j) {
            eliminate<T>(source.get(i), source.get(j), target.get(i), target.get(j), i,
                         0);
        }
    }

    for (uint32_t i = 0; i < size; ++i) {
        eliminate<T>(source.get(i), source.get(i), target.get(i), target.get(i), i, 1);
    };
}

template <typename T> num::matrix<T> inverse(const num::matrix<T> mat) {
    num::shape_t shape = mat.shape();
    if (shape.rows != shape.cols) {
        throw std::invalid_argument("Matrix must be square");
    }
    num::matrix<T> mat_copy(mat);
    num::matrix<T> res = num::eye<T>(shape.rows);
    gaussian_elemination<T>(mat_copy, res);
    return res;
}
