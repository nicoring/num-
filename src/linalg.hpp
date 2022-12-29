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

template <typename T> T eigenvalue(const matrix<T>& A, const matrix<T>& v) {
    return (A.dot(v) / v)[0][0];
}

} // namespace num

/*

#https://jeremykun.com/2016/05/16/singular-value-decomposition-part-2-theorem-proof-algorithm/
#https://github.com/j2kun/svd/blob/main/svd.py
def eigenvalue(A, v):
    val = A @ v / v
    return val[0]

def svd_dominant_eigen(A, epsilon=0.01):
    """returns dominant eigenvalue and dominant eigenvector of matrix A"""
    n, m = A.shape
    k=min(n,m)
    v = np.ones(k) / np.sqrt(k)
    if n > m:
        A = A.T @ A
    elif n < m:
        A = A @ A.T

    ev = eigenvalue(A, v)

    while True:
        Av = A@ v
        v_new = Av / np.linalg.norm(Av)
        ev_new = eigenvalue(A, v_new)
        if np.abs(ev - ev_new) < epsilon:
            break

        v = v_new
        ev = ev_new

    return ev_new, v_new

def svd(A, k=None, epsilon=1e-10):
    """returns k dominant eigenvalues and eigenvectors of matrix A"""
    A = np.array(A, dtype=float)
    n, m = A.shape

    svd_so_far = []
    if k is None:
        k = min(n, m)

    for i in range(k):
        matrix_for_1d = A.copy()

        for singular_value, u, v in svd_so_far[:i]:
            matrix_for_1d -= singular_value * np.outer(u, v)

        if n > m:
            _, v = svd_dominant_eigen(matrix_for_1d, epsilon=epsilon)  # next singular
vector u_unnormalized = A @ v sigma = np.linalg.norm(u_unnormalized)  # next singular
value u = u_unnormalized / sigma else:
            _, u = svd_dominant_eigen(matrix_for_1d, epsilon=epsilon)  # next singular
vector v_unnormalized = A.T @ u sigma = np.linalg.norm(v_unnormalized)  # next singular
value v = v_unnormalized / sigma

        svd_so_far.append((sigma, u, v))

    singular_values, us, vs = [np.array(x) for x in zip(*svd_so_far)]
    return singular_values, us.T, vs

*/