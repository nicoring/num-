#include <chrono>
#include <iostream>

#include "matrix.hpp"

template <typename T>
uint64_t benchmark_dot_product(const num::matrix<T>& m1, const num::matrix<T>& m2) {
    const auto insert_start = std::chrono::steady_clock::now();
    num::matrix<T> m3 = m1.dot(m2);
    const auto insert_end = std::chrono::steady_clock::now();
    uint64_t time_milli =
        std::chrono::duration_cast<std::chrono::milliseconds>(insert_end - insert_start)
            .count();
    uint64_t ops = m1.shape().rows * m1.shape().cols * m2.shape().cols;
    uint64_t ops_per_millisecond = ops / time_milli;
    return ops_per_millisecond;
}

int main() {
    const num::matrix<double> m1(100000, 1000, 0.0);
    const num::matrix<double> m2(1000, 100, 0.0);
    uint64_t ops_per_millisecond = benchmark_dot_product(m1, m2);
    std::cout << "ops per millisecond: " << ops_per_millisecond << std::endl;
}
