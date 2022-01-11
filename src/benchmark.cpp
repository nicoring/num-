#include <chrono>
#include <iostream>

#include "matrix.hpp"

int main() {
    num::matrix<double> m1(1000, 1000);
    num::matrix<double> m2(1000, 1000);
    const auto insert_start = std::chrono::steady_clock::now();
    num::matrix<double> m3 = m1.dot(m2);
    const auto insert_end = std::chrono::steady_clock::now();
    uint64_t time = std::chrono::duration_cast<std::chrono::milliseconds>(
                        insert_end - insert_start)
                        .count();
    std::cout << time << std::endl;
}
