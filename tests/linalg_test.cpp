#include <gtest/gtest.h>

#include "../src/linalg.hpp"

TEST(SwapRows, Values) {
    num::matrix<int> mat({{1, 1}, {2, 2}, {3, 3}});
    num::matrix<int> expected({{2, 2}, {1, 1}, {3, 3}});

    swap_rows<int>(mat[0], mat[1]);

    EXPECT_TRUE((mat == expected).all());
}

TEST(Inverse, Values) {
    num::matrix<double> m1 = num::eye<double>(3);
    m1[2][2] = 0.5;
    num::matrix<double> r1 = inverse(m1);
    num::matrix<double> expected1({{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.5}});
    EXPECT_TRUE((m1 == expected1).all());

    num::matrix<double> m2({{1, -2, 5}, {3, 1, 0}, {1, 0, 2}});
    num::matrix<double> r2 = inverse(m2);
    num::matrix<double> expected2(
        {{0.2222222222222221, 0.44444444444444453, -0.5555555555555558},
         {-0.666666666666667, -0.3333333333333335, 1.6666666666666679},
         {-0.11111111111111116, -0.2222222222222222, 0.7777777777777779}});
    double err = (r2 - expected2).sum();
    EXPECT_NEAR(err, 0, 0.000001);
}
