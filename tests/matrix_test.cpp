#include <gtest/gtest.h>

#include "../src/matrix.hpp"

TEST(Initialization, Empty) {
    num::matrix<double> empty_matrix;
    EXPECT_TRUE(empty_matrix.empty());
    EXPECT_EQ(empty_matrix.shape().rows, 0);
    EXPECT_EQ(empty_matrix.shape().cols, 0);
}

TEST(Initialization, NonEmpty) {
    num::matrix<double> mat(10, 20);
    EXPECT_FALSE(mat.empty());
    EXPECT_EQ(mat.shape().rows, 10);
    EXPECT_EQ(mat.shape().cols, 20);
}

TEST(Initialization, WithValue) {
    num::matrix<double> mat(2, 2, 1.);
    EXPECT_FALSE(mat.empty());
    EXPECT_EQ(mat.shape().rows, 2);
    EXPECT_EQ(mat.shape().cols, 2);
    EXPECT_DOUBLE_EQ(mat.get(0, 0), 1.);
    EXPECT_DOUBLE_EQ(mat.get(0, 1), 1.);
    EXPECT_DOUBLE_EQ(mat.get(1, 0), 1.);
    EXPECT_DOUBLE_EQ(mat.get(1, 1), 1.);
}

TEST(Initialization, WithValuesVector1) {
    num::matrix<int> mat(2, 3, {1, 2, 3, 4, 5, 6});
    EXPECT_FALSE(mat.empty());
    EXPECT_EQ(mat.shape().rows, 2);
    EXPECT_EQ(mat.shape().cols, 3);
    EXPECT_EQ(mat.get(0, 0), 1);
    EXPECT_EQ(mat.get(0, 1), 2);
    EXPECT_EQ(mat.get(0, 2), 3);
    EXPECT_EQ(mat.get(1, 0), 4);
    EXPECT_EQ(mat.get(1, 1), 5);
    EXPECT_EQ(mat.get(1, 2), 6);
}

TEST(Initialization, WithValuesVector2) {
    num::matrix<int> mat(3, 2, {1, 2, 3, 4, 5, 6});
    EXPECT_FALSE(mat.empty());
    EXPECT_EQ(mat.shape().rows, 3);
    EXPECT_EQ(mat.shape().cols, 2);
    EXPECT_EQ(mat.get(0, 0), 1);
    EXPECT_EQ(mat.get(0, 1), 2);
    EXPECT_EQ(mat.get(1, 0), 3);
    EXPECT_EQ(mat.get(1, 1), 4);
    EXPECT_EQ(mat.get(2, 0), 5);
    EXPECT_EQ(mat.get(2, 1), 6);
}

TEST(Initialization, WrongValuesVector) {
    EXPECT_THROW(num::matrix<int> mat(1, 100, {{1, 2, 3, 4, 5, 6}}),
                 std::invalid_argument);
}

TEST(Initialization, WithValuesVectorOfVectors) {
    num::matrix<int> mat({{1, 2, 3}, {4, 5, 6}});
    EXPECT_FALSE(mat.empty());
    EXPECT_EQ(mat.shape().rows, 2);
    EXPECT_EQ(mat.shape().cols, 3);
    EXPECT_EQ(mat.get(0, 0), 1);
    EXPECT_EQ(mat.get(0, 1), 2);
    EXPECT_EQ(mat.get(0, 2), 3);
    EXPECT_EQ(mat.get(1, 0), 4);
    EXPECT_EQ(mat.get(1, 1), 5);
    EXPECT_EQ(mat.get(1, 2), 6);
}

TEST(Initialization, WrongValuesVectorOfVectors) {
    EXPECT_THROW(num::matrix<int> mat({{1, 2}, {4, 5, 6}}), std::invalid_argument);
}

TEST(Initialization, RandUniform) {
    num::matrix<double> m = num::rand(1000, 1000);
    double mean = m.sum() / (double)m.size();
    EXPECT_NEAR(mean, 0.5, 0.01);
}

TEST(Initialization, RandNormal) {
    num::matrix<double> m = num::randn(1000, 1000);
    double mean = m.sum() / (double)m.size();
    EXPECT_NEAR(mean, 0, 0.01);
}

TEST(Transpose, values) {
    num::matrix<int> mat({{1, 2, 3}, {4, 5, 6}});
    num::matrix<int> expected({{1, 4}, {2, 5}, {3, 6}});

    num::matrix<int> mat_transposed = mat.transpose();

    EXPECT_TRUE((mat_transposed == expected).all());
}

TEST(Operators, PlusScalarInt) {
    num::matrix<int> mat(2, 2, 0);
    auto mat2 = mat + 2;
    EXPECT_EQ(mat2.get(0, 0), 2);
    auto mat3 = mat + 2.0;
    EXPECT_EQ(mat3.get(0, 0), 2.0);
    auto mat4 = 2 + mat;
    EXPECT_EQ(mat4.get(0, 0), 2);
    auto mat5 = 2.0 + mat;
    EXPECT_EQ(mat5.get(0, 0), 2);
    auto mat6 = 1 + mat + 2.0;
    EXPECT_EQ(mat6.get(0, 0), 3);
}

TEST(Operators, PlusScalarFloat) {
    num::matrix<float> mat(2, 2, 0.0);
    auto mat2 = mat + 2;
    EXPECT_FLOAT_EQ(mat2.get(0, 0), 2.0);
    auto mat3 = mat + 2.5;
    EXPECT_FLOAT_EQ(mat3.get(0, 0), 2.5);
    auto mat4 = 2 + mat;
    EXPECT_FLOAT_EQ(mat4.get(0, 0), 2.0);
    auto mat5 = 2.2 + mat;
    EXPECT_FLOAT_EQ(mat5.get(0, 0), 2.2);
    auto mat6 = 1 + mat + 2.4;
    EXPECT_FLOAT_EQ(mat6.get(0, 0), 3.4);
}

TEST(Operators, DivideInt) {
    num::matrix<int> mat(2, 2, 10);
    auto mat2 = mat / 2;
    EXPECT_EQ(mat2.get(0, 0), 5);
    auto mat3 = 1 / mat;
    EXPECT_EQ(mat3.get(0, 0), 0);
}

TEST(Operators, DivideFloat) {
    num::matrix<float> mat(2, 2, 2.0);
    auto mat2 = mat / 10;
    EXPECT_FLOAT_EQ(mat2.get(0, 0), 0.2);
    auto mat3 = 1 / mat;
    EXPECT_FLOAT_EQ(mat3.get(0, 0), 0.5);
    auto mat4 = mat / 0.5;
    EXPECT_FLOAT_EQ(mat4.get(0, 0), 4);
}

TEST(AggregationMethods, SumInt) {
    num::matrix<int> mat(2, 2, 2);
    EXPECT_EQ(mat.sum(), 8.0);
    num::matrix<int> empty;
    EXPECT_EQ(empty.sum(), 0);
}

TEST(AggregationMethods, SumSmallDouble) {
    num::matrix<double> mat(2, 2, 0.2);
    EXPECT_EQ(mat.sum(), 0.8);
    num::matrix<double> empty;
    EXPECT_EQ(empty.sum(), 0.0);
}

TEST(AggregationMethods, All) {
    num::matrix<bool> mat(2, 2, false);
    EXPECT_FALSE(mat.all());
    num::matrix<int> empty;
    EXPECT_TRUE(empty.all());
    num::matrix<bool> mat2(2, 2, true);
    EXPECT_TRUE(mat2.all());
    mat2.set(0, 0, false);
    EXPECT_FALSE(mat2.all());
}

TEST(AggregationMethods, Any) {
    num::matrix<bool> mat(2, 2, false);
    EXPECT_FALSE(mat.any());
    mat.set(0, 0, true);
    EXPECT_TRUE(mat.any());
    num::matrix<int> empty;
    EXPECT_FALSE(empty.any());
}

TEST(AggregationMethods, Max) {
    num::matrix<int> mat({{1, 2, 3}, {100, 5, 6}});
    EXPECT_EQ(mat.max(), 100);
}

TEST(AggregationMethods, Min) {
    num::matrix<int> mat({{1, 2, 3}, {0, 5, 6}});
    EXPECT_EQ(mat.min(), 0);
}

TEST(AggregationMethods, MinMaxThrow) {
    num::matrix<int> mat;
    EXPECT_THROW(mat.min(), std::invalid_argument);
    EXPECT_THROW(mat.max(), std::invalid_argument);
}

TEST(ComparisonOperators, Eq) {
    num::matrix<int> mat1(2, 2, 2);
    num::matrix<int> mat2(2, 2, 2);
    auto res = mat1 == mat2;
    EXPECT_TRUE(res.get(0, 0));
    EXPECT_TRUE(res.all());

    mat2.set(0, 0, 42);
    auto res2 = mat1 == mat2;
    EXPECT_FALSE(res2.get(0, 0));
    EXPECT_TRUE(res2.get(0, 1));
    EXPECT_FALSE(res2.all());
}

TEST(ComparisonOperators, UnEq) {
    num::matrix<int> mat1(2, 2, 1);
    num::matrix<int> mat2(2, 2, 2);
    auto res = mat1 != mat2;
    EXPECT_TRUE(res.get(0, 0));
    EXPECT_TRUE(res.all());

    mat2.set(0, 0, 1);
    auto res2 = mat1 != mat2;
    EXPECT_FALSE(res2.get(0, 0));
    EXPECT_TRUE(res2.get(0, 1));
    EXPECT_FALSE(res2.all());
}

TEST(ComparisonOperators, GreaterEq) {
    num::matrix<int> mat1(2, 2, 2);
    num::matrix<int> mat2(2, 2, 1);
    auto res = mat1 >= mat2;
    EXPECT_TRUE(res.all());
    mat1.set(0, 0, 0);
    mat1.set(0, 1, 1);
    auto res2 = mat1 >= mat2;
    EXPECT_FALSE(res2.get(0, 0));
    EXPECT_TRUE(res2.get(0, 1));
    EXPECT_FALSE(res2.all());
}

TEST(RowAccess, GetRow) {
    num::matrix<int> mat(3, 2, 2);
    num::matrix<int> expected({{2, 2}, {4, 4}, {2, 2}});
    num::matrix<int> expected2({{2, 2}, {4, 4}, {3, 3}});

    num::matrix_row<int> row = mat.get(1);
    row += 2;
    EXPECT_TRUE((mat == expected).all());

    mat[2] += 1;
    EXPECT_TRUE((mat == expected2).all());

    EXPECT_EQ(mat.get(1, 0), 4);

    row[0] += 1;
    EXPECT_EQ(mat.get(1, 0), 5);

    int value = row[0];
    value += 1;
    EXPECT_EQ(mat.get(1, 0), 5);
}

TEST(RowAccess, Assignment) {
    num::matrix<int> mat(3, 2, 2);

    mat[1] = 4;

    num::matrix<int> expected({{2, 2}, {4, 4}, {2, 2}});

    EXPECT_TRUE((mat == expected).all());
}

TEST(Access, Retrieval) {
    num::matrix<int> mat(3, 2, 2);
    int value = mat[0][1];
    EXPECT_EQ(mat[0][0], 2);
    EXPECT_EQ(value, 2);
}

TEST(Access, InplaceEdit) {
    num::matrix<int> mat(3, 2, 2);

    mat[0][0] += 1;
    EXPECT_EQ(mat.get(0, 0), 3);
}

TEST(Access, Assignment) {
    num::matrix<int> mat(3, 2, 2);

    mat[0][0] = 5;
    EXPECT_EQ(mat.get(0, 0), 5);
}

TEST(Stacking, hstack) {
    const num::matrix<int> a({{1, 2}, {5, 6}});
    const num::matrix<int> b({{3, 4}, {7, 8}});
    num::matrix<int> expected({{1, 2, 3, 4}, {5, 6, 7, 8}});

    const num::matrix<int> result = num::hstack(a, b);

    EXPECT_TRUE((result == expected).all());
}