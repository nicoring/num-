#include <gtest/gtest.h>

#include "../src/linear.hpp"

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
    EXPECT_THROW(num::matrix<int> mat({{1, 2}, {4, 5, 6}}),
                 std::invalid_argument);
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

TEST(AggregationMethods, Sum) {
    num::matrix<int> mat(2, 2, 2);
    EXPECT_EQ(mat.sum(), 2 * 2 * 2);
    num::matrix<int> empty;
    EXPECT_EQ(empty.sum(), 0);
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

TEST(ComparisonOperators, Eq) {

}