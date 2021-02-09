#pragma once

#include <string>
#include <list>
#include <Eigen/Core>
#include <gtest/gtest.h>

template <typename Derived1, typename Derived2>
::testing::AssertionResult CmpHelperMatrixEQ(const char* lhs_expression,
                                             const char* rhs_expression,
                                             const Eigen::MatrixBase<Derived1> & lhs_value,
                                             const Eigen::MatrixBase<Derived2> & rhs_value) {
    auto max_error = ((lhs_value) - (rhs_value)).cwiseAbs().maxCoeff();

    if (::testing::internal::FloatingPoint<double>(0).AlmostEquals(::testing::internal::FloatingPoint<double>(max_error))) {
        return ::testing::AssertionSuccess();
    }

    Eigen::IOFormat clean(std::numeric_limits<double>::digits10 + 2, 0, " ", "\n", "    |", "|");

    ::std::stringstream lhs_ss;
    lhs_ss << lhs_value.format(clean);
    ::std::stringstream rhs_ss;
    rhs_ss << rhs_value.format(clean);


    return ::testing::internal::EqFailure(
        lhs_expression,
        rhs_expression,
        ::testing::internal::StringStreamToString(&lhs_ss),
        ::testing::internal::StringStreamToString(&rhs_ss),
        false);
}

#define EXPECT_MATRIX_EQUAL(val1, val2)\
  EXPECT_PRED_FORMAT2(::CmpHelperMatrixEQ, val1, val2)

template <typename Derived1, typename Derived2>
::testing::AssertionResult CmpHelperMatrixEQ(const char* lhs_expression,
                                             const char* rhs_expression,
                                             const char* /*abs_error_expr*/,
                                             const Eigen::MatrixBase<Derived1> & lhs_value,
                                             const Eigen::MatrixBase<Derived2> & rhs_value,
                                             double abs_error) {
    auto max_error = ((lhs_value) - (rhs_value)).cwiseAbs().maxCoeff();

    if (max_error < abs_error) {
        return ::testing::AssertionSuccess();
    }

    Eigen::IOFormat clean(static_cast<int>(std::abs(std::log10(abs_error))), 0, " ", "\n", "    |", "|");

    return ::testing::AssertionFailure()
        << "The difference between " << lhs_expression << " and " << rhs_expression
        << " is " << max_error << ", which exceeds " << abs_error << ", where\n"
        << lhs_expression << " evaluates to \n" << lhs_value.format(clean) << "\n"
        << rhs_expression << " evaluates to \n" << rhs_value.format(clean) << "\n";
}

#define EXPECT_MATRIX_NEAR(val1, val2, abs_error)\
  EXPECT_PRED_FORMAT3(::CmpHelperMatrixEQ, val1, val2, abs_error)

template<template<typename, typename...> typename H, typename T, typename... Ts>
bool list_are_equals (const std::list<T> & l1, const H<T, Ts...> & l2)
{
    if (l1.size() != l2.size())
        return false;

    typename std::list<T>::const_iterator it1;
    typename H<T, Ts...>::const_iterator it2;

    for (it1 = l1.begin(), it2 = l2.begin(); it1 != l1.end(); ++it1, ++it2)
        if (*it1 != *it2)
            return false;
    return true;
}

extern std::string executable_directory_path;