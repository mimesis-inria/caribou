#pragma once

#include <numeric>
#include <string>
#include <list>
#include <Eigen/Core>
#include <gtest/gtest.h>

template <typename Derived1, typename Derived2>
::testing::AssertionResult CmpHelperMatrixEQ(const char* lhs_expression,
                                             const char* rhs_expression,
                                             const Eigen::MatrixBase<Derived1> & lhs_value,
                                             const Eigen::MatrixBase<Derived2> & rhs_value) {
    auto max_error = static_cast<double>(((lhs_value) - (rhs_value)).cwiseAbs().maxCoeff());

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
    auto max_error = static_cast<double>(((lhs_value) - (rhs_value)).cwiseAbs().maxCoeff());

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

template<template<typename, typename...> typename H, typename T1, typename T2, typename... Ts>
::testing::AssertionResult list_are_equals_predicate (const char* lhs_expression,
                                                      const char* rhs_expression,
                                                      const std::list<T1> & lhs_value,
                                                      const H<T2, Ts...> & rhs_value)
{
    typename std::list<T1>::const_iterator it1;
    typename H<T2, Ts...>::const_iterator it2;

    if (lhs_value.size() == rhs_value.size()) {
        bool all_equals = true;
        for (it1 = lhs_value.begin(), it2 = rhs_value.begin(); it1 != lhs_value.end(); ++it1, ++it2) {
            if (*it1 != *it2) {
                all_equals = false;
                break;
            }
        }
        if (all_equals) {
            return ::testing::AssertionSuccess();
        }
    }

    auto lhs_string = (lhs_value.size() == 0) ? "" :
                      std::accumulate(std::next(lhs_value.begin()), lhs_value.end(), std::to_string(*lhs_value.begin()),
                      [](std::string a, const T1 & b) {
                          return a + ", " + std::to_string(b);
                      })
    ;
    auto rhs_string = (rhs_value.size() == 0) ? "" :
                      std::accumulate(std::next(rhs_value.begin()), rhs_value.end(), std::to_string(*rhs_value.begin()),
                      [](std::string a, const T2 & b) {
                          return a + ", " + std::to_string(b);
                      })
    ;

    return ::testing::internal::EqFailure(
            lhs_expression,
            rhs_expression,
            "[" + lhs_string + "]",
            "[" + rhs_string + "]",
            false);
}

#define EXPECT_LIST_EQUAL(val1, val2)\
  EXPECT_PRED_FORMAT2(::list_are_equals_predicate, val1, val2)

extern std::string executable_directory_path;