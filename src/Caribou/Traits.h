#ifndef CARIBOU_TRAITS_H
#define CARIBOU_TRAITS_H

#include <type_traits>

namespace caribou {

/**
 * REQUIRES(condition)
 *
 * Same thing as std::enable_if, just easier to read for the user.
 */
# define REQUIRES(...)                                      \
  typename std::enable_if<(__VA_ARGS__), bool>::type = true


/**
 * is_detect
 * Detects if a given type implements a given method.
 *
 * @example
 * \code{.cpp}
 * // Detects if a given type implements the foo() method
 * template <class T>
 * using has_foo_t = decltype( std::declval<T>().foo() );
 *
 * struct A {void foo() {}};
 *
 * constexpr bool has_foo = is_detected<has_foo_t, A>::value;
 * \endcode
 *
 */
namespace detail
{
template<template<class...> class Trait, class Enabler, class... Args>
struct is_detected : std::false_type
{
};

template<template<class...> class Trait, class... Args>
struct is_detected<Trait, std::void_t < Trait<Args...>>, Args...> : std::true_type
{
};
}

template<template<class...> class Trait, class... Args>
using is_detected = typename detail::is_detected<Trait, void, Args...>::type;


/**
 * has_subscript
 * Detects if a given type implements the subscript [] operator.
 *
 * @example
 * \code{.cpp}
 * // Detects if a given type implements the subscript operator
 *
 * struct A {int operator[](const std::string & key) {...} };
 *
 * constexpr bool has_string_subscript = has_subscript<A, std::string>::value; // true
 * constexpr bool has_double_subscript = has_subscript<A, double>::value; // false
 * \endcode
 */
template<class T, class Index>
using subscript_t = decltype(std::declval<T>()[std::declval<Index>()]);

template<class T, class Index>
using has_subscript = is_detected<subscript_t, T, Index>;

/**
 * is_vector
 * Detects if a given type is a vector-like class (ie. implements the subscript [] operator with an int as parameter).
 *
 * @example
 * \code{.cpp}
 * // Detects if a given type is a vector-like class
 *
 * struct A {int operator[](const std::string & key) {...} };
 *
 * constexpr bool array_is_a_vector = is_vector<std::array<double, 3>>::value; // true
 * constexpr bool a_is_a_vector = is_vector<A>::value; // false
 * \endcode
 */

template<class T>
using is_vector = is_detected<subscript_t, T, std::size_t>;

template<class T>
inline constexpr bool is_vector_v = is_vector<T>::value;
/**
 * remove_extent_of_vector
 * Similar to std::remove_extent with c-arrays, this will return the type contained in a vector-like class
 * (see caribou::is_vector).
 *
 * @example
 * \code{.cpp}
 * // Detects if the type contained in a vector is numeric
 *
 * constexpr bool contains_numeric = std::is_arithmetic_v< remove_extent_of_vector<std::vector<int>>::type >; // true
 * constexpr bool contains_numeric = std::is_arithmetic_v< remove_extent_of_vector<std::vector<string>>::type >; // false
 * \endcode
 */
template <class T>
struct remove_extent_of_vector {
    static_assert(is_vector<T>::value, "T must be a vector (eg. implements the T[Index] operator).");
    using type = typename std::remove_reference<decltype(std::declval<T>()[std::declval<std::size_t>()])>::type ;
};

template <class T>
using remove_extent_of_vector_t = typename remove_extent_of_vector<T>::type;

} // namespace caribou
#endif //CARIBOU_TRAITS_H
