#ifndef CARIBOU_TRAITS_H
#define CARIBOU_TRAITS_H

#include <type_traits>

/**
 * REQUIRES(condition)
 * Usage for function template parameter
 *
 * Same thing as std::enable_if, just easier to read for the user.
 *
 * Example:
 * \code{.cpp}
 * #include <Caribou/Traits.h>
 *
 * template <typename Scalar, REQUIRES(std::is_integral_v<Scalar>)>
 * Scalar add(const Scalar & x, const Scalar & y) {
 *     // Integer addition
 *     return x + y;
 * }
 *
 * void main() {
 *    int result = add(1, 1); // Ok, returns 2
 *    // float result = add(1.0, 1.0); // Compile error, the function isn't define for non-integer types
 * }
 * \endcode
 */
# define REQUIRES(...)                                      \
  typename std::enable_if<(__VA_ARGS__), bool>::type = true

/**
 * CLASS_REQUIRES(condition)
 *
 * Same thing as std::enable_if, just easier to read for the user.
 *
 * Example:
 * \code{.cpp}
 * #include <iostream>
 * #include <Caribou/Traits.h>
 *
 * template <typename T, typename Enable = void>
 * class MyObject {
 *     void print() const {std::cout << "Generic object\n";}
 * };
 *
 * template <typename T>
 * class MyObject<T, CLASS_REQUIRES(std::is_integral_v<Scalar>)> {
 *     void print() const {std::cout << "Integral object\n";}
 * };
 *
 * void main() {
 *    MyObject<float> o1; o1.print(); // Prints "Generic object"
 *    MyObject<int>   o2; o2.print(); // Prints "Integral object"
 * }
 * \endcode
 */
# define CLASS_REQUIRES(...)                                      \
  typename std::enable_if<(__VA_ARGS__)>::type

namespace caribou::internal {
/**
 * Detector infrastructure
 * http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/n4502.pdf
 */
// primary template handles all types not supporting the archetypal Op:
        template<class Default, class// always void; supplied externally
                , template<class...> class Op, class... Args
        >
        struct detector {
            using value_t = std::false_type;
            using type = Default;
        };

// the specialization recognizes and handles only types supporting Op:
        template<class Default, template<class...> class Op, class... Args
        >
        struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> {
            using value_t = std::true_type;
            using type = Op<Args...>;
        };

// nonesuch
        struct nonesuch {
            nonesuch() = delete;

            ~nonesuch() = delete;

            nonesuch(nonesuch const &) = delete;

            void operator=(nonesuch const &) = delete;
        };

// is_detected
        template<template<class...> class Op, class... Args>
        using is_detected = typename detector<nonesuch, void, Op, Args...>::value_t;

        template<template<class...> class Op, class... Args>
        constexpr bool is_detected_v = is_detected<Op, Args...>::value;

        template<template<class...> class Op, class... Args>
        using detected_t= typename detector<nonesuch, void, Op, Args...>::type;


// detected_or
        template<class Default, template<class...> class Op, class... Args>
        using detected_or= detector<Default, void, Op, Args...>;

        template<class Default, template<class...> class Op, class... Args>
        using detected_or_t= typename detected_or<Default, Op, Args...>::type;

// is_detected_exact
        template<class Expected, template<class...> class Op, class... Args>
        using is_detected_exact= std::is_same<Expected, detected_t<Op, Args...> >;

        template<class Expected, template<class...> class Op, class... Args>
        constexpr bool is_detected_exact_v = is_detected_exact<Expected, Op, Args...>::value;

// is_detected_convertible
        template<class To, template<class...> class Op, class... Args>
        using is_detected_convertible = std::is_convertible<detected_t<Op, Args...>, To>;

        template<class To, template<class...> class Op, class... Args>
        constexpr bool is_detected_convertible_v = is_detected_convertible<To, Op, Args...>::value;
} // namespace caribou::internal
#endif //CARIBOU_TRAITS_H
