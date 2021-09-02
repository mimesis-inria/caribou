#ifndef CARIBOU_CONFIG_H
#define CARIBOU_CONFIG_H

#include <cstdint>
#include <iso646.h>

/* #undef CARIBOU_USE_FLOAT */

#ifndef NDEBUG
# ifndef CARIBOU_DEBUG
#  define CARIBOU_DEBUG
# endif
# include <cassert>
#endif


#ifdef CARIBOU_USE_FLOAT
#define FLOATING_POINT_TYPE float
#define EPSILON  ((float) 1.0e-15)
#else
#define FLOATING_POINT_TYPE double
#define EPSILON  ((double)1.0e-19)
#endif

#define INTEGER_TYPE signed long long int
#define UNSIGNED_INTEGER_TYPE unsigned long long int

#ifdef CARIBOU_DEBUG
    #define caribou_assert(x) assert(x)
#else
    #define caribou_assert(x)
#endif

#ifdef __cpp_if_constexpr
#define CONSTEXPR_IF constexpr
#else
#define CONSTEXPR_IF
#endif

#endif //CARIBOU_CONFIG_H
