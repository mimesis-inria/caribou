#ifndef CARIBOU_CONFIG_H
#define CARIBOU_CONFIG_H

#include <cstdint>

#define CARIBOU_USE_DOUBLE

#ifdef CARIBOU_USE_DOUBLE
#define FLOATING_POINT_TYPE double
#define INTEGER_TYPE signed long long int
#define UNSIGNED_INTEGER_TYPE unsigned long long int
#define EPSILON  ((double)1.0e-19)
#else
#define FLOATING_POINT_TYPE float
#define INTEGER_TYPE int32_t
#define UNSIGNED_INTEGER_TYPE uint32_t
#define EPSILON  ((float) 1.0e-15)
#endif

#ifdef __cpp_if_constexpr
#define CONSTEXPR_IF constexpr
#else
#define CONSTEXPR_IF
#endif

#endif //CARIBOU_CONFIG_H
