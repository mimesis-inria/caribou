#ifndef CARIBOU_CONFIG_H
#define CARIBOU_CONFIG_H

#include <cstdint>

#ifdef CARIBOU_USE_DOUBLE
#define FLOATING_POINT_TYPE double
#define INTEGER_TYPE int64_t
#define UNSIGNED_INTEGER_TYPE uint64_t
#define EPSILON  ((float)1.0e-100)
#else
#define FLOATING_POINT_TYPE float
#define INTEGER_TYPE int32_t
#define UNSIGNED_INTEGER_TYPE uint32_t
#define EPSILON  ((float) 1.0e-20)
#endif

#ifdef __cpp_if_constexpr
#define CONSTEXPR_IF constexpr
#else
#define CONSTEXPR_IF
#endif

#endif //CARIBOU_CONFIG_H
