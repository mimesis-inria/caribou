#ifndef CARIBOU_CONFIG_H
#define CARIBOU_CONFIG_H

#ifdef CARIBOU_USE_DOUBLE
#define FLOATING_POINT_TYPE double
#define EPSILON  ((float)1.0e-100)
#else
#define FLOATING_POINT_TYPE float
#define EPSILON  ((float) 1.0e-20)
#endif

#endif //CARIBOU_CONFIG_H
