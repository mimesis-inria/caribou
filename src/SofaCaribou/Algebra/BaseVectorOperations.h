#pragma once

#include <SofaCaribou/config.h>

// Various utilities to perform numerical operations on SOFA's BaseVector.
// These utilities are responsible to automatically find the type of vector
// and perform the optimal operations on them.
DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h> // for SOFA_VERSION
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
namespace sofa::defaulttype {
    class BaseVector;
}
#else
namespace sofa::linearalgebra {
    class BaseVector;
}
#endif // (defined(SOFA_VERSION) && SOFA_VERSION < 211299)

namespace SofaCaribou::Algebra {

#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
using BaseVector = sofa::defaulttype::BaseVector;
#else
using BaseVector = sofa::linearalgebra::BaseVector;
#endif

/**
 * Compute the dot product between two BaseVector, i.e. scalar = v1.dot(v2) *
 *
 * @todo (jnbrunet2000@gmail.com) : this could be much more optimal using delegate functions. For example, we could first
 *                                  get a callback to an optimal dot function for these types of v1 and v2, and later
 *                                  simply reuse this callback. Hence we avoid the need to do multiple virtual checks
 *                                  each time this function is called
 */
 double dot(const BaseVector * v1, const BaseVector * v2);

} // namespace SofaCaribou::Algebra
