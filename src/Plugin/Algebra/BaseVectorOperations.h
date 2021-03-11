#pragma once

// Various utilities to perform numerical operations on SOFA's BaseVector.
// These utilities are responsible to automatically find the type of vector
// and perform the optimal operations on them.

namespace sofa::defaulttype {
class BaseVector;
}

namespace SofaCaribou::Algebra {


/**
 * Compute the dot product between two BaseVector, i.e. scalar = v1.dot(v2) *
 *
 * @todo (jnbrunet2000@gmail.com) : this could be much more optimal using delegate functions. For example, we could first
 *                                  get a callback to an optimal dot function for these types of v1 and v2, and later
 *                                  simply reuse this callback. Hence we avoid the need to do multiple virtual checks
 *                                  each time this function is called
 */
double dot(const sofa::defaulttype::BaseVector * v1, const sofa::defaulttype::BaseVector * v2);

} // namespace SofaCaribou::Algebra