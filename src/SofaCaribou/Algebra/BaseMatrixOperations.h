#pragma once

#include <SofaCaribou/config.h>

// Various utilities to perform numerical operations on SOFA's BaseMatrix.
DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h> // for SOFA_VERSION
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
#include <sofa/defaulttype/BaseMatrix.h>
#else
#include <sofa/linearalgebra/BaseMatrix.h>
#endif // (defined(SOFA_VERSION) && SOFA_VERSION < 211299)

namespace SofaCaribou::Algebra
{

#if (defined(SOFA_VERSION) && SOFA_VERSION < 211200)
using BaseMatrix = sofa::defaulttype::BaseMatrix;
#else
using BaseMatrix = sofa::linearalgebra::BaseMatrix;
#endif // (defined(SOFA_VERSION) && SOFA_VERSION < 211299)


} // namespace SofaCaribou::Algebra
