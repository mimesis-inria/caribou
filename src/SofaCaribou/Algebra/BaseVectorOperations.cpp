#include <SofaCaribou/config.h>
#include <SofaCaribou/Algebra/BaseVectorOperations.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/version.h>
#if (defined(SOFA_VERSION) && SOFA_VERSION < 211299)
#include <SofaBaseLinearSolver/FullVector.h>
#else
#include <sofa/linearalgebra/FullVector.h>
#endif // #if (defined(SOFA_VERSION) && SOFA_VERSION < 211299)
DISABLE_ALL_WARNINGS_END

#if (defined(SOFA_VERSION) && SOFA_VERSION < 201200)
namespace sofa {
using Size = sofa::defaulttype::BaseVector::Index;
using Index = sofa::defaulttype::BaseVector::Index;
}
#endif

namespace SofaCaribou::Algebra {

namespace { // Anonymous
#if (defined(SOFA_VERSION) && SOFA_VERSION < 211299)
using namespace sofa::component::linearsolver;
#else
using namespace sofa::linearalgebra;
#endif

template <typename Real>
auto to_full_vector(const BaseVector* v) -> const FullVector<Real> * {
    return dynamic_cast<const FullVector<Real> *>(v);
}

// Dot product between a SOFA full vector and any BaseVector
template <typename ReturnValue, typename Real>
ReturnValue dot(const FullVector<Real> * v1, const BaseVector * v2) {
    const auto n = static_cast<sofa::Size>(v1->size());
    auto value = static_cast<ReturnValue>(0);
    for (sofa::Index i = 0; i < n; ++i) {
        value +=  ( v1->element(i) * v2->element(i) );
    }

    return value;
}

// Dot product between a two SOFA full vectors
template <typename ReturnValue, typename Real1, typename Real2>
ReturnValue dot(const FullVector<Real1> * v1, const FullVector<Real2> * v2) {
    const auto n = static_cast<sofa::Size>(v1->size());
    auto value = static_cast<ReturnValue>(0);
    for (sofa::Index i = 0; i < n; ++i) {
        value +=  ( v1->element(i) * v2->element(i) );
    }

    return value;
}
}

/** Compute the dot product between two BaseVector, i.e. scalar = v1.dot(v2) */
 double dot(const BaseVector * v1, const BaseVector * v2) {
    caribou_assert(v1->size() == v2->size());

    // todo(jnbrunet2000@gmail.com): This function is awful. I'm sure there is a much much cleaner and faster way
    //                               to do this. But I'm running out of time.

    // Most common case first FullVector<SReal> - FullVector<SReal> :
    {
        const auto * full_v1 = to_full_vector<SReal >(v1);
        const auto * full_v2 = to_full_vector<SReal >(v2);
        if (full_v1 && full_v2) {
            return full_v1->dot(*full_v2);
        }
    }

    // Next, we need to find the good virtual type....

    // FullVector<double> - *
    if (const auto * full_v1 = to_full_vector<double>(v1)) {
        if (const auto * full_v2 = to_full_vector<double>(v2)) {
            // FullVector<double> - FullVector<double>
            return full_v1->dot(*full_v2);
        } else if (const auto * full_v2 = to_full_vector<float>(v2)) {
            // FullVector<double> - FullVector<float>
            return dot<double>(full_v1, full_v2);
        } else {
            // FullVector<double> - BaseVector
            return dot<double>(full_v1, v2);
        }
    }

    // FullVector<float> - *
    if (const auto * full_v1 = to_full_vector<float>(v1)) {
        if (const auto * full_v2 = to_full_vector<float>(v2)) {
            // FullVector<float> - FullVector<float>
            return full_v1->dot(*full_v2);
        } else if (const auto * full_v2 = to_full_vector<double>(v2)) {
            // FullVector<float> - FullVector<double>
            return dot<double>(full_v1, full_v2);
        } else {
            // FullVector<float> - BaseVector
            return dot<double>(full_v1, v2);
        }
    }

    // v1 isn't a full vector... Maybe v2 is?
    // * - FullVector<double>
    if (const auto * full_v2 = to_full_vector<double>(v2)) {
        return dot<double>(full_v2, v1);
    }

    // * - FullVector<float>
    if (const auto * full_v2 = to_full_vector<float>(v2)) {
        return dot<double>(full_v2, v1);
    }

    // Generic case (unoptimized!)
    // * - *
    const auto n = static_cast<sofa::Size>(v1->size());
    auto value = static_cast<SReal>(0);
    for (sofa::Index i = 0; i < n; ++i) {
        value +=  ( v1->element(i) * v2->element(i) );
    }

    return static_cast<double> (value);
}

} // namespace SofaCaribou::Algebra
