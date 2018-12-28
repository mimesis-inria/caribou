#ifndef CARIBOU_ALGEBRA_VECTOR_H
#define CARIBOU_ALGEBRA_VECTOR_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Matrix.h>
#include <Caribou/Algebra/Internal/BaseVector.h>

namespace caribou {
namespace algebra {

/**
 * Generic row vector (Rx1 matrix).
 */
template <size_t R_, typename ValueType>
struct Matrix<R_,1, ValueType> : public internal::BaseVector<Matrix, R_, 1, ValueType>
{
    using Base = internal::BaseVector<Matrix, R_, 1, ValueType>;
    using Base::Base;
};

/**
 * Generic column vector (1xC matrix).
 */
template <size_t C_, typename ValueType>
struct Matrix<1,C_, ValueType> : public internal::BaseVector<Matrix, 1, C_, ValueType>
{
    using Base = internal::BaseVector<Matrix, 1, C_, ValueType>;
    using Base::Base;

    /** Alias for inner_product (dot product). */
    using Base::operator*;
    template <typename OtherValueType>
    constexpr ValueType
    operator*(const Matrix<C_, 1, OtherValueType> & other) const
    {
        return this->inner_product(other);
    }
};

/**
 * Row vector (3x1 matrix).
 */
template <typename ValueType>
struct Matrix<3,1, ValueType> : public internal::BaseVector3D<Matrix, 3, 1, ValueType>
{
    using Base = internal::BaseVector3D<Matrix, 3, 1, ValueType>;
    using Base::Base;
};

/**
 * Column vector (1x3 matrix).
 */
template <typename ValueType>
struct Matrix<1,3, ValueType> : public internal::BaseVector3D<Matrix, 1, 3, ValueType>
{
    using Base = internal::BaseVector3D<Matrix, 1, 3, ValueType>;
    using Base::Base;

    /** Alias for inner_product (dot product). */
    using Base::operator*;
    template <typename OtherValueType>
    constexpr ValueType
    operator*(const Matrix<3, 1, OtherValueType> & other) const
    {
        return this->inner_product(other);
    }
};

/**
 * 1D vector (1x1 matrix).
 */
template <typename ValueType>
struct Matrix<1,1, ValueType> : public internal::BaseVector<Matrix, 1, 1, ValueType>
{
    using Base = internal::BaseVector<Matrix, 1, 1, ValueType>;
    using Index = typename Base::Index;

    using Base::Base;
};

/**
 * Vector Alias for Rx1 matrix.
 */
template <size_t Dimension, typename ValueType = FLOATING_POINT_TYPE>
using Vector = Matrix<Dimension, 1, ValueType>;

} // namespace algebra
} // namespace caribou

#endif //CARIBOU_ALGEBRA_VECTOR_H
