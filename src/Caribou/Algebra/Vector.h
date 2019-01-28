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

    /**
     * Constructor by c-array or initializer list.
     */
    constexpr
    Matrix(ValueType const (&components)[R_]) noexcept
    : Base (components)
    {}

    /**
     * Copy constructor from another column-vector of the same data type
     */
    constexpr
    Matrix(const Matrix<R_, 1, ValueType> & other) noexcept
    : Base(other)
    {}

    /** Constructor from a list of parameters (each parameter is a scalar component of the matrix) **/
    template<
            typename ...Args,
            typename std::enable_if<R_ == sizeof...(Args) + 1, int>::type = 0,
            typename std::enable_if<std::is_integral<ValueType>::value or std::is_floating_point<ValueType>::value ,int>::type = 0>
    constexpr
    Matrix(ValueType first_value, Args&&...e) noexcept
    : Base(first_value, std::forward<Args>(e)...)
    {}
};

/**
 * Generic column vector (1xC matrix).
 */
template <size_t C_, typename ValueType>
struct Matrix<1,C_, ValueType> : public internal::BaseVector<Matrix, 1, C_, ValueType>
{
    using Base = internal::BaseVector<Matrix, 1, C_, ValueType>;
    using Base::Base;

    /**
     * Constructor by c-array or initializer list.
     */
    constexpr
    Matrix(ValueType const (&components)[C_]) noexcept
    : Base (components)
    {}

    /**
     * Copy constructor from another row-vector of the same data type
     */
    constexpr
    Matrix(const Matrix<1, C_, ValueType> & other) noexcept
    : Base (other)
    {}

    /** Constructor from a list of parameters (each parameter is a scalar component of the matrix) **/
    template<
            typename ...Args,
            typename std::enable_if<C_ == sizeof...(Args) + 1, int>::type = 0,
            typename std::enable_if<std::is_integral<ValueType>::value or std::is_floating_point<ValueType>::value ,int>::type = 0>
    constexpr
    Matrix(ValueType first_value, Args&&...e) noexcept
    : Base(first_value, std::forward<Args>(e)...)
    {}

    /** Alias for inner_product (dot product). */
    using Base::operator*;
    template <typename OtherValueType>
    constexpr ValueType
    operator*(const Matrix<C_, 1, OtherValueType> & other) const noexcept
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

    /**
     * Copy constructor from another column-vector of the same data type
     */
    constexpr
    Matrix(const Matrix<3, 1, ValueType> & other) noexcept
    : Base(other)
    {}

    /** Constructor from a list of parameters (each parameter is a scalar component of the matrix) **/
    template<
            typename ...Args,
            typename std::enable_if<3 == sizeof...(Args) + 1, int>::type = 0,
            typename std::enable_if<std::is_integral<ValueType>::value or std::is_floating_point<ValueType>::value ,int>::type = 0>
    constexpr
    Matrix(ValueType first_value, Args&&...e) noexcept
    : Base(first_value, std::forward<Args>(e)...)
    {}
};

/**
 * Column vector (1x3 matrix).
 */
template <typename ValueType>
struct Matrix<1,3, ValueType> : public internal::BaseVector3D<Matrix, 1, 3, ValueType>
{
    using Base = internal::BaseVector3D<Matrix, 1, 3, ValueType>;
    using Base::Base;

    /**
     * Copy constructor from another row-vector of the same data type
     */
    constexpr
    Matrix(const Matrix<1, 3, ValueType> & other) noexcept
    : Base(other)
    {}

    /** Constructor from a list of parameters (each parameter is a scalar component of the matrix) **/
    template<
            typename ...Args,
            typename std::enable_if<3 == sizeof...(Args) + 1, int>::type = 0,
            typename std::enable_if<std::is_integral<ValueType>::value or std::is_floating_point<ValueType>::value ,int>::type = 0>
    constexpr
    Matrix(ValueType first_value, Args&&...e) noexcept
    : Base(first_value, std::forward<Args>(e)...)
    {}

    /** Alias for inner_product (dot product). */
    using Base::operator*;
    template <typename OtherValueType>
    constexpr ValueType
    operator*(const Matrix<3, 1, OtherValueType> & other) const noexcept
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

    /**
     * Copy constructor from another row-vector of the same data type
     */
    constexpr
    Matrix(const Matrix<1, 1, ValueType> & other) noexcept
    : Base(other)
    {}


    constexpr
    Matrix (const ValueType & value) noexcept {
        (*this)[0] = value;
    }

    operator ValueType() const noexcept
    { return (*this)[0]; }
};

/**
 * Vector Alias for Rx1 matrix.
 */
template <size_t Dimension, typename ValueType=FLOATING_POINT_TYPE>
using Vector = Matrix<Dimension, 1, ValueType>;

/**
 * Vector Alias for Rx1 column vector.
 */
template <size_t Dimension, typename ValueType=FLOATING_POINT_TYPE>
using ColumnVector = Matrix<Dimension, 1, ValueType>;

/**
 * Vector Alias for 1xC row vector.
 */
template <size_t Dimension, typename ValueType=FLOATING_POINT_TYPE>
using RowVector = Matrix<1, Dimension, ValueType>;

} // namespace algebra
} // namespace caribou

#endif //CARIBOU_ALGEBRA_VECTOR_H
