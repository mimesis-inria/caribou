#ifndef CARIBOU_ALGEBRA_VECTOR_H
#define CARIBOU_ALGEBRA_VECTOR_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Matrix.h>

namespace caribou {
namespace algebra {

/**
 * Vector Rx1 matrix.
 *
 * ** Do not use this class directly. Use instead caribou::algebra::Vector<R> which is an alias to matrix<R,1>. **
 *
 * This class extend the BaseMatrix by adding functions only available on vectors. It is meant to be derived by
 * a partial specialization of the class Matrix<R, C> (see later in this file for the row-vector, column-vector
 * and 3D specializations).
 *
 * The functions declared in this class ca be used with any type of vectors (1D, 2D, 3D, ...).
 *
 * @tparam MatrixType_ <R_, C_, ValueType_> See BaseMatrix template MatrixType_.
 * @tparam R_ The dimension of the vector (number of rows in a Rx1 vector-matrix).
 * @tparam C_ The dimension of the vector (number of columns in a 1xC vector-matrix).
 * @tparam ValueType_ The data type of the vector's components (default to float)/
 */
template <template <size_t, size_t, typename> class MatrixType_, size_t R_, size_t C_, typename ValueType_=FLOATING_POINT_TYPE>
struct BaseVector : public BaseMatrix<MatrixType_, R_, C_, ValueType_>
{
    static constexpr size_t R = R_; ///< Number of rows in a Rx1 vector-matrix
    static constexpr size_t C = C_; ///< Number of columns in a 1xC vector-matrix
    static constexpr size_t N = R*C; ///< Number of elements

    static_assert((R >= 1 and C == 1) or (R == 1 and C >= 1), "A vector must be of dimension Rx1 (row vector) or 1xC (column vector).");

    template<size_t R__ = R, size_t C__ = C, typename ValueType__= ValueType_>
    using MatrixType = MatrixType_<R__, C__, ValueType__>;

    using ValueType = ValueType_;

    using Base = BaseMatrix<MatrixType_, R_, C_, ValueType_>;
    using Index = typename Base::Index;

    using Base::Base;

    //////////////////////
    //// Constructors ////
    //////////////////////


    /**
     * Constructor by c-array or initializer list
     */
    template <typename OtherValueType>
    BaseVector(OtherValueType const (&components)[N]) {
        for (size_t i = 0; i < N; ++i) {
            (*this)[i] = static_cast<ValueType>(components[i]);
        }
    }

    /**
     * Constructor by a list initializer (ex. Vector<3> A {0, 100, 0}).
     */
    template <typename OtherValueType>
    BaseVector(std::initializer_list<OtherValueType> il) {
        std::copy(il.begin(), il.end(), this->begin());
    }

    /**
     * Constructor by another vector type. The vector type must implements the [] operator to access its components.
     * It must also have a template signature of vector<N, ValueType> where N is its dimension (number of components) and
     * ValueType is the type of its components.
     */
    template <template <int, typename> class OtherVectorType, typename OtherValueType>
    BaseVector(const OtherVectorType<N, OtherValueType> & other) {
        for (size_t i = 0; i < N; ++i)
            (*this)[i] = static_cast<ValueType> (other[i]);
    }

    /**
     * Constructor by variadic arguments.
     *
     * Example:
     * \code{.cpp}
     * Vector<3> v (0,1,0);
     * \endcode
     *
     * Warning: The data type of the components will the one declared on the class template
     * specifier (float if nothing is specified), hence it will not be the one of the arguments passed.
     *
     * For example, the following vector:
     *
     * \code{.cpp}
     * Vector<3, char> v ((float) 0.2, (float) 1.1, (float) 1.4);
     * \endcode
     *
     * will have its coordinates as char values, and not as float values.
     */
    template<typename... Args>
    BaseVector(ValueType first_component, Args&&... other_components)
    {
        recursive_set<0> (first_component, std::forward<Args>(other_components)...);
    }


    /////////////////////////////////
    //// Mathematical operations ////
    /////////////////////////////////

    /** Compute the inner product with another vector (aliases are dot, scalar_product and operator*). */
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == N and OtherC == 1) or (OtherC == N and OtherR == 1), int>::type = 0
            >
    inline constexpr ValueType
    inner_product(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const
    {
        return std::inner_product(std::begin(other), std::end(other), std::begin(*this), (ValueType) 0);
    }

    /** Alias for inner_product. */
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == N and OtherC == 1) or (OtherC == N and OtherR == 1), int>::type = 0
    >
    inline constexpr ValueType
    dot(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const
    {
        return inner_product(other);
    }

    /** Alias for inner_product. */
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == N and OtherC == 1) or (OtherC == N and OtherR == 1), int>::type = 0
    >
    inline constexpr ValueType
    scalar_product(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const
    {
        return inner_product(other);
    }

    /////////////////////////////////
    //// Mathematical properties ////
    /////////////////////////////////

    /** Compute the length of the vector. **/
    inline constexpr ValueType
    length() const
    {
        return sqrt(length_squared());
    };

    /** Compute the length^2 of the vector. **/
    inline constexpr ValueType
    length_squared() const
    {
        using Vec = const MatrixType<R,C,ValueType>&;
        const auto & self = static_cast<Vec>(*this);

        return self.T() * self;
    }

    /** Get the unit vector (current vector normalized to unit length) **/
    inline constexpr MatrixType<R, C>
    unit() const
    {
        using Vec = const MatrixType<R,C,ValueType>&;
        const auto & self = static_cast<Vec>(*this);

        return self/length();
    }

    /////////////////////////////////
    ////      Tool functions     ////
    /////////////////////////////////
    inline std::string
    to_string() const
    {
        return std::string("(") +
               std::accumulate(std::next(this->begin()), this->end(), std::to_string(this->at(0)), [](const std::string & s, const ValueType & component) -> std::string {
                   return s + std::string(", ") + std::to_string(component);
               }) +
               std::string(")");
    }

private:
    template<size_t current_index, typename... Args>
    void recursive_set(ValueType component, Args&&... other_components) {
        (*this)[current_index] = component;
        recursive_set<current_index+1>(std::forward<Args>(other_components)...);
    }

    template<size_t current_index>
    void recursive_set(ValueType component) {
        (*this)[current_index] = component;
    }

};

/** 3-dimensional vectors */
template <template <size_t, size_t, typename> class MatrixType_, size_t R_, size_t C_, typename ValueType_=FLOATING_POINT_TYPE>
struct BaseVector3D : public BaseVector<MatrixType_, R_, C_, ValueType_>
{
    static constexpr size_t R = R_; ///< Number of rows in a Rx1 vector-matrix
    static constexpr size_t C = C_; ///< Number of columns in a 1xC vector-matrix
    static constexpr size_t N = R * C; ///< Number of elements

    static_assert((R == 3 and C == 1) or (R == 1 and C == 3),
                  "A 3D vector must be of dimension 3x1 (row vector) or 1x3 (column vector).");

    template<size_t R__ = R, size_t C__ = C, typename ValueType__= ValueType_>
    using MatrixType = MatrixType_<R__, C__, ValueType__>;
    using ValueType = ValueType_;
    using Base = BaseVector<MatrixType_, R_, C_, ValueType_>;
    using Base::Base;

    /** Compute the cross product with another vector. Only supported for vectors of dimension 3. */
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == 3 and OtherC == 1) or (OtherR == 1 and OtherC == 3), int>::type = 0
    >
    inline constexpr MatrixType<3, 1>
    cross_product(const OtherVectorType<OtherR, OtherC, OtherValueType> & b) const
    {
        return {
                (*this)[1]*b[2] -(*this)[2]*b[1],
                (*this)[2]*b[0] -(*this)[0]*b[2],
                (*this)[0]*b[1] -(*this)[1]*b[0],
        };
    }

    /** Alias for cross_product **/
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == 3 and OtherC == 1) or (OtherR == 1 and OtherC == 3), int>::type = 0
    >
    inline constexpr MatrixType<3, 1>
    cross(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const
    {
        return cross_product(other);
    }

    /** Alias for cross_product */
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == 3 and OtherC == 1) or (OtherR == 1 and OtherC == 3), int>::type = 0
    >
    inline constexpr MatrixType<3, 1>
    operator^(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const
    {
        return cross_product(other);
    }

};

/**
 * Row vector (Rx1 matrix).
 */
template <size_t R_, typename ValueType>
struct Matrix<R_,1, ValueType> : public BaseVector<Matrix, R_, 1, ValueType>
{
    using Base = BaseVector<Matrix, R_, 1, ValueType>;
    using Base::Base;
};

/**
 * Column vector (1xC matrix).
 */
template <size_t C_, typename ValueType>
struct Matrix<1,C_, ValueType> : public BaseVector<Matrix, 1, C_, ValueType>
{
    using Base = BaseVector<Matrix, 1, C_, ValueType>;
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
struct Matrix<3,1, ValueType> : public BaseVector3D<Matrix, 3, 1, ValueType>
{
    using Base = BaseVector3D<Matrix, 3, 1, ValueType>;
    using Base::Base;
};

/**
 * Column vector (1x3 matrix).
 */
template <typename ValueType>
struct Matrix<1,3, ValueType> : public BaseVector3D<Matrix, 1, 3, ValueType>
{
    using Base = BaseVector3D<Matrix, 1, 3, ValueType>;
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
struct Matrix<1,1, ValueType> : public BaseVector<Matrix, 1, 1, ValueType>
{
    using Base = BaseVector<Matrix, 1, 1, ValueType>;
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
