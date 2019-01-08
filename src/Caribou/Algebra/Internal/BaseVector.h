#ifndef CARIBOU_ALGEBRA_INTERNAL_BASEVECTOR_H
#define CARIBOU_ALGEBRA_INTERNAL_BASEVECTOR_H

#include <numeric>
#include <cmath>

#include <Caribou/config.h>
#include <Caribou/Algebra/Internal/BaseMatrix.h>

namespace caribou {
namespace algebra {
namespace internal {

/**
 * Vector Rx1 matrix.
 *
 * ** Do not use this class directly. Use instead caribou::algebra::Vector<R> which is an alias to matrix<R,1>. **
 *
 * This class extend the BaseMatrix by adding functions only available on vectors. It is meant to be derived by
 * a partial specialization of the class Matrix<R, C> (see later in this file for the row-vector, column-vector
 * and 3D implementations).
 *
 * The functions declared in this class ca be used with any type of vectors (1D, 2D, 3D, ...).
 *
 * @tparam MatrixType_ <R_, C_, ValueType_> See BaseMatrix template MatrixType_.
 * @tparam R_ The dimension of the vector (number of rows in a Rx1 vector-matrix).
 * @tparam C_ The dimension of the vector (number of columns in a 1xC vector-matrix).
 * @tparam ValueType_ The data type of the vector's components (default to float)/
 */
template <template <size_t, size_t, typename> class MatrixType_, size_t R_, size_t C_, typename ValueType_>
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

    //////////////////////
    //// Constructors ////
    //////////////////////

    using Base::Base; // Importing BaseMatrix constructors

    ///////////////////
    //// Operators ////
    ///////////////////

    /** Comparison operator with a vector of a same dimension and with component data type of OtherValueType **/
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == N and OtherC == 1) or (OtherC == N and OtherR == 1), int>::type = 0
    >
    inline constexpr bool
    operator==(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const noexcept
    {
        return std::equal(this->begin(), this->end(), other.begin(), [](const ValueType & v1, const OtherValueType & v2) -> bool {
            return (-EPSILON <= (v1 - v2) and (v1 - v2) <= +EPSILON);
        });
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
    inner_product(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const noexcept
    {
        return std::inner_product(std::begin(other), std::end(other), std::begin(*this), (ValueType) 0);
    }

    /** Alias for inner_product. */
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == N and OtherC == 1) or (OtherC == N and OtherR == 1), int>::type = 0
    >
    inline constexpr ValueType
    dot(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const noexcept
    {
        return inner_product(other);
    }

    /** Alias for inner_product. */
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == N and OtherC == 1) or (OtherC == N and OtherR == 1), int>::type = 0
    >
    inline constexpr ValueType
    scalar_product(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const noexcept
    {
        return inner_product(other);
    }

    /////////////////////////////////
    //// Mathematical properties ////
    /////////////////////////////////

    /** Compute the length of the vector. **/
    inline constexpr ValueType
    length() const noexcept
    {
        return sqrt(length_squared());
    };

    /** Compute the length^2 of the vector. **/
    inline constexpr ValueType
    length_squared() const noexcept
    {
        using Vec = const MatrixType<R,C,ValueType>&;
        const auto & self = static_cast<Vec>(*this);

        return self.T() * self;
    }

    /** Get the unit vector (current vector normalized to unit length) **/
    inline constexpr MatrixType<R, C>
    unit() const noexcept
    {
        using Vec = const MatrixType<R,C,ValueType>&;
        const auto & self = static_cast<Vec>(*this);

        return self/length();
    }

    /////////////////////////////////
    ////      Tool functions     ////
    /////////////////////////////////
    inline std::string
    to_string() const noexcept
    {
        return std::string("(") +
               std::accumulate(std::next(this->begin()), this->end(), std::to_string(this->at(0)), [](const std::string & s, const ValueType & component) -> std::string {
                   return s + std::string(", ") + std::to_string(component);
               }) +
               std::string(")");
    }

private:
    template<size_t current_index, typename... Args>
    void recursive_set(ValueType component, Args&&... other_components) noexcept {
        (*this)[current_index] = component;
        recursive_set<current_index+1>(std::forward<Args>(other_components)...);
    }

    template<size_t current_index>
    void recursive_set(ValueType component) noexcept {
        (*this)[current_index] = component;
    }
};

/** 3-dimensional vectors */
template <template <size_t, size_t, typename> class MatrixType_, size_t R_, size_t C_, typename ValueType_>
struct BaseVector3D : public internal::BaseVector<MatrixType_, R_, C_, ValueType_>
{
    static constexpr size_t R = R_; ///< Number of rows in a Rx1 vector-matrix
    static constexpr size_t C = C_; ///< Number of columns in a 1xC vector-matrix
    static constexpr size_t N = R * C; ///< Number of elements

    static_assert((R == 3 and C == 1) or (R == 1 and C == 3),
                  "A 3D vector must be of dimension 3x1 (row vector) or 1x3 (column vector).");

    template<size_t R__ = R, size_t C__ = C, typename ValueType__= ValueType_>
    using MatrixType = MatrixType_<R__, C__, ValueType__>;
    using ValueType = ValueType_;
    using Base = internal::BaseVector<MatrixType_, R_, C_, ValueType_>;
    using Base::Base;

    /** Compute the cross product with another vector. Only supported for vectors of dimension 3. */
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == 3 and OtherC == 1) or (OtherR == 1 and OtherC == 3), int>::type = 0
    >
    inline constexpr MatrixType<3, 1>
    cross_product(const OtherVectorType<OtherR, OtherC, OtherValueType> & b) const noexcept
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
    cross(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const noexcept
    {
        return cross_product(other);
    }

    /** Alias for cross_product */
    template <
            template <size_t, size_t, typename> class OtherVectorType, size_t OtherR, size_t OtherC, typename OtherValueType,
            typename std::enable_if<(OtherR == 3 and OtherC == 1) or (OtherR == 1 and OtherC == 3), int>::type = 0
    >
    inline constexpr MatrixType<3, 1>
    operator^(const OtherVectorType<OtherR, OtherC, OtherValueType> & other) const noexcept
    {
        return cross_product(other);
    }

};

} // namespace internal
} // namespace algebra
} // namespace caribou

#endif //CARIBOU_ALGEBRA_INTERNAL_BASEVECTOR_H