#ifndef CARIBOU_ALGEBRA_MATRIX_H
#define CARIBOU_ALGEBRA_MATRIX_H

#include <ostream>
#include <cstddef>
#include <array>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <cmath>

#include <Caribou/config.h>

namespace caribou {
namespace algebra {

/**
 * A simple representation of a matrix.
 *
 * ** Do not use this class directly. Use instead caribou::algebra::Matrix. **
 *
 * The functions declared in this class can be used with any type of matrices (squared, vectors, rectangular).
 *
 * Do to so, it uses the Curiously recurring template pattern (CRTP) :
 *    https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
 *
 * @tparam MatrixType_ <R_, C_, ValueType_> The class type that will specialize or add more functions to this class.
 * @tparam R_ The number of rows
 * @tparam C_ The number of columns
 * @tparam ValueType_ The data type of the matrix's components (default to float)
 */
template <template <size_t, size_t, typename> class MatrixType_, size_t R_, size_t C_, typename ValueType_=FLOATING_POINT_TYPE>
struct BaseMatrix : public std::array<ValueType_, R_*C_>
{
    static constexpr size_t R = R_; ///< Number of rows
    static constexpr size_t C = C_; ///< Number of columns per row
    static constexpr size_t N = R*C; ///< Number of elements


    template<size_t R__ = R, size_t C__= C, typename ValueType__= ValueType_>
    using MatrixType = MatrixType_<R__, C__, ValueType__>;

    using ValueType = ValueType_;
    using Self = MatrixType<R, C, ValueType>;
    using Row = MatrixType<1, C, ValueType>;
    using Column = MatrixType<R, 1, ValueType>;
    using Index = size_t;

    //////////////////////
    //// Constructors ////
    //////////////////////

    /**
     * Main constructor
     * @param initialize_to_zero If true, initialize the scalar components of the vector to zero
     */
    explicit
    BaseMatrix(bool initialize_to_zero = false) {
        if (initialize_to_zero) {
            this->fill(0);
        }
    }

    /**
     * Constructor by c-array or initializer list.
     * Ex:
     * Matrix<3,3> A (
     *   {
     *     {1,2,3}, // Row 0
     *     {4,5,6}, // Row 1
     *     {7,8,9} //  Row 2
     *   }
     * );.
     */
    BaseMatrix(ValueType const (&components)[R][C]) {
        for (size_t row = 0; row < R; ++row) {
            for (size_t column = 0; column < C; ++column) {
                (*this)[row*C + column] = components[row][column];
            }
        }
    }

    /**
     * Constructor by c-array or initializer list.
     * Ex:
     * Matrix<3,3> A (
     *   {
     *     {1,2,3}, // Row 0
     *     {4,5,6}, // Row 1
     *     {7,8,9} //  Row 2
     *   }
     * );.
     */
    template <typename OtherValueType>
    BaseMatrix(OtherValueType const (&components)[R][C]) {
        for (size_t row = 0; row < R; ++row) {
            for (size_t column = 0; column < C; ++column) {
                (*this)[row*C + column] = static_cast<ValueType>(components[row][column]);
            }
        }
    }

    /**
     * Copy constructor from another matrix of a different data type
     */
    template <template <size_t, size_t, typename> class OtherMatrixType, typename OtherValueType>
    BaseMatrix(const OtherMatrixType<R, C, OtherValueType> & other) {
        for (size_t row = 0; row < R; ++row) {
            for (size_t column = 0; column < C; ++column) {
                (*this)[row*C + column] = static_cast<ValueType>(other(row, column));
            }
        }
    }

    ///////////////////
    //// Accessors ////
    ///////////////////

    inline constexpr ValueType &
    operator () (const Index & row, const Index & column)
    {
        return (*this)[row*C + column];
    }

    inline constexpr const ValueType &
    operator () (const Index & row, const Index & column) const
    {
        return (*this)[row*C + column];
    }

    ///////////////////
    //// Operators ////
    ///////////////////

    /** Assignment operator **/
    template<typename OtherValueType>
    inline Self &
    operator = (const MatrixType<R, C, OtherValueType> & other)
    {
        for (size_t i = 0; i < N; ++i)
            (*this)[i] = other[i];
        return static_cast<Self&>(*this);
    }


    /** Comparison operator with a matrix of a same dimension and with component data type of OtherValueType **/
    template<typename OtherValueType>
    inline constexpr bool
    operator==(const MatrixType<C, R, OtherValueType> & other) const
    {
        return std::equal(this->begin(), this->end(), other.begin(), [](const ValueType & v1, const OtherValueType & v2) -> bool {
            return (-EPSILON <= (v1 - v2) and (v1 - v2) <= +EPSILON);
        });
    }

    /** Matrix multiplication **/
    template<size_t M, typename OtherValueType>
    inline MatrixType<R, M, ValueType>
    operator*(const MatrixType<C, M, OtherValueType> & other) const
    {
        MatrixType<R, M, ValueType> r;
        for (size_t i = 0; i < R; ++i) {
            for (size_t j = 0; j < M; ++j) {
                r(i, j)=(*this)(i,0) * other(0, j);
                for(size_t k=1; k<C; k++)
                    r(i, j) += (*this)(i, k) * other(k, j);
            }
        }
        return r;
    }

    /** Matrix-scalar multiplication **/
    template<typename ScalarType>
    inline MatrixType<R, C, ValueType>
    operator*(const ScalarType & scalar) const
    {
        MatrixType<R, C, ValueType> result = static_cast<const MatrixType<R, C, ValueType> &> (*this);
        std::transform(std::begin(*this), std::end(*this), std::begin(result), [scalar] (const ValueType & component) {
            return component*scalar;
        });
        return result;
    }

    /** Matrix multiplication-assignment **/
    template<typename ScalarType>
    inline MatrixType<R, C, ValueType> &
    operator*= (const ScalarType & scalar)
    {
        std::transform(std::begin(*this), std::end(*this), std::begin(*this), [scalar] (const ValueType & component) {
            return component*scalar;
        });
        return static_cast<MatrixType<R, C, ValueType> &> (*this);
    }

    /** Matrix-scalar division **/
    template<typename ScalarType>
    inline MatrixType<R, C, ValueType>
    operator/(const ScalarType & scalar) const
    {
        MatrixType<R, C, ValueType> result = static_cast<const MatrixType<R, C, ValueType> &> (*this);
        std::transform(std::begin(*this), std::end(*this), std::begin(result), [scalar] (const ValueType & component) {
            return component/scalar;
        });
        return result;
    }

    /** Matrix division-assignment **/
    template<typename ScalarType>
    inline MatrixType<R, C, ValueType> &
    operator/= (const ScalarType & scalar)
    {
        std::transform(std::begin(*this), std::end(*this), std::begin(*this), [scalar] (const ValueType & component) {
            return component/scalar;
        });
        return static_cast<MatrixType<R, C, ValueType> &> (*this);
    }

    /** Matrix addition **/
    template<typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    operator+ (const MatrixType<R, C, OtherValueType> & other) const
    {
        MatrixType<R, C, ValueType> result = static_cast<const MatrixType<R, C, ValueType> &> (*this);
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(result), std::plus<ValueType>());
        return result;
    }

    /** Matrix addition-assignment **/
    template<typename OtherValueType>
    inline MatrixType<R, C, ValueType> &
    operator+= (const MatrixType<R, C, OtherValueType> & other)
    {
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(*this), std::plus<ValueType>());
        return static_cast<MatrixType<R, C, ValueType> &> (*this);
    }

    /** Matrix subtraction **/
    template<typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    operator- (const MatrixType<R, C, OtherValueType> & other) const
    {
        MatrixType<R, C, ValueType> result;
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(result), std::minus<ValueType>());
        return result;
    }

    /** Matrix subtraction-assignment **/
    template<typename OtherValueType>
    inline MatrixType<R, C, ValueType> &
    operator-= (const MatrixType<R, C, OtherValueType> & other)
    {
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(*this), std::minus<ValueType>());
        return static_cast<MatrixType<R, C, ValueType> &> (*this);
    }


    /////////////////////////////////
    //// Mathematical operations ////
    /////////////////////////////////

    /** Get the transposed matrix */
    inline MatrixType<C, R, ValueType>
    T() const
    {
        const auto self = static_cast<const MatrixType<R, C, ValueType> &> (*this);
        return self.transposed();
    }

    /** Get the transposed matrix */
    inline MatrixType<C, R, ValueType>
    transposed() const
    {
        MatrixType<C, R, ValueType> Mt;
        for (size_t row = 0; row < C; ++row) {
            for (size_t column = 0; column < R; ++column) {
                Mt[row*R + column] = (*this)(column, row);
            }
        }

        return Mt;
    }

    /**
     * Compute the direct sum with another vector (sum between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_summation(const MatrixType<R, C, OtherValueType> & other) const
    {
        MatrixType<R, C, ValueType> result(false);
        std::transform(std::begin(*this), std::end(*this), std::begin(other), std::begin(result), std::plus<ValueType >());
        return result;
    }

    /** Alias to direct_summation */
    template <typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_sum(const MatrixType<R, C, OtherValueType> & other) const
    {
        return direct_summation(other);
    }

    /**
     * Compute the direct sub with another vector (subtraction between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_substraction(const MatrixType<R, C, OtherValueType> & other) const
    {
        MatrixType<R, C, ValueType> result(false);
        std::transform(std::begin(*this), std::end(*this), std::begin(other), std::begin(result), std::minus<ValueType >());
        return result;
    }

    /** Alias to direct_substraction */
    template <typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_sub(const MatrixType<R, C, OtherValueType> & other) const
    {
        return direct_substraction(other);
    }

    /**
     * Compute the direct multiplication with another vector (multiplication between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_multiplication(const MatrixType<R, C, OtherValueType> & other) const
    {
        MatrixType<R, C, ValueType> result(false);
        std::transform(std::begin(*this), std::end(*this), std::begin(other), std::begin(result), std::multiplies<ValueType >());
        return result;
    }

    /** Alias to direct_multiplication */
    template <typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_mult(const MatrixType<R, C, OtherValueType> & other) const
    {
        return direct_multiplication(other);
    }

    /**
     * Compute the direct division with another vector (division between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_division(const MatrixType<R, C, OtherValueType> & other) const
    {
        MatrixType<R, C, ValueType> result(false);
        std::transform(std::begin(*this), std::end(*this), std::begin(other), std::begin(result), std::divides<ValueType >());
        return result;
    }

    /** Alias to direct_division */
    template <typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_div(const MatrixType<R, C, OtherValueType> & other) const
    {
        return direct_division(other);
    }
};

/**
 * A simple representation of a matrix.
 *
 * This class can be extended by partial specializations for special matrices :
 *    Square matrix, rectangular matrix or vector matrix.
 *
 * It derives BaseMatrix to get its core functionalities (operations and functions that apply to all types of matrices).
 *
 * @tparam R_ The number of rows
 * @tparam C_ The number of columns
 * @tparam ValueType_ The data type of the matrix's components (default to float)
 */
template <size_t R, size_t C, typename ValueType = double>
struct Matrix : public BaseMatrix<Matrix, R, C, ValueType>
{
    using BaseMatrix<Matrix, R, C, ValueType>::BaseMatrix;
};

} // namespace algebra

} // namespace caribou

#include <Caribou/Algebra/SquareMatrix.h>
#include <Caribou/Algebra/Vector.h>

/**
 * Stream format of a matrix
 */
template <size_t R, size_t C, typename TComponent=FLOATING_POINT_TYPE>
inline std::ostream&
operator<<(std::ostream& os, const caribou::algebra::Matrix<R, C, TComponent>& m)
{
    const auto precision = os.precision();
    const auto width = os.width();
    os << std::string("[\n");
    for (size_t i = 0; i < R; ++i) {
        os << std::string("  [");

        for (size_t j = 0; j < C; ++j) {
            os.width(precision+1);
            os << m(i, j);
            os.width(width);
            if (j < C-1)
                os << ", ";
        }

        os << std::string("]");
        if (i < R-1)
            os << std::string(",");
        os << '\n';
    }
    os << std::string("]\n");
    return os;
}

/**
 * Multiplication of a scalar s with a matrix m : s * m
 */
template <size_t R, size_t C, typename TComponent, typename TScalar, typename std::enable_if<std::is_arithmetic<TScalar>::value, int>::type = 0>
inline caribou::algebra::Matrix<R, C, TComponent>
operator * (const TScalar & s, const caribou::algebra::Matrix<R, C, TComponent>& m)
{
    return m*s;
}

// The following templates will be instantiated in the algebra library. If you are using only those, you don't need to
// include the matrix.inl in your code, simply link the algebra library to your program/library.
//extern template struct Matrix<3, 3, FLOATING_POINT_TYPE>;
//extern template struct Matrix<2, 2, FLOATING_POINT_TYPE>;
//extern template struct Matrix<3, 1, FLOATING_POINT_TYPE>;
//extern template struct Matrix<1, 3, FLOATING_POINT_TYPE>;

#endif //CARIBOU_ALGEBRA_MATRIX_H
