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
#include <Caribou/Algebra/Vector.h>

namespace caribou
{
namespace algebra
{

/**
 * A simple representation of a matrix.
 *
 * This class is a minimum memory container of a matrix that offers many operations and tool functions. It is designed
 * to be used as a heap memory structure and is based on std::array for its storage.
 *
 *
 * @tparam R_ The number of rows
 * @tparam C_ The number of columns
 * @tparam ValueType_ The data type of the matrix's components (default to float)
 */
template <size_t R_, size_t C_, typename ValueType_=FLOATING_POINT_TYPE>
struct Matrix : public std::array<ValueType_, R_*C_>
{
    static constexpr size_t R = R_; ///< Number of rows
    static constexpr size_t C = C_; ///< Number of columns per row
    static constexpr size_t N = R*C; ///< Number of elements

    using ValueType = ValueType_;
    using Self = Matrix<R, C, ValueType>;
    using Row = Matrix<1, C, ValueType>;
    using Column = Matrix<R, 1, ValueType>;
    using Index = size_t;

    //////////////////////
    //// Constructors ////
    //////////////////////

    /**
     * Main constructor
     * @param initialize_to_zero If true, initialize the scalar components of the vector to zero
     */
    explicit
    Matrix(bool initialize_to_zero = false) {
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
    template <typename OtherValueType>
    Matrix(OtherValueType const (&components)[R][C]) {
        for (size_t row = 0; row < R; ++row) {
            for (size_t column = 0; column < C; ++column) {
                (*this)[row*C + column] = static_cast<ValueType>(components[row][column]);
            }
        }
    }

    /**
     * Copy constructor from another matrix of a different data type
     */
    template <typename OtherValueType>
    Matrix(const Matrix<R, C, OtherValueType> & other) {
        for (size_t row = 0; row < R; ++row) {
            for (size_t column = 0; column < C; ++column) {
                (*this)[row*C + column] = static_cast<ValueType>(other(row, column));
            }
        }
    }

    ///////////////////
    //// Accessors ////
    ///////////////////

    inline ValueType &
    operator () (const Index & row, const Index & column)
    {
        return (*this)[row*C + column];
    }

    inline const ValueType &
    operator () (const Index & row, const Index & column) const
    {
        return (*this)[row*C + column];
    }

    ///////////////////
    //// Operators ////
    ///////////////////

    /** Assignment operator **/
    template<typename OtherValueType>
    Matrix &
    operator = (const Matrix<R, C, OtherValueType> & other)
    {
        for (size_t i = 0; i < N; ++i)
            (*this)[i] = other[i];
        return (*this);
    }


    /** Comparison operator with a matrix of a same dimension and with component data type of OtherValueType **/
    template<typename OtherValueType>
    bool
    operator==(const Matrix<C, R, OtherValueType> & other) const
    {
        return std::equal(this->begin(), this->end(), other.begin());
    }

    /** Matrix multiplication **/
    template<size_t M, typename OtherValueType>
    Matrix<R, M, ValueType>
    operator*(const Matrix<C, M, OtherValueType> & other) const
    {
        Matrix<R, M, ValueType> r;
        for (size_t i = 0; i < R; ++i) {
            for (size_t j = 0; j < M; ++j) {
                r(i, j)=(*this)(i,0) * other(0, j);
                for(size_t k=1; k<C; k++)
                    r(i, j) += (*this)(i, k) * other(k, j);
            }
        }
        return r;
    }

    /** Matrix-vector multiplication **/
    template<typename OtherValueType>
    Matrix<R, 1, ValueType>
    operator*(const Vector<C, OtherValueType> & other) const
    {
        Matrix<R, 1, ValueType> r;
        for (size_t i = 0; i < R; ++i) {
            r[i] = (*this)(i,0) * other[0];
            for (size_t j = 1; j < C; ++j) {
                r[i]+=(*this)(i,j) * other[j];
            }
        }
        return r;
    }

    /** Matrix addition **/
    template<typename OtherValueType>
    Matrix<R, C, ValueType>
    operator+ (const Matrix<R, C, OtherValueType> & other) const
    {
        Matrix<R, C, ValueType> result;
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(result), std::plus<ValueType>());
        return result;
    }

    /** Matrix addition-assignment **/
    template<typename OtherValueType>
    Matrix<R, C, ValueType> &
    operator+= (const Matrix<R, C, OtherValueType> & other)
    {
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(*this), std::plus<ValueType>());
        return *this;
    }

    /** Matrix subtraction **/
    template<typename OtherValueType>
    Matrix<R, C, ValueType>
    operator- (const Matrix<R, C, OtherValueType> & other) const
    {
        Matrix<R, C, ValueType> result;
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(result), std::minus<ValueType>());
        return result;
    }

    /** Matrix subtraction-assignment **/
    template<typename OtherValueType>
    Matrix<R, C, ValueType> &
    operator-= (const Matrix<R, C, OtherValueType> & other)
    {
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(*this), std::minus<ValueType>());
        return *this;
    }

    /** Compute the exponent of a matrix. If the exponent is negative, we first inverse the matrix. **/
    template<typename Integer>
    Matrix<R, C, ValueType>
    operator^ (const Integer & exponent)
    {
        static_assert(std::is_integral<Integer>::value, "The exponent must be of integral type (ex. int, char, long, ...)");

        if (R != C)
            throw std::logic_error("Only the exponent of squared matrix can be computed.");

        if (exponent == 0)
            return Matrix<R, C, ValueType>::Identity();

        Integer N = (exponent < 0 ? -exponent : exponent);
        Matrix<R, C, ValueType> M = (exponent < 0 ? (*this).inverted() : (*this));
        Matrix<R, C, ValueType> res = M;

        for (Integer i = 1; i < N; ++i)
            res = res*M;

        return res;
    }



    /////////////////////////////////
    //// Mathematical operations ////
    /////////////////////////////////

    /** Get the transposed matrix */
    inline Matrix<C, R, ValueType>
    transposed() const
    {
        Matrix<C, R, ValueType> Mt;
        for (size_t row = 0; row < C; ++row) {
            for (size_t column = 0; column < R; ++column) {
                Mt[row*R + column] = (*this)(column, row);
            }
        }

        return Mt;
    }

    /** Returns the identity matrix */
    static Matrix<R,C, ValueType> Identity()
    {
        Matrix<R,C, ValueType> id(true /*initialize to zero*/);

        if (R != C)
            return id;

        for (Index i = 0; i < C; i++)
            id(i, i) = 1;

        return id;
    }

    /** Compute the determinant */
    ValueType
    determinant() const;

    /** Compute the inverted matrix */
    Self inverted() const;
};

/**
 * Compute the inverse of a matrix
 * @throws std::overflow_error when the matrix is singular
 */
template <size_t R_, size_t C_, typename ValueType_=FLOATING_POINT_TYPE>
Matrix<R_, C_, ValueType_> inverse(const Matrix<R_, C_, ValueType_> & m);

/**
 * Compute the determinant of a matrix
 */
template <size_t R_, size_t C_, typename ValueType_=FLOATING_POINT_TYPE>
ValueType_ determinant(const Matrix<R_, C_, ValueType_> & m);

// The following templates will be instantiated in the algebra library. If you are using only those, you don't need to
// include the matrix.inl in your code, simply link the algebra library to your program/library.
extern template struct Matrix<3, 3, FLOATING_POINT_TYPE>;
extern template struct Matrix<2, 2, FLOATING_POINT_TYPE>;
extern template struct Matrix<3, 1, FLOATING_POINT_TYPE>;
extern template struct Matrix<1, 3, FLOATING_POINT_TYPE>;

} // namespace algebra

} // namespace caribou

#endif //CARIBOU_ALGEBRA_MATRIX_H
