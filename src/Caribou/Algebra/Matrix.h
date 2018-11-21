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
 * @tparam N_ The number of lines
 * @tparam M_ The number of columns
 * @tparam TransposedData_ If true, the data will be stored by columns ([c1, c2, ..., cm]) instead of rows ([r1, r2, ..., rn]).
 *                        Set this to true if you are doing a lot of matrix*vector operations.
 * @tparam ValueType_ The data type of the matrix's components (default to float)
 */
template <size_t N_, size_t M_, bool TransposedData_ = true, typename ValueType_=FLOATING_POINT_TYPE>
struct Matrix : public std::array<ValueType_, N_*M_>
{
    static constexpr size_t N = N_; ///< Number of rows
    static constexpr size_t M = M_; ///< Number of columns per line
    static constexpr bool TransposedData = TransposedData_; ///< Internal data stored by columns

    using ValueType = ValueType_;
    using Self = Matrix<N, M, TransposedData, ValueType>;
    using Line = Matrix<1, M, TransposedData, ValueType>;
    using Column = Matrix<N, 1, TransposedData, ValueType>;
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
    Matrix(OtherValueType const (&components)[N][M]) {
        if (TransposedData) {
            for (size_t column = 0; column < M; ++column) {
                for (size_t row = 0; row < N; ++row) {
                    (*this)[column*N + row] = static_cast<ValueType>(components[row][column]);
                }
            }
        } else {
            for (size_t row = 0; row < N; ++row) {
                for (size_t column = 0; column < M; ++column) {
                    (*this)[row*M + column] = static_cast<ValueType>(components[row][column]);
                }
            }
        }
    }

    /**
     * Copy constructor from another matrix of a different data type
     */
    template <typename OtherValueType, bool OtherTransposedData>
    Matrix(const Matrix<N, M, OtherTransposedData, OtherValueType> & other) {
        if (TransposedData) {
            for (size_t column = 0; column < M; ++column) {
                for (size_t row = 0; row < N; ++row) {
                    (*this)[column*N + row] = static_cast<ValueType>(other(row, column));
                }
            }
        } else {
            for (size_t row = 0; row < N; ++row) {
                for (size_t column = 0; column < M; ++column) {
                    (*this)[row*M + column] = static_cast<ValueType>(other(row, column));
                }
            }
        }
    }

    ///////////////////
    //// Accessors ////
    ///////////////////

    inline ValueType &
    operator () (const Index & row, const Index & column)
    {
        return (TransposedData ? (*this)[column*N + row] : (*this)[row*M + column]);
    }

    inline const ValueType &
    operator () (const Index & row, const Index & column) const
    {
        return (TransposedData ? (*this)[column*N + row] : (*this)[row*M + column]);
    }

    ///////////////////
    //// Operators ////
    ///////////////////

    /** Comparison operator with a matrix of a same dimension and with component data type of OtherValueType **/
    template<bool OtherTransposedData, typename OtherValueType>
    constexpr bool
    operator==(const Matrix<M, N, OtherTransposedData, OtherValueType> & other) const
    {
        if (TransposedData == OtherTransposedData)
            return std::equal(this->begin(), this->end(), other.begin());
        else
            return std::equal(this->begin(), this->end(), other.transposed().begin());
    }

    /////////////////////////////////
    //// Mathematical operations ////
    /////////////////////////////////

    /** Get the transposed matrix */
    inline Matrix<M, N, TransposedData, ValueType>
    transposed() const
    {
        Matrix<M, N, TransposedData, ValueType> Mt;
        if (TransposedData) {
            for (size_t column = 0; column < N; ++column) {
                for (size_t row = 0; row < M; ++row) {
                    Mt[column*M + row] = (*this)(column, row);
                }
            }
        } else {
            for (size_t row = 0; row < M; ++row) {
                for (size_t column = 0; column < N; ++column) {
                    Mt[row*N + column] = (*this)(column, row);
                }
            }
        }

        return Mt;
    }

    /** Returns the identity matrix */
    static Matrix<N,M,TransposedData, ValueType> Identity()
    {
        Matrix<N,M,TransposedData, ValueType> id(true /*initialize to zero*/);

        if (N != M)
            return id;

        for (Index i = 0; i < M; i++)
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
template <size_t N_, size_t M_, bool TransposedData_ = true, typename ValueType_=FLOATING_POINT_TYPE>
Matrix<N_, M_, TransposedData_, ValueType_> inverse(const Matrix<N_, M_, TransposedData_, ValueType_> & m);

/**
 * Compute the determinant of a matrix
 */
template <size_t N_, size_t M_, bool TransposedData_ = true, typename ValueType_=FLOATING_POINT_TYPE>
ValueType_ determinant(const Matrix<N_, M_, TransposedData_, ValueType_> & m);

// The following templates will be instantiated in the algebra library. If you are using only those, you don't need to
// include the matrix.inl in your code, simply link the algebra library to your program/library.
extern template struct Matrix<3, 3, true, FLOATING_POINT_TYPE>;
extern template struct Matrix<3, 3, false, FLOATING_POINT_TYPE>;

extern template struct Matrix<2, 2, true, FLOATING_POINT_TYPE>;
extern template struct Matrix<2, 2, false, FLOATING_POINT_TYPE>;

extern template struct Matrix<3, 1, true, FLOATING_POINT_TYPE>;
extern template struct Matrix<3, 1, false, FLOATING_POINT_TYPE>;

extern template struct Matrix<1, 3, true, FLOATING_POINT_TYPE>;
extern template struct Matrix<1, 3, false, FLOATING_POINT_TYPE>;

} // namespace algebra

} // namespace caribou

#endif //CARIBOU_ALGEBRA_MATRIX_H
