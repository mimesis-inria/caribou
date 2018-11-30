#ifndef CARIBOU_ALGEBRA_SQUAREMATRIX_H
#define CARIBOU_ALGEBRA_SQUAREMATRIX_H

#include <ostream>
#include <cstddef>
#include <array>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <cmath>

#include <Caribou/config.h>
#include <Caribou/Algebra/Matrix.h>

namespace caribou {
namespace algebra {

/**
 * Square RxR matrix.
 *
 * ** Do not use this class directly. Use instead caribou::algebra::Matrix<R, C> with R == C. **
 *
 * This class extend the BaseMatrix by adding functions only available on square matrices. It is meant to be derived by
 * a partial specialization of the class Matrix<R, C> (see later in this file for the 2x2 and 3x3 specializations)
 *
 * The functions declared in this class ca be used with any type of square matrices (2x2, 3x3, ...).
 *
 * @tparam MatrixType_ <R_, C_, ValueType_> See BaseMatrix template MatrixType_.
 * @tparam R_ The number of rows and columns
 * @tparam ValueType_ The data type of the matrix's components (default to float)
 */
template <template <size_t, size_t, typename> class MatrixType_, size_t R_, typename ValueType_=FLOATING_POINT_TYPE>
struct BaseSquareMatrix : public BaseMatrix<MatrixType_, R_, R_, ValueType_>
{
    static constexpr size_t R = R_; ///< Number of rows
    static constexpr size_t C = R_; ///< Number of columns per row
    static constexpr size_t N = R*C; ///< Number of elements

    template<size_t R__ = R, size_t C__ = C, typename ValueType__= ValueType_>
    using MatrixType = MatrixType_<R__, C__, ValueType__>;

    using ValueType = ValueType_;

    using Index = typename BaseMatrix<MatrixType_, R, R, ValueType>::Index;

    using BaseMatrix<MatrixType_, R, R, ValueType>::BaseMatrix;

    /** Returns the identity matrix */
    static MatrixType<R, R, ValueType> Identity()
    {
        MatrixType<R, R, ValueType> id(true /*initialize to zero*/);

        for (Index i = 0; i < R; i++)
            id(i, i) = 1;

        return id;
    }

    /** Compute the exponent of a matrix. If the exponent is negative, we first inverse the matrix. **/
    template<typename Integer>
    MatrixType<R, R, ValueType>
    operator^ (const Integer & exponent) const
    {
        static_assert(std::is_integral<Integer>::value, "The exponent must be of integral type (ex. int, char, long, ...)");

        if (exponent == 0)
            return MatrixType<R, R, ValueType>::Identity();

        auto self = static_cast<const MatrixType<R, R, ValueType>&>(*this);
        Integer N = (exponent < 0 ? -exponent : exponent);
        MatrixType<R, R, ValueType> M = (exponent < 0 ? self.inverted() : self);
        MatrixType<R, R, ValueType> res = M;

        for (Integer i = 1; i < N; ++i)
            res = res*M;

        return res;
    }
};

/**
 * Square NxN matrix.
 */
 template <size_t N, typename ValueType>
struct Matrix<N,N, ValueType> : public BaseSquareMatrix<Matrix, N, ValueType>
{
    using BaseSquareMatrix<Matrix, N, ValueType>::BaseSquareMatrix;
};

/**
 * Square 2x2 matrix.
 */
template <typename ValueType>
struct Matrix<2,2, ValueType> : public BaseSquareMatrix<Matrix, 2, ValueType>
{
    using BaseSquareMatrix<Matrix, 2, ValueType>::BaseSquareMatrix;

    /** Returns the identity matrix */
    static constexpr Matrix<2,2, ValueType> Identity()
    {
        return Matrix<2,2, ValueType>({{
            {1, 0},
            {0, 1}
        }});
    }

    /** Compute the determinant */
    inline constexpr ValueType
    determinant() const {
        return (*this)(0,0)*(*this)(1,1) - (*this)(1,0)*(*this)(0,1);
    }

    /**
     * Compute the inverted matrix
     *
     * @throws std::overflow_error when the matrix is singular (determinant is close to zero)
     */
    inline Matrix<2,2, ValueType>
    inverted() const {
        Matrix<2, 2, ValueType> dest( true /* initialize_to_zero */);
        ValueType det=determinant();

        if ( -EPSILON<=det && det<=EPSILON)
        {
            throw std::overflow_error("Trying to inverse a singular 2x2 matrix (determinant is close to zero).");
        }

        const Matrix<2,2, ValueType> & m = *this;
        dest(0,0)=  m(1,1)/det;
        dest(0,1)= -m(0,1)/det;
        dest(1,0)= -m(1,0)/det;
        dest(1,1)=  m(0,0)/det;

        return dest;
    }
};

/**
 * Square 3x3 matrix.
 */
template <typename ValueType>
struct Matrix<3,3, ValueType> : public BaseSquareMatrix<Matrix, 3, ValueType>
{
    using BaseSquareMatrix<Matrix, 3, ValueType>::BaseSquareMatrix;

    /** Returns the identity matrix */
    static constexpr Matrix<3,3, ValueType> Identity()
    {
        return Matrix<3,3, ValueType>({{
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        }});
    }

    /** Compute the determinant */
    inline constexpr ValueType
    determinant() const {
        return     (*this)(0,0)*(*this)(1,1)*(*this)(2,2)
                 + (*this)(1,0)*(*this)(2,1)*(*this)(0,2)
                 + (*this)(2,0)*(*this)(0,1)*(*this)(1,2)
                 - (*this)(0,0)*(*this)(2,1)*(*this)(1,2)
                 - (*this)(1,0)*(*this)(0,1)*(*this)(2,2)
                 - (*this)(2,0)*(*this)(1,1)*(*this)(0,2);
    }

    /**
     * Compute the inverted matrix
     *
     * @throws std::overflow_error when the matrix is singular (determinant is close to zero)
     */
    inline Matrix<3,3, ValueType>
    inverted() const {
        Matrix<3, 3, ValueType> dest( true /* initialize_to_zero */);
        ValueType det=determinant();

        if ( -EPSILON<=det && det<=EPSILON)
        {
            throw std::overflow_error("Trying to inverse a singular 3x3 matrix (determinant is close to zero).");
        }

        const Matrix<3,3, ValueType> & m = *this;
        dest(0,0)= (m(1,1)*m(2,2) - m(2,1)*m(1,2))/det;
        dest(1,0)= (m(1,2)*m(2,0) - m(2,2)*m(1,0))/det;
        dest(2,0)= (m(1,0)*m(2,1) - m(2,0)*m(1,1))/det;
        dest(0,1)= (m(2,1)*m(0,2) - m(0,1)*m(2,2))/det;
        dest(1,1)= (m(2,2)*m(0,0) - m(0,2)*m(2,0))/det;
        dest(2,1)= (m(2,0)*m(0,1) - m(0,0)*m(2,1))/det;
        dest(0,2)= (m(0,1)*m(1,2) - m(1,1)*m(0,2))/det;
        dest(1,2)= (m(0,2)*m(1,0) - m(1,2)*m(0,0))/det;
        dest(2,2)= (m(0,0)*m(1,1) - m(1,0)*m(0,1))/det;

        return dest;
    }
};

} // namespace algebra
} // namespace caribou

#endif //CARIBOU_ALGEBRA_SQUAREMATRIX_H