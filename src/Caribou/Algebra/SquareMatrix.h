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
#include <Caribou/Algebra/Internal/BaseSquareMatrix.h>

namespace caribou {
namespace algebra {

/**
 * Generic square NxN matrix.
 */
 template <size_t N, typename ValueType>
struct Matrix<N,N, ValueType> : public internal::BaseSquareMatrix<Matrix, N, ValueType>
{
    using internal::BaseSquareMatrix<Matrix, N, ValueType>::BaseSquareMatrix;
};

/**
 * Square 2x2 matrix.
 */
template <typename ValueType>
struct Matrix<2,2, ValueType> : public internal::BaseSquareMatrix<Matrix, 2, ValueType>
{
    using internal::BaseSquareMatrix<Matrix, 2, ValueType>::BaseSquareMatrix;

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
struct Matrix<3,3, ValueType> : public internal::BaseSquareMatrix<Matrix, 3, ValueType>
{
    using internal::BaseSquareMatrix<Matrix, 3, ValueType>::BaseSquareMatrix;

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