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
#include <Caribou/Algebra/Internal/BaseMatrix.h>

namespace caribou {
namespace algebra {



/**
 * Generic RxC matrix.
 *
 * This class is the generic implementation, which means that any R and C values that aren't specialized elsewhere will
 * fall into this class implementation. It derives BaseMatrix to get its core functionalities (operations and functions
 * that apply to all types of matrices).
 *
 * In the following tree, the Matrix<R, C> are specific implementations of this base class following the R and C template
 * arguments choosen.
 *
 * BaseMatrix
 * ├── Matrix<R, C> ................Implementation of general matrix RxC
 * ├── BaseSquareMatrix
 * │         ├─ Matrix<R, R> .......Implementation of square matrix RxR
 * │         ├─ Matrix<2, 2> .......Implementation of square matrix 2x2
 * │         └─ Matrix<3, 3> .......Implementation of square matrix 3x3
 * └── BaseVector
 *           ├─Matrix<R, 1> ........Implementation of row-vector matrix Rx1
 *           ├─Matrix<1, C> ........Implementation of col-vector matrix 1xC
 *           ├─Matrix<1, 1> ........Implementation of unit matrix 1x1
 *           └─BaseVector3D
 *                ├─ Matrix<3, 1> ..Implementation of row-vector matrix 3x1
 *                └─ Matrix<1, 3> ..Implementation of col-vector matrix 1x3
 *
 * @tparam R_ The number of rows
 * @tparam C_ The number of columns
 * @tparam ValueType_ The data type of the matrix's components (default to float)
 */
template <size_t R, size_t C, typename ValueType = double>
struct Matrix : public internal::BaseMatrix<Matrix, R, C, ValueType>
{
    using internal::BaseMatrix<Matrix, R, C, ValueType>::BaseMatrix;
};

} // namespace algebra
} // namespace caribou

#include <Caribou/Algebra/SquareMatrix.h>

/**
 * Output stream of a generic RxC matrix.
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
 * Multiplication of a scalar s with a generic RxC matrix m : s * m
 */
template <size_t R, size_t C, typename TComponent, typename TScalar, typename std::enable_if<std::is_arithmetic<TScalar>::value, int>::type = 0>
inline caribou::algebra::Matrix<R, C, TComponent>
operator * (const TScalar & s, const caribou::algebra::Matrix<R, C, TComponent>& m)
{
    return m*s;
}

#endif //CARIBOU_ALGEBRA_MATRIX_H
