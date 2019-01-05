#ifndef CARIBOU_ALGEBRA_MATRIX_H
#define CARIBOU_ALGEBRA_MATRIX_H

#include <ostream>
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
template <size_t R, size_t C, typename ValueType = FLOATING_POINT_TYPE>
struct Matrix : public internal::BaseMatrix<Matrix, R, C, ValueType>
{
    using Base = internal::BaseMatrix<Matrix, R, C, ValueType>;
    using Base::Base;

    /**
     * Constructor by c-array or initializer list.
     * @example
     * \code{.cpp}
     * Matrix<3,3> A (
     *   {
     *     {1,2,3}, // Row 0
     *     {4,5,6}, // Row 1
     *     {7,8,9} //  Row 2
     *   }
     * );
     * \endcode
     */
    constexpr
    Matrix(ValueType const (&components)[R][C]) : Base (components) {}

    /**
     * Constructor by c-array or initializer list.
     * @example
     * \code{.cpp}
     * Matrix<3,3> A (
     *   {
     *     1,2,3, // Row 0
     *     4,5,6, // Row 1
     *     7,8,9 //  Row 2
     *   }
     * );
     * \endcode
     */
    constexpr
    Matrix(ValueType const (&components)[R*C]) : Base (components) {}

    /**
     * Copy constructor from another matrix of the same data type
     */
    constexpr
    Matrix(const Matrix<R, C, ValueType> & other) : Base(other) {}

    /** Constructor from a list of parameters (each parameter is a scalar component of the matrix) **/
    template<
            typename ...Args,
            REQUIRES(R*C == sizeof...(Args) + 1)
    >
    constexpr
    Matrix(ValueType first_value, Args&&...e)
    : Base(first_value, std::forward<Args>(e)...) {}
};


// Deduction guides

/** Constructor of a matrix RxC from a list of R rows (each parameter is a vector of size Cx1) **/
template <
        template <size_t, size_t, typename> class OtherMatrixType,
        size_t C,
        typename OtherValueType,
        typename ...Args,
        REQUIRES(std::is_arithmetic_v<OtherValueType>)
>
Matrix(const OtherMatrixType<C, 1, OtherValueType> & first_row, Args&&...rows) -> Matrix<sizeof...(rows)+1, C, OtherValueType>;

} // namespace algebra
} // namespace caribou

#include <Caribou/Algebra/SquareMatrix.h>

/**
 * Output stream of a generic RxC matrix.
 */
template <size_t R, size_t C, typename TComponent>
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
