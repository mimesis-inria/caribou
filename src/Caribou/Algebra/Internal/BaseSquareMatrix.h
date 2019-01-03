#ifndef CARIBOU_ALGEBRA_INTERNAL_BASESQUAREMATRIX_H
#define CARIBOU_ALGEBRA_INTERNAL_BASESQUAREMATRIX_H

#include <Caribou/config.h>
#include <Caribou/Algebra/Internal/BaseMatrix.h>

namespace caribou {
namespace algebra {
namespace internal {

/**
 * Square RxR matrix.
 *
 * ** Do not use this class directly. Use instead caribou::algebra::Matrix<R, C> with R == C. **
 *
 * This class extend the BaseMatrix by adding functions only available on square matrices. It is meant to be derived by
 * a partial specialization of the class Matrix<R, C> (see later in this file for the 2x2 and 3x3 implementations)
 *
 * The functions declared in this class ca be used with any type of square matrices (2x2, 3x3, ...).
 *
 * @tparam MatrixType_ <R_, C_, ValueType_> See BaseMatrix template MatrixType_.
 * @tparam R_ The number of rows and columns
 * @tparam ValueType_ The data type of the matrix's components (default to float)
 */
template <template <size_t, size_t, typename> class MatrixType_, size_t R_, typename ValueType_>
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

} // namespace internal
} // namespace algebra
} // namespace caribou

#endif //CARIBOU_ALGEBRA_INTERNAL_BASESQUAREMATRIX_H
