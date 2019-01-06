#ifndef CARIBOU_ALGEBRA_INTERNAL_BASEMATRIX_H
#define CARIBOU_ALGEBRA_INTERNAL_BASEMATRIX_H

#include <array>

#include <Caribou/config.h>
#include <Caribou/Traits.h>

namespace caribou {
namespace algebra {
namespace internal {

// Can be used to statically detect if a BaseMatrix is base of any type
struct CaribouMatrix {};

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
template <template <size_t, size_t, typename> class MatrixType_, size_t R_, size_t C_, typename ValueType_>
struct BaseMatrix : public std::array<ValueType_, R_*C_>, public CaribouMatrix
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
     * Constructor by c-array or initializer list of another type (must be an integral or floating-point).
     * @example
     * \code{.cpp}
     * Matrix<3,3, float> A (
     *   {
     *     {1,2,3}, // Row 0
     *     {4,5,6}, // Row 1
     *     {7,8,9} //  Row 2
     *   }
     * );.
     * \endcode
     */
    template <typename OtherValueType, typename std::enable_if<std::is_integral<OtherValueType>::value or std::is_floating_point<OtherValueType>::value ,int>::type = 0>
    constexpr
    BaseMatrix(OtherValueType const (&components)[R][C]) {
        copy_from<0,0, OtherValueType>(components);
    }

    /**
     * Constructor by c-array or initializer list of another type (must be an integral or floating-point).
     * @example
     * \code{.cpp}
     * Matrix<3,3, float> A (
     *   {
     *     1,2,3, // Row 0
     *     4,5,6, // Row 1
     *     7,8,9 //  Row 2
     *   }
     * );.
     * \endcode
     */
    template <typename OtherValueType, typename std::enable_if<std::is_integral<OtherValueType>::value or std::is_floating_point<OtherValueType>::value ,int>::type = 0>
    constexpr
    BaseMatrix(OtherValueType const (&components)[R*C]) {
        copy_from<0,0, OtherValueType>(components);
    }

    /**
     * Copy constructor from another matrix of a different data type
     */
    template <template <size_t, size_t, typename> class OtherMatrixType, typename OtherValueType>
    constexpr
    BaseMatrix(const OtherMatrixType<R, C, OtherValueType> & other) {
        copy_from<0, ValueType> (other);
    }

    /** Constructor from a list of parameters (each parameter is a scalar component of the matrix) **/
    template<
            typename OtherValueType, typename ...Args,
            REQUIRES(R_*C_ == sizeof...(Args) + 1),
            REQUIRES(std::is_arithmetic_v<OtherValueType>)
    >
    constexpr
    BaseMatrix(OtherValueType first_value, Args&&...e)
    : std::array<ValueType, R_*C_> {{
        static_cast<ValueType>(first_value),
        static_cast<ValueType>(std::forward<Args>(e))...
    }}
    {}

    /** Constructor from a list of rows (each parameter is a vector of size Rx1) **/
    template <
            template <size_t, size_t, typename> class OtherMatrixType,
            typename OtherValueType,
            typename ...Args,
            REQUIRES(std::is_arithmetic_v<OtherValueType>),
            REQUIRES(sizeof...(Args)+1 == R_)
    >
    constexpr
    BaseMatrix(const OtherMatrixType<C_, 1, OtherValueType> & first_vector, Args&&...remaining_vectors)
    {
        copy_from<0>(std::make_index_sequence<C_>{}, first_vector, std::forward<Args>(remaining_vectors)...);
    }

    ///////////////////
    //// Accessors ////
    ///////////////////

    inline constexpr ValueType &
    operator () (const Index & row, const Index & column)
    {
        return static_cast<Self&>(*this) [row*C + column];
    }

    inline constexpr const ValueType &
    operator () (const Index & row, const Index & column) const
    {
        return static_cast<const Self &>(*this) [row*C + column];
    }

    inline
    MatrixType<1, C, ValueType>
    row(const Index & row) const
    {
        MatrixType<1, C, ValueType> r;
        for (std::size_t c = 0; c < C; ++c)
            r(0, c) = static_cast<const Self&>(*this) (row, c);
        return r;
    }

    inline
    MatrixType<R, 1, ValueType>
    column(const Index & column) const
    {
        MatrixType<R, 1, ValueType> c;
        for (std::size_t r = 0; r < R; ++r)
            c(r, 0) = static_cast<const Self&>(*this) (r, column);
        return c;
    }

    ///////////////////
    //// Operators ////
    ///////////////////

    /** Assignment operator **/
    template<typename OtherValueType>
    constexpr
    Self &
    operator = (const MatrixType<R, C, OtherValueType> & other)
    {
        copy_from<0> (other);
        return static_cast<Self&>(*this);
    }


    /** Comparison operator with a matrix of a same dimension and with component data type of OtherValueType **/
    template<typename OtherValueType>
    inline constexpr bool
    operator==(const MatrixType<R, C, OtherValueType> & other) const
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
     * Compute the direct summation with :
     *
     * - Same sized Matrix
     *
     * |A0 A1 A2|                     |B0 B1 B2|       |A0+B0 A1+B1 A2+B2|
     * |A3 A4 A5| . direct_summation( |B3 B4 B5| ) ==  |A3+B3 A4+B4 A5+B5|
     * |A6 A7 A8|                     |B6 B7 B8|       |A6+B6 A7+B7 A8+B8|
     *
     * - Row vector
     *
     * |A0 A1 A2|                                      |A0+B0 A1+B1 A2+B2|
     * |A3 A4 A5| . direct_summation( |B0 B1 B2| ) ==  |A3+B0 A4+B1 A5+B2|
     * |A6 A7 A8|                                      |A6+B0 A7+B1 A8+B2|
     *
     * - Column vector
     *
     * |A0 A1 A2|                     |B0|       |A0+B0 A1+B0 A2+B0|
     * |A3 A4 A5| . direct_summation( |B1| ) ==  |A3+B1 A4+B1 A5+B1|
     * |A6 A7 A8|                     |B2|       |A6+B2 A7+B2 A8+B2|
     *
     */
    template <size_t OtherR, size_t OtherC, typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_summation(const MatrixType<OtherR, OtherC, OtherValueType> & other) const
    {
        static_assert(
                (R == OtherR and C == OtherC) or (R == OtherR and OtherC == 1) or (C == OtherC and OtherR == 1),
                "Direct operations must be done with a same sized matrix, or with a column or row vector."
        );

        MatrixType<R, C, ValueType> result(false);

        if CONSTEXPR_IF (R == OtherR and C == OtherC) {
            // Direct operation with another matrix of the same size
            std::transform(std::begin(*this), std::end(*this), std::begin(other), std::begin(result), std::plus<ValueType >());

        } else if CONSTEXPR_IF(C == OtherC and OtherR == 1) {
            // Direct operation with a row-vector
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t c = 0; c < C; ++c) {
                    result(r,c) = static_cast<const Self&>(*this) (r, c) + other(0, c);
                }
            }
        } else if CONSTEXPR_IF(R == OtherR and OtherC == 1) {
            // Direct operation with a column-vector
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t c = 0; c < C; ++c) {
                    result(r,c) = static_cast<const Self&>(*this) (r, c) + other(r, 0);
                }
            }
        }

        return result;
    }

    /**
     * Compute the direct substraction with :
     *
     * - Same sized Matrix
     *
     * |A0 A1 A2|                        |B0 B1 B2|       |A0-B0 A1-B1 A2-B2|
     * |A3 A4 A5| . direct_substraction( |B3 B4 B5| ) ==  |A3-B3 A4-B4 A5-B5|
     * |A6 A7 A8|                        |B6 B7 B8|       |A6-B6 A7-B7 A8-B8|
     *
     * - Row vector
     *
     * |A0 A1 A2|                                         |A0-B0 A1-B1 A2-B2|
     * |A3 A4 A5| . direct_substraction( |B0 B1 B2| ) ==  |A3-B0 A4-B1 A5-B2|
     * |A6 A7 A8|                                         |A6-B0 A7-B1 A8-B2|
     *
     * - Column vector
     *
     * |A0 A1 A2|                        |B0|       |A0-B0 A1-B0 A2-B0|
     * |A3 A4 A5| . direct_substraction( |B1| ) ==  |A3-B1 A4-B1 A5-B1|
     * |A6 A7 A8|                        |B2|       |A6-B2 A7-B2 A8-B2|
     *
     */
    template <size_t OtherR, size_t OtherC, typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_substraction(const MatrixType<OtherR, OtherC, OtherValueType> & other) const
    {
        static_assert(
                (R == OtherR and C == OtherC) or (R == OtherR and OtherC == 1) or (C == OtherC and OtherR == 1),
                "Direct operations must be done with a same sized matrix, or with a column or row vector."
        );

        MatrixType<R, C, ValueType> result(false);

        if CONSTEXPR_IF (R == OtherR and C == OtherC) {
            // Direct operation with another matrix of the same size
            std::transform(std::begin(*this), std::end(*this), std::begin(other), std::begin(result), std::minus<ValueType >());

        } else if CONSTEXPR_IF(C == OtherC and OtherR == 1) {
            // Direct operation with a row-vector
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t c = 0; c < C; ++c) {
                    result(r,c) = static_cast<const Self&>(*this) (r, c) - other(0, c);
                }
            }
        } else if CONSTEXPR_IF(R == OtherR and OtherC == 1) {
            // Direct operation with a column-vector
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t c = 0; c < C; ++c) {
                    result(r,c) = static_cast<const Self&>(*this) (r, c) - other(r, 0);
                }
            }
        }

        return result;
    }

    /**
     * Compute the direct multiplication with :
     *
     * - Same sized Matrix
     *
     * |A0 A1 A2|                          |B0 B1 B2|       |A0*B0 A1*B1 A2*B2|
     * |A3 A4 A5| . direct_multiplication( |B3 B4 B5| ) ==  |A3*B3 A4*B4 A5*B5|
     * |A6 A7 A8|                          |B6 B7 B8|       |A6*B6 A7*B7 A8*B8|
     *
     * - Row vector
     *
     * |A0 A1 A2|                                           |A0*B0 A1*B1 A2*B2|
     * |A3 A4 A5| . direct_multiplication( |B0 B1 B2| ) ==  |A3*B0 A4*B1 A5*B2|
     * |A6 A7 A8|                                           |A6*B0 A7*B1 A8*B2|
     *
     * - Column vector
     *
     * |A0 A1 A2|                          |B0|       |A0*B0 A1*B0 A2*B0|
     * |A3 A4 A5| . direct_multiplication( |B1| ) ==  |A3*B1 A4*B1 A5*B1|
     * |A6 A7 A8|                          |B2|       |A6*B2 A7*B2 A8*B2|
     *
     */
    template <size_t OtherR, size_t OtherC, typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_multiplication(const MatrixType<OtherR, OtherC, OtherValueType> & other) const
    {
        static_assert(
                (R == OtherR and C == OtherC) or (R == OtherR and OtherC == 1) or (C == OtherC and OtherR == 1),
                "Direct operations must be done with a same sized matrix, or with a column or row vector."
        );

        MatrixType<R, C, ValueType> result(false);

        if CONSTEXPR_IF (R == OtherR and C == OtherC) {
            // Direct operation with another matrix of the same size
            std::transform(std::begin(*this), std::end(*this), std::begin(other), std::begin(result), std::multiplies<ValueType >());

        } else if CONSTEXPR_IF(C == OtherC and OtherR == 1) {
            // Direct operation with a row-vector
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t c = 0; c < C; ++c) {
                    result(r,c) = static_cast<const Self&>(*this) (r, c) * other(0, c);
                }
            }
        } else if CONSTEXPR_IF(R == OtherR and OtherC == 1) {
            // Direct operation with a column-vector
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t c = 0; c < C; ++c) {
                    result(r,c) = static_cast<const Self&>(*this) (r, c) * other(r, 0);
                }
            }
        }

        return result;
    }

    /**
     * Compute the direct division with :
     *
     * - Same sized Matrix
     *
     * |A0 A1 A2|                    |B0 B1 B2|       |A0/B0 A1/B1 A2/B2|
     * |A3 A4 A5| . direct_division( |B3 B4 B5| ) ==  |A3/B3 A4/B4 A5/B5|
     * |A6 A7 A8|                    |B6 B7 B8|       |A6/B6 A7/B7 A8/B8|
     *
     * - Row vector
     *
     * |A0 A1 A2|                                     |A0/B0 A1/B1 A2/B2|
     * |A3 A4 A5| . direct_division( |B0 B1 B2| ) ==  |A3/B0 A4/B1 A5/B2|
     * |A6 A7 A8|                                     |A6/B0 A7/B1 A8/B2|
     *
     * - Column vector
     *
     * |A0 A1 A2|                    |B0|       |A0/B0 A1/B0 A2/B0|
     * |A3 A4 A5| . direct_division( |B1| ) ==  |A3/B1 A4/B1 A5/B1|
     * |A6 A7 A8|                    |B2|       |A6/B2 A7/B2 A8/B2|
     *
     */
    template <size_t OtherR, size_t OtherC, typename OtherValueType>
    inline MatrixType<R, C, ValueType>
    direct_division(const MatrixType<OtherR, OtherC, OtherValueType> & other) const
    {
        static_assert(
                (R == OtherR and C == OtherC) or (R == OtherR and OtherC == 1) or (C == OtherC and OtherR == 1),
                "Direct operations must be done with a same sized matrix, or with a column or row vector."
        );

        MatrixType<R, C, ValueType> result(false);

        if CONSTEXPR_IF (R == OtherR and C == OtherC) {
            // Direct operation with another matrix of the same size
            std::transform(std::begin(*this), std::end(*this), std::begin(other), std::begin(result), std::divides<ValueType >());

        } else if CONSTEXPR_IF(C == OtherC and OtherR == 1) {
            // Direct operation with a row-vector
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t c = 0; c < C; ++c) {
                    result(r,c) = static_cast<const Self&>(*this) (r, c) / other(0, c);
                }
            }
        } else if CONSTEXPR_IF(R == OtherR and OtherC == 1) {
            // Direct operation with a column-vector
            for (std::size_t r = 0; r < R; ++r) {
                for (std::size_t c = 0; c < C; ++c) {
                    result(r,c) = static_cast<const Self&>(*this) (r, c) / other(r, 0);
                }
            }
        }

        return result;
    }

private:
    /** Static copying the values of another Matrix. The copy will be done entirely during the compilation by using template recursion. */
    template<size_t index, typename OtherValueType>
    constexpr
    void
    copy_from (const MatrixType<R, C, OtherValueType> & other)
    {
        static_assert(index >= 0,  "Cannot copy a value from an index outside of the matrix data array.");
        static_assert(index < R*C, "Cannot copy a value from an index outside of the matrix data array.");
        (*this)[index] = static_cast<ValueType> (other[index]);

        if CONSTEXPR_IF (index+1 < R*C)
            copy_from<index+1, OtherValueType>(other);
    }

    /** Static copying the values of RxC C-array. The copy will be done entirely during the compilation by using template recursion. */
    template<size_t row_index, size_t column_index, typename OtherValueType>
    constexpr
    void
    copy_from (OtherValueType const (&components)[R][C])
    {
        static_assert(row_index >= 0,     "Cannot copy a value from an index outside of the matrix data array.");
        static_assert(column_index >= 0,  "Cannot copy a value from an index outside of the matrix data array.");
        static_assert(row_index < R,      "Cannot copy a value from an index outside of the matrix data array.");
        static_assert(column_index < C,   "Cannot copy a value from an index outside of the matrix data array.");

        (*this)[C*row_index + column_index] = static_cast<ValueType> (components[row_index][column_index]);

        if CONSTEXPR_IF (column_index+1 < C)
            copy_from<row_index, column_index+1, OtherValueType>(components);
        else if CONSTEXPR_IF (row_index+1 < R)
            copy_from<row_index+1, 0, OtherValueType>(components);
    }

    /** Static copying the values of RxC C-array. The copy will be done entirely during the compilation by using template recursion. */
    template<size_t row_index, size_t column_index, typename OtherValueType>
    constexpr
    void
    copy_from (OtherValueType const (&components)[R*C])
    {
        static_assert(row_index >= 0,     "Cannot copy a value from an index outside of the matrix data array.");
        static_assert(column_index >= 0,  "Cannot copy a value from an index outside of the matrix data array.");
        static_assert(row_index < R,      "Cannot copy a value from an index outside of the matrix data array.");
        static_assert(column_index < C,   "Cannot copy a value from an index outside of the matrix data array.");

        (*this)[C*row_index + column_index] = static_cast<ValueType> (components[C*row_index + column_index]);

        if CONSTEXPR_IF (column_index+1 < C)
            copy_from<row_index, column_index+1, OtherValueType>(components);
        else if CONSTEXPR_IF (row_index+1 < R)
            copy_from<row_index+1, 0, OtherValueType>(components);
    }

    /** Static copying from a set of rows (vectors of size Cx1) */
    template <
            size_t row_id,
            std::size_t... Ix,
            template <size_t, size_t, typename> class OtherMatrixType,
            typename OtherValueType,
            typename ...Args
    >
    constexpr
    void
    copy_from( const std::index_sequence<Ix...> & indices, const OtherMatrixType<C_, 1, OtherValueType> & current_vector, Args&&...remaining_vectors)
    {
        (void ((*this)[C*row_id + Ix] = static_cast<ValueType> (current_vector[Ix])), ...);
        if constexpr (row_id< R-1)
            copy_from<row_id+1>(indices, std::forward<Args>(remaining_vectors)...);
    }
};

} // namespace internal
} // namespace algebra
} // namespace caribou

#endif //CARIBOU_ALGEBRA_INTERNAL_BASEMATRIX_H
