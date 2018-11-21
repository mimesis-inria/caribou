#ifndef CARIBOU_ALGEBRA_MATRIX_INL
#define CARIBOU_ALGEBRA_MATRIX_INL

#include <Caribou/config.h>
#include <Caribou/Algebra/Matrix.h>

namespace caribou
{
namespace algebra
{


template <size_t R_, size_t C_, typename ValueType_=FLOATING_POINT_TYPE>
Matrix<R_, C_, ValueType_>
inverse(const Matrix<R_, C_, ValueType_> & /* m */)
{
    throw std::logic_error("Inverse of Matrix " + std::to_string(R_) + "x" + std::to_string(C_) + " isn't implemented.");
}

template <typename ValueType_=FLOATING_POINT_TYPE>
Matrix<2, 2, ValueType_>
inverse(const Matrix<2, 2, ValueType_> & from)
{
    Matrix<2, 2, ValueType_> dest( true /* initialize_to_zero */);

    ValueType_ det=determinant(from);

    if ( -EPSILON<=det && det<=EPSILON)
    {
        throw std::overflow_error("Trying to inverse a singular 2x2 matrix (determinant is close to zero).");
    }

    dest(0,0)=  from(1,1)/det;
    dest(0,1)= -from(0,1)/det;
    dest(1,0)= -from(1,0)/det;
    dest(1,1)=  from(0,0)/det;

    return dest;
}

template <typename ValueType_=FLOATING_POINT_TYPE>
Matrix<3, 3, ValueType_>
inverse(const Matrix<3, 3, ValueType_> & from)
{
    Matrix<3, 3, ValueType_> dest( true /* initialize_to_zero */);

    ValueType_ det=determinant(from);

    if ( -EPSILON<=det && det<=EPSILON)
    {
        throw std::overflow_error("Trying to inverse a singular 3x3 matrix (determinant is close to zero).");
    }

    dest(0,0)= (from(1,1)*from(2,2) - from(2,1)*from(1,2))/det;
    dest(1,0)= (from(1,2)*from(2,0) - from(2,2)*from(1,0))/det;
    dest(2,0)= (from(1,0)*from(2,1) - from(2,0)*from(1,1))/det;
    dest(0,1)= (from(2,1)*from(0,2) - from(0,1)*from(2,2))/det;
    dest(1,1)= (from(2,2)*from(0,0) - from(0,2)*from(2,0))/det;
    dest(2,1)= (from(2,0)*from(0,1) - from(0,0)*from(2,1))/det;
    dest(0,2)= (from(0,1)*from(1,2) - from(1,1)*from(0,2))/det;
    dest(1,2)= (from(0,2)*from(1,0) - from(1,2)*from(0,0))/det;
    dest(2,2)= (from(0,0)*from(1,1) - from(1,0)*from(0,1))/det;

    return dest;
}

template <size_t R_, size_t C_, typename ValueType_=FLOATING_POINT_TYPE>
ValueType_
determinant(const Matrix<R_, C_, ValueType_> & /* m */)
{
    throw std::logic_error("Determinant of Matrix " + std::to_string(R_) + "x" + std::to_string(C_) + " isn't implemented.");
}

template <typename ValueType_=FLOATING_POINT_TYPE>
ValueType_
determinant(const Matrix<2, 2, ValueType_> & m)
{
    return m(0,0)*m(1,1) - m(1,0)*m(0,1);
}

template <typename ValueType_=FLOATING_POINT_TYPE>
ValueType_
determinant(const Matrix<3, 3, ValueType_> & m)
{
    return   m(0,0)*m(1,1)*m(2,2)
           + m(1,0)*m(2,1)*m(0,2)
           + m(2,0)*m(0,1)*m(1,2)
           - m(0,0)*m(2,1)*m(1,2)
           - m(1,0)*m(0,1)*m(2,2)
           - m(2,0)*m(1,1)*m(0,2);
}

template <size_t R_, size_t C_, typename ValueType_>
typename Matrix<R_, C_, ValueType_>::Self
Matrix<R_, C_, ValueType_>::inverted() const
{
    return caribou::algebra::inverse(*this);
}

template <size_t R_, size_t C_, typename ValueType_>
ValueType_
Matrix<R_, C_, ValueType_>::determinant() const
{
    return caribou::algebra::determinant(*this);
}

} // namespace algebra

} // namespace caribou

#endif //CARIBOU_ALGEBRA_MATRIX_INL
