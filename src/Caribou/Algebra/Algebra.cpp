
#include <Caribou/Algebra/Matrix.inl>

namespace caribou
{
namespace algebra
{

template struct Matrix<3, 3, true, FLOATING_POINT_TYPE>;
template struct Matrix<3, 3, false, FLOATING_POINT_TYPE>;

template struct Matrix<2, 2, true, FLOATING_POINT_TYPE>;
template struct Matrix<2, 2, false, FLOATING_POINT_TYPE>;

template struct Matrix<3, 1, true, FLOATING_POINT_TYPE>;
template struct Matrix<3, 1, false, FLOATING_POINT_TYPE>;

template struct Matrix<1, 3, true, FLOATING_POINT_TYPE>;
template struct Matrix<1, 3, false, FLOATING_POINT_TYPE>;

} // namespace algebra

} // namespace caribou
