#include <SofaCaribou/Solver/EigenSolver.inl>

namespace SofaCaribou::solver {
template class EigenSolver<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>>;
template class EigenSolver<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>;
}
