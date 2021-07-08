#include <SofaCaribou/Solver/ConjugateGradientSolver.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

#include <Eigen/Sparse>

namespace SofaCaribou::solver {
int CGLinearSolverClass = sofa::core::RegisterObject("Linear system solver using the conjugate gradient iterative algorithm")
        .add< ConjugateGradientSolver<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>>(true);

template class ConjugateGradientSolver<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>;
} // namespace SofaCaribou::solver