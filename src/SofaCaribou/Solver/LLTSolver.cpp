#include <SofaCaribou/Solver/LLTSolver.inl>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/core/ObjectFactory.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::solver {

static int SparseLLTSolverClass = sofa::core::RegisterObject("Caribou Sparse LLT linear solver")
    .add< LLTSolver<Eigen::SimplicialLLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>, Eigen::Lower, Eigen::AMDOrdering<int>>> >(true)
#ifdef CARIBOU_WITH_MKL
    .add< LLTSolver<Eigen::PardisoLLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>> >()
#endif
;

} // namespace SofaCaribou::solver