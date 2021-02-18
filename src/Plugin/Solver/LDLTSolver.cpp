#include <SofaCaribou/Solver/LDLTSolver.inl>

#include <sofa/core/ObjectFactory.h>

namespace SofaCaribou::solver {

static int SparseLDLTSolverClass = sofa::core::RegisterObject("Caribou Sparse LDLT linear solver")
    .add< LDLTSolver<Eigen::SimplicialLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>, Eigen::Lower, Eigen::AMDOrdering<int>>> >(true)
#ifdef CARIBOU_WITH_MKL
    .add< LDLTSolver<Eigen::PardisoLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>> >()
#endif
;

} // namespace SofaCaribou::solver