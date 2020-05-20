#include "SparseLDLTSolver.h"
#include <SofaCaribou/Solver/EigenSparseSolver.inl>

#include<Eigen/SparseCholesky>

#ifdef CARIBOU_WITH_MKL
#include <Eigen/PardisoSupport>
#endif

#include <sofa/core/ObjectFactory.h>

namespace SofaCaribou::solver {


namespace internal {
template<typename MatrixType, int UpLo, typename Ordering>
struct solver_traits <Eigen::SimplicialLDLT<MatrixType,UpLo,Ordering>> {
    static auto TemplateName() -> std::string {return "SimplicialLDLT";}
    static auto BackendName() -> std::string {return "Eigen";}
};

#ifdef CARIBOU_WITH_MKL
template<typename MatrixType, int UpLo>
struct solver_traits <Eigen::PardisoLDLT< MatrixType, UpLo >> {
    static auto TemplateName() -> std::string {return "PardisoLDLT";}
    static auto BackendName() -> std::string {return "Pardiso";}
};
#endif
}

template<typename EigenSolver>
SparseLDLTSolver<EigenSolver>::SparseLDLTSolver()
    : d_backend(initData(&d_backend,
        "backend",
        R"(
            Solver backend used.

            Available backends are:
                Eigen:   Eigen LDLT solver (SimplicialLDLT) [default].
                Pardiso: Pardiso LDLT solver.
        )",
        true /*displayed_in_GUI*/, true /*read_only_in_GUI*/))
{
    d_backend.setValue(sofa::helper::OptionsGroup(std::vector<std::string> {
        "Eigen", "Pardiso"
    }));

    sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> backend = d_backend;
    if (internal::solver_traits<EigenSolver>::TemplateName() == "PardisoLDLT") {
        backend->setSelectedItem(static_cast<unsigned int>(1));
    } else {
        backend->setSelectedItem(static_cast<unsigned int>(0));
    }
}

static int SparseLDLTSolverClass = sofa::core::RegisterObject("Caribou Sparse LDLT linear solver")
    .add< SparseLDLTSolver<Eigen::SimplicialLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor>>> >(true)
#ifdef CARIBOU_WITH_MKL
    .add< SparseLDLTSolver<Eigen::PardisoLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, long long int>>> >()
#endif
;

} // namespace SofaCaribou::solver