#include "SparseLLTSolver.h"
#include <SofaCaribou/Solver/EigenSparseSolver.inl>

#include<Eigen/SparseCholesky>

#ifdef CARIBOU_WITH_MKL
#include <Eigen/PardisoSupport>
#endif

#include <sofa/core/ObjectFactory.h>

namespace SofaCaribou::solver {

namespace internal {

template<typename MatrixType, int UpLo, typename Ordering>
struct solver_traits<Eigen::SimplicialLLT<MatrixType, UpLo, Ordering>> {
    static auto TemplateName() -> std::string { return "SimplicialLLT"; }

    static auto BackendName() -> std::string { return "Eigen"; }
};

#ifdef CARIBOU_WITH_MKL

template<typename MatrixType, int UpLo>
struct solver_traits<Eigen::PardisoLLT<MatrixType, UpLo >> {
    static auto TemplateName() -> std::string { return "PardisoLLT"; }

    static auto BackendName() -> std::string { return "Pardiso"; }
};

#endif
}

template<typename EigenSolver>
SparseLLTSolver<EigenSolver>::SparseLLTSolver()
    : d_backend(initData(&d_backend,
        "backend",
        R"(
            Solver backend used.

            Available backends are:
                Eigen:   Eigen LLT solver (SimplicialLLT) [default].
                Pardiso: Pardiso LLT solver.
        )",
        true /*displayed_in_GUI*/, true /*read_only_in_GUI*/))
{
    d_backend.setValue(sofa::helper::OptionsGroup(std::vector<std::string> {
        "Eigen", "Pardiso"
    }));

    sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> backend = d_backend;
    if (internal::solver_traits<EigenSolver>::TemplateName() == "PardisoLLT") {
        backend->setSelectedItem(static_cast<unsigned int>(1));
    } else {
        backend->setSelectedItem(static_cast<unsigned int>(0));
    }
}

static int SparseLLTSolverClass = sofa::core::RegisterObject("Caribou Sparse LLT linear solver")
    .add< SparseLLTSolver<Eigen::SimplicialLLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor>>> >(true)
#ifdef CARIBOU_WITH_MKL
    .add< SparseLLTSolver<Eigen::PardisoLLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, long long int>>> >()
#endif
;

} // namespace SofaCaribou::solver