#include "LDLTSolver.h"
#include <SofaCaribou/Solver/EigenSparseSolver.inl>

#include <algorithm>
#include <cctype>
#include <string>

#include<Eigen/SparseCholesky>
#ifdef CARIBOU_WITH_MKL

// Bug introduced in Eigen 3.3.8, fixed in bfdd4a9
#ifndef EIGEN_USING_STD
#define EIGEN_USING_STD(a) EIGEN_USING_STD_MATH(a)
#endif

#include <Eigen/PardisoSupport>
#endif

#include <sofa/core/ObjectFactory.h>

namespace SofaCaribou::solver {


namespace internal {
template<typename MatrixType, int UpLo, typename Ordering>
struct solver_traits <Eigen::SimplicialLDLT<MatrixType,UpLo,Ordering>> {
    static auto BackendName() -> std::string {return "Eigen";}
};

#ifdef CARIBOU_WITH_MKL
template<typename MatrixType, int UpLo>
struct solver_traits <Eigen::PardisoLDLT< MatrixType, UpLo >> {
    static auto BackendName() -> std::string {return "Pardiso";}
};
#endif
}

template<typename EigenSolver>
LDLTSolver<EigenSolver>::LDLTSolver()
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


    // Put the backend name in lower case
    std::string backend_str = internal::solver_traits<EigenSolver>::BackendName();
    std::transform(backend_str.begin(), backend_str.end(), backend_str.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> backend = d_backend;
    if (backend_str == "pardiso") { // Case insensitive
        backend->setSelectedItem(static_cast<unsigned int>(1));
    } else {
        backend->setSelectedItem(static_cast<unsigned int>(0));
    }
}

static int SparseLDLTSolverClass = sofa::core::RegisterObject("Caribou Sparse LDLT linear solver")
    .add< LDLTSolver<Eigen::SimplicialLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>, Eigen::Lower, Eigen::AMDOrdering<int>>> >(true)
#ifdef CARIBOU_WITH_MKL
    .add< LDLTSolver<Eigen::PardisoLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>> >()
#endif
;

} // namespace SofaCaribou::solver