#include "SparseLUSolver.h"
#include <SofaCaribou/Solver/EigenSparseSolver.inl>

#include<Eigen/SparseLU>

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
template<typename MatrixType, typename Ordering>
struct solver_traits <Eigen::SparseLU<MatrixType,Ordering>> {
    static auto TemplateName() -> std::string {return "SparseLU";}
    static auto BackendName() -> std::string {return "Eigen";}
};

#ifdef CARIBOU_WITH_MKL
template<typename MatrixType>
struct solver_traits <Eigen::PardisoLU< MatrixType >> {
    static auto TemplateName() -> std::string {return "PardisoLU";}
    static auto BackendName() -> std::string {return "Pardiso";}
};
#endif
}

template<typename EigenSolver>
SparseLUSolver<EigenSolver>::SparseLUSolver()
    : d_backend(initData(&d_backend,
        "backend",
        R"(
            Solver backend used.

            Available backends are:
                Eigen:   Eigen LU solver (SparseLU) [default].
                Pardiso: Pardiso LU solver.
        )",
        true /*displayed_in_GUI*/, true /*read_only_in_GUI*/))
{
    d_backend.setValue(sofa::helper::OptionsGroup(std::vector<std::string> {
        "Eigen", "Pardiso"
    }));

    sofa::helper::WriteAccessor<Data< sofa::helper::OptionsGroup >> backend = d_backend;
    if (internal::solver_traits<EigenSolver>::TemplateName() == "PardisoLU") {
        backend->setSelectedItem(static_cast<unsigned int>(1));
    } else {
        backend->setSelectedItem(static_cast<unsigned int>(0));
    }
}

static int SparseLUSolverClass = sofa::core::RegisterObject("Caribou Sparse LU linear solver")
    .add< SparseLUSolver<Eigen::SparseLU<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor>>> >(true)
#ifdef CARIBOU_WITH_MKL
    .add< SparseLUSolver<Eigen::PardisoLU<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, long long int>>> >()
#endif
;

} // namespace SofaCaribou::solver