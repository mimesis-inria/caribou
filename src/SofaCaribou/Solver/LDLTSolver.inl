#pragma once

#include <SofaCaribou/Algebra/EigenMatrix.h>
#include <SofaCaribou/Solver/LDLTSolver.h>
#include <SofaCaribou/Solver/EigenSolver.inl>

#include<Eigen/SparseCholesky>

#include <algorithm>
#include <cctype>
#include <string>

#ifdef CARIBOU_WITH_MKL

// Bug introduced in Eigen 3.3.8, fixed in bfdd4a9
#ifndef EIGEN_USING_STD
#define EIGEN_USING_STD(a) EIGEN_USING_STD_MATH(a)
#endif

#include <Eigen/PardisoSupport>
#endif

namespace SofaCaribou::solver {

namespace {
template <typename T>
struct solver_traits {};

template<typename MatrixType, int UpLo, typename Ordering>
struct solver_traits<Eigen::SimplicialLDLT < MatrixType, UpLo, Ordering>> {
static auto BackendName() -> std::string { return "Eigen"; }
};

#ifdef CARIBOU_WITH_MKL
template<typename MatrixType, int UpLo>
struct solver_traits <Eigen::PardisoLDLT< MatrixType, UpLo >> {
    static auto BackendName() -> std::string {return "Pardiso";}
};
#endif
}

template<class EigenSolver_t>
bool LDLTSolver<EigenSolver_t>::analyze_pattern() {
    auto A_ = this->A();
    if (not A_) {
        throw std::runtime_error("Tried to analyze an incompatible matrix (not an Eigen matrix).");
    }

    p_solver.analyzePattern(A_->matrix());

    return (p_solver.info() == Eigen::Success);
}

template<class EigenSolver_t>
bool LDLTSolver<EigenSolver_t>::factorize() {
    auto A_ = this->A();
    if (not A_) {
        throw std::runtime_error("Tried to analyze an incompatible matrix (not an Eigen matrix).");
    }

    p_solver.factorize(A_->matrix());

    return (p_solver.info() == Eigen::Success);
}

template<class EigenSolver_t>
bool LDLTSolver<EigenSolver_t>::solve(const SofaCaribou::Algebra::BaseVector * F,
                                      SofaCaribou::Algebra::BaseVector *X) {
    auto F_ = dynamic_cast<const SofaCaribou::Algebra::EigenVector<Vector> *>(F);
    auto X_ = dynamic_cast<SofaCaribou::Algebra::EigenVector<Vector> *>(X);

    X_->vector() = p_solver.solve(F_->vector());
    return (p_solver.info() == Eigen::Success);
}

template<class EigenSolver_t>
std::string LDLTSolver<EigenSolver_t>::BackendName() {
    return solver_traits<EigenSolver_t>::BackendName();
}

template<typename EigenSolver_t>
LDLTSolver<EigenSolver_t>::LDLTSolver()
: d_backend(initData(&d_backend
, "backend"
, R"(
    Solver backend used.

    Available backends are:
    Eigen:   Eigen LDLT solver (SimplicialLDLT) [default].
    Pardiso: Pardiso LDLT solver.
  )" , true /*displayed_in_GUI*/, true /*read_only_in_GUI*/))
{
    d_backend.setValue(sofa::helper::OptionsGroup(std::vector < std::string > {
            "Eigen", "Pardiso"
    }));


    // Put the backend name in lower case
    std::string backend_str = solver_traits<EigenSolver_t>::BackendName();
    std::transform(backend_str.begin(), backend_str.end(), backend_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    sofa::helper::WriteAccessor <Data<sofa::helper::OptionsGroup >> backend = d_backend;
    if (backend_str == "pardiso") { // Case insensitive
        backend->setSelectedItem(static_cast<unsigned int>(1));
    } else {
        backend->setSelectedItem(static_cast<unsigned int>(0));
    }

    // Explicitly state that the matrix is symmetric (would not be possible to do an LDLT decomposition otherwise)
    this->set_symmetric(true);
}


} // namespace SofaCaribou::solver
