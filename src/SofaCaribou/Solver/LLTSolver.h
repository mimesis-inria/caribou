#pragma once

#include <SofaCaribou/Solver/EigenSparseSolver.h>
#include <sofa/helper/OptionsGroup.h>

namespace SofaCaribou::solver {

/**
 * Implementation of a direct LLT sparse linear solver.
 *
 * This class provides a LL^T Cholesky factorizations of sparse matrices that are selfadjoint and positive definite.
 * In order to reduce the fill-in, a symmetric permutation P is applied prior to the factorization such that the
 * factorized matrix is P A P^-1.
 *
 * The component uses the Eigen SimplicialLLT class as the solver backend.
 *
 * @tparam EigenSolver Eigen direct solver type
 */
template <class EigenSolver>
class LLTSolver : public EigenSparseSolver<EigenSolver> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(LLTSolver, EigenSolver), SOFA_TEMPLATE(EigenSparseSolver, EigenSolver));

    template <typename T>
    using Data = sofa::Data<T>;

    LLTSolver();

private:
    ///< Solver backend used (Eigen or Pardiso)
    Data<sofa::helper::OptionsGroup> d_backend;
};

} // namespace SofaCaribou::solver
