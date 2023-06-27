#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Solver/EigenSolver.h>

DISABLE_ALL_WARNINGS_BEGIN
#include <sofa/helper/OptionsGroup.h>
DISABLE_ALL_WARNINGS_END

namespace SofaCaribou::solver {

/**
 * Implementation of a direct LDLT sparse linear solver.
 *
 * This class provides a LDL^T Cholesky factorizations of sparse matrices that are selfadjoint and positive definite.
 * In order to reduce the fill-in, a symmetric permutation P is applied prior to the factorization such that the
 * factorized matrix is P A P^-1.
 *
 * The component uses the Eigen SimplicialLDLT class as the solver backend.
 *
 * @tparam EigenSolver_t
 */
template <class EigenSolver_t>
class LDLTSolver : public EigenSolver<typename EigenSolver_t::MatrixType> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(LDLTSolver, EigenSolver_t), SOFA_TEMPLATE(EigenSolver, typename EigenSolver_t::MatrixType));

    template <typename T>
    using Data = sofa::Data<T>;

    using Base = EigenSolver<typename EigenSolver_t::MatrixType>;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

    
    LDLTSolver();

    /** @see LinearSolver::analyze_pattern */
    
    bool analyze_pattern() override;

    /** @see LinearSolver::factorize */
    
    bool factorize() override;

    /**
     * @see SofaCaribou::solver::LinearSolver::solve
     */
    
    bool solve(const SofaCaribou::Algebra::BaseVector * F, SofaCaribou::Algebra::BaseVector * X) override;

    // Get the backend name of the class derived from the EigenSolver template parameter
    
    static std::string BackendName();
private:
    /// Solver backend used (Eigen or Pardiso)
    Data<sofa::helper::OptionsGroup> d_backend;

    /// The actual Eigen solver used (its type is passed as a template parameter and must be derived from Eigen::SparseSolverBase)
    EigenSolver_t p_solver;
};

} // namespace SofaCaribou::solver
