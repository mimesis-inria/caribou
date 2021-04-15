#pragma once

#include <SofaCaribou/config.h>
#include <SofaCaribou/Solver/EigenSolver.h>
#include <sofa/helper/OptionsGroup.h>

namespace SofaCaribou::solver {

template <class EigenSolver_t>
class CARIBOU_API LUSolver : public EigenSolver<typename EigenSolver_t::MatrixType> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(LUSolver, EigenSolver_t), SOFA_TEMPLATE(EigenSolver, typename EigenSolver_t::MatrixType));

    template <typename T>
    using Data = sofa::Data<T>;

    using Base = EigenSolver<typename EigenSolver_t::MatrixType>;
    using Matrix = typename Base::Matrix;
    using Vector = typename Base::Vector;

    LUSolver();

    /** @see LinearSolver::analyze_pattern */
    bool analyze_pattern(const sofa::defaulttype::BaseMatrix * A) override;

    /** @see LinearSolver::factorize */
    bool factorize(const sofa::defaulttype::BaseMatrix * A) override;

    /**
     * @see SofaCaribou::solver::LinearSolver::solve
     */
    bool solve(const sofa::defaulttype::BaseVector * F, sofa::defaulttype::BaseVector * X) const override;

    /**
     * States if the system matrix is symmetric. Note that this value isn't set automatically, the user must
     * explicitly specify it using set_symmetric(true). When it is true, some optimizations will be enabled.
     */
    inline bool symmetric() const override {return d_is_symmetric.getValue();}

    /**
     * Explicitly states if this matrix is symmetric.
     */
    inline void set_symmetric(bool is_symmetric) override { d_is_symmetric.setValue(is_symmetric); }

    // Get the backend name of the class derived from the EigenSolver template parameter
    static std::string BackendName();
private:
    /// Solver backend used (Eigen or Pardiso)
    Data<sofa::helper::OptionsGroup> d_backend;

    /// States if the system matrix is symmetric. This will enable some optimizations.
    Data<bool> d_is_symmetric;

    /// The actual Eigen solver used (its type is passed as a template parameter and must be derived from Eigen::SparseSolverBase)
    EigenSolver_t p_solver;
};

} // namespace SofaCaribou::solver
