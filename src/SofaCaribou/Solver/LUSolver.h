#pragma once

#include <SofaCaribou/Solver/EigenSparseSolver.h>
#include <sofa/helper/OptionsGroup.h>

namespace SofaCaribou::solver {

template <class EigenSolver>
class LUSolver : public EigenSparseSolver<EigenSolver> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(LUSolver, EigenSolver), SOFA_TEMPLATE(EigenSparseSolver, EigenSolver));

    template <typename T>
    using Data = sofa::Data<T>;

    LUSolver();

    /**
     * States if the system matrix is symmetric. Note that this value isn't set automatically, the user must
     * explicitly specify it using set_symmetric(true). When it is true, some optimizations will be enabled.
     */
    inline bool symmetric() const override {return d_is_symmetric.getValue();}

    /**
     * Explicitly states if this matrix is symmetric.
     */
    inline void set_symmetric(bool is_symmetric) override { d_is_symmetric.setValue(is_symmetric); }

    /**
     * Assemble the system matrix A = (mM + bB + kK) inside the SparseMatrix p_A.
     * @param mparams Mechanical parameters containing the m, b and k factors.
     */
    void assemble (const sofa::core::MechanicalParams* mparams) override;
private:
    ///< Solver backend used (Eigen or Pardiso)
    Data<sofa::helper::OptionsGroup> d_backend;

    ///< States if the system matrix is symmetric. This will enable some optimizations.
    Data<bool> d_is_symmetric;
};

} // namespace SofaCaribou::solver
