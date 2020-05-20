#pragma once

#include <SofaCaribou/Solver/EigenSparseSolver.h>
#include <sofa/helper/OptionsGroup.h>

namespace SofaCaribou::solver {

template <class EigenSolver>
class SparseLUSolver : public EigenSparseSolver<EigenSolver> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(SparseLUSolver, EigenSolver), SOFA_TEMPLATE(EigenSparseSolver, EigenSolver));

    template <typename T>
    using Data = sofa::Data<T>;

    SparseLUSolver();
private:
    ///< Solver backend used (Eigen or Pardiso)
    Data<sofa::helper::OptionsGroup> d_backend;
};

} // namespace SofaCaribou::solver
