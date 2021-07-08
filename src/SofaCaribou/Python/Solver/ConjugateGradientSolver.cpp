#include "ConjugateGradientSolver.h"

#include <SofaCaribou/Solver/ConjugateGradientSolver.h>

namespace py = pybind11;

namespace SofaCaribou::solver::python {

void addConjugateGradientSolver(py::module & m) {
    bind_ConjugateGradientSolver<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::RowMajor, int>>(m);
}

} // namespace SofaCaribou::solver::python