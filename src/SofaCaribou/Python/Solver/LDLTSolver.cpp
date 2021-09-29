#include "LDLTSolver.h"


namespace py = pybind11;

namespace SofaCaribou::solver::python {

    void addLDLTSolver(py::module & m) {
        bind_LDLTSolver<Eigen::SimplicialLDLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>, Eigen::Lower, Eigen::AMDOrdering<int>>>(m);

    }

} // namespace SofaCaribou::solver::python