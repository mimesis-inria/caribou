#include "LLTSolver.h"

namespace py = pybind11;

namespace SofaCaribou::solver::python {

    void addLLTSolver(py::module & m) {
        bind_LLTSolver<Eigen::SimplicialLLT<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>, Eigen::Lower, Eigen::AMDOrdering<int>>>(m);
    }

} // namespace SofaCaribou::solver::python