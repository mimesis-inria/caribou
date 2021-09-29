#include "LUSolver.h"

namespace py = pybind11;

namespace SofaCaribou::solver::python {

    void addLUSolver(py::module & m) {
        bind_LUSolver<Eigen::SparseLU<Eigen::SparseMatrix<FLOATING_POINT_TYPE, Eigen::ColMajor, int>, Eigen::AMDOrdering<int>>>(m);

    }

} // namespace SofaCaribou::solver::python