#include "NewtonRaphsonSolver.h"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#define NDEBUG
#include <SofaCaribou/Ode/NewtonRaphsonSolver.h>
#include <SofaCaribou/Algebra/EigenMatrix.h>
#include <SofaCaribou/Algebra/EigenVector.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>


namespace py = pybind11;

namespace SofaCaribou::ode::python {

void addNewtonRaphsonSolver(pybind11::module &m) {
    using namespace sofa::core::objectmodel;
    py::class_<NewtonRaphsonSolver, BaseObject, sofapython3::py_shared_ptr<NewtonRaphsonSolver>> c(m, "NewtonRaphsonSolver");

    py::enum_<NewtonRaphsonSolver::PatternAnalysisStrategy>(
            c, "PatternAnalysisStrategy",
            "Different strategies to determine when the pattern of the system matrix should be "
            "analyzed in order to, for example, compute a permutation matrix before factorizing it.")
    .value("ALWAYS", NewtonRaphsonSolver::PatternAnalysisStrategy::ALWAYS, "Analyze the pattern of the matrix at each Newton iterations.")
    .value("NEVER", NewtonRaphsonSolver::PatternAnalysisStrategy::NEVER, "Never analyze the pattern of the matrix.")
    .value("BEGINNING_OF_THE_SIMULATION", NewtonRaphsonSolver::PatternAnalysisStrategy::BEGINNING_OF_THE_SIMULATION, "Analyze the pattern of the matrix once at the beginning of the simulation.")
    .value("BEGINNING_OF_THE_TIME_STEP", NewtonRaphsonSolver::PatternAnalysisStrategy::BEGINNING_OF_THE_TIME_STEP, "Analyze the pattern of the matrix once per time step.")
    .export_values();

    py::enum_<NewtonRaphsonSolver::Event>(
            c, "Event",
            "Events associated to registered callback function for different time in the NR iteration process.")
    .value("ITERATION_BEGIN", NewtonRaphsonSolver::Event::ITERATION_BEGIN, "Event triggered at the very beginning of the Newton iteration.")
    .value("MATRIX_ASSEMBLED", NewtonRaphsonSolver::Event::MATRIX_ASSEMBLED, "Event triggered after the system matrix has been assembled.")
    .value("MATRIX_ANALYZED", NewtonRaphsonSolver::Event::MATRIX_ANALYZED, "Event triggered after the system matrix has been analyzed.")
    .value("MATRIX_FACTORIZED", NewtonRaphsonSolver::Event::MATRIX_FACTORIZED, "Event triggered after the system matrix has been factorized.")
    .value("INCREMENT_SOLVED", NewtonRaphsonSolver::Event::INCREMENT_SOLVED, "Event triggered after the system solution vector has been solved.")
    .value("INCREMENT_PROPAGATED", NewtonRaphsonSolver::Event::INCREMENT_PROPAGATED, "Event triggered after the system solution vector has been propagated through mappings.")
    .value("RESIDUAL_UPDATED", NewtonRaphsonSolver::Event::RESIDUAL_UPDATED, "Event triggered after the system force vector has been updated.")
    .value("ITERATION_END", NewtonRaphsonSolver::Event::ITERATION_END, "Event triggered at the very end of the Newton iteration.")
    .export_values();

    c.def_property_readonly("iteration_times", &NewtonRaphsonSolver::iteration_times, "List of times (in nanoseconds) that each Newton-Raphson iteration took to compute in the last call to Solve().");
    c.def_property_readonly("squared_residuals", &NewtonRaphsonSolver::squared_residuals, "The list of squared residual norms (||r||^2) of every newton iterations of the last solve call.");
    c.def_property_readonly("squared_initial_residual", &NewtonRaphsonSolver::squared_initial_residual, "The initial squared residual (||r0||^2) of the last solve call.");
    c.def_property("pattern_analysis_strategy", &NewtonRaphsonSolver::pattern_analysis_strategy, &NewtonRaphsonSolver::set_pattern_analysis_strategy, "Get the current strategy that determine when the pattern of the system matrix should be analyzed.");
    c.def_property_readonly("current_iteration", &NewtonRaphsonSolver::current_iteration, "Get the current newton iteration within the time step solve call. First iteration is 1 (index starts from 1).");

    c.def_property_readonly("A", [](const NewtonRaphsonSolver & solver) -> py::object {
        if (solver.A() == nullptr) {
            return py::none();
        }

        if (const auto * a = dynamic_cast<const SofaCaribou::Algebra::EigenMatrix<Eigen::SparseMatrix<float, Eigen::ColMajor, int>> *>(solver.A())) {
            return py::cast(a->matrix());
        }

        if (const auto * a = dynamic_cast<const SofaCaribou::Algebra::EigenMatrix<Eigen::SparseMatrix<double, Eigen::ColMajor, int>> *>(solver.A())) {
            return py::cast(a->matrix());
        }

        if (const auto * a = dynamic_cast<const SofaCaribou::Algebra::EigenMatrix<Eigen::SparseMatrix<float, Eigen::RowMajor, int>> *>(solver.A())) {
            return py::cast(a->matrix());
        }

        if (const auto * a = dynamic_cast<const SofaCaribou::Algebra::EigenMatrix<Eigen::SparseMatrix<double, Eigen::RowMajor, int>> *>(solver.A())) {
            return py::cast(a->matrix());
        }

        throw std::runtime_error("Cannot bind the system matrix type to python.");
    }, "Complete matrix of the linearized part of system");

    c.def_property_readonly("dx", [](const NewtonRaphsonSolver & solver) -> py::object {
        if (solver.dx() == nullptr) {
            return py::none();
        }

        if (const auto * a = dynamic_cast<const SofaCaribou::Algebra::EigenVector<Eigen::Matrix<double, Eigen::Dynamic, 1>> *>(solver.dx())) {
            return py::cast(a->vector());
        }

        if (const auto * a = dynamic_cast<const SofaCaribou::Algebra::EigenVector<Eigen::Matrix<float, Eigen::Dynamic, 1>> *>(solver.dx())) {
            return py::cast(a->vector());
        }

        throw std::runtime_error("Cannot bind the system increment vector type to python.");
        }, "Solution of the linearized system");

    c.def_property_readonly("F", [](const NewtonRaphsonSolver & solver) -> py::object {
        if (solver.F() == nullptr) {
            return py::none();
        }

        if (const auto * a = dynamic_cast<const SofaCaribou::Algebra::EigenVector<Eigen::Matrix<double, Eigen::Dynamic, 1>> *>(solver.F())) {
            return py::cast(a->vector());
        }

        if (const auto * a = dynamic_cast<const SofaCaribou::Algebra::EigenVector<Eigen::Matrix<float, Eigen::Dynamic, 1>> *>(solver.F())) {
            return py::cast(a->vector());
        }

        throw std::runtime_error("Cannot bind the system RHS vector type to python.");
        }, "Force vector (right-hand side) of the system");

    c.def("register_callback", [](NewtonRaphsonSolver & solver, const NewtonRaphsonSolver::Event & e, py::function callback) {
        solver.register_callback(e, [callback](const NewtonRaphsonSolver & s) {
            callback(py::cast(s, py::return_value_policy::reference));
        });
    }, "Register a callback function to be called at the specific NR event (see NewtonRaphsonSolver::Event for the list of events).");
}
}