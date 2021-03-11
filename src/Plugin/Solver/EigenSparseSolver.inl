#pragma once

#include <SofaCaribou/Solver/EigenSparseSolver.h>

#include <SofaCaribou/Solver/EigenSparseSolver.h>
#include<SofaCaribou/Algebra/EigenMatrix.h>
#include <SofaCaribou/Visitor/AssembleGlobalMatrix.h>
#include <SofaCaribou/Visitor/ConstrainGlobalMatrix.h>

#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <SofaEigen2Solver/EigenVectorWrapper.h>

#include <cxxabi.h>

namespace SofaCaribou::solver {

template <class EigenSolver>
void EigenSparseSolver<EigenSolver>::resetSystem() {
    p_A.resize(0, 0);
    p_x.resize(0);
    p_b.resize(0);
    p_accessor.clear();
}

template <class EigenSolver>
void EigenSparseSolver<EigenSolver>::assemble (const sofa::core::MechanicalParams* mparams) {
    using Timer = sofa::helper::AdvancedTimer;

    // Step 1. Preparation stage
    //         This stage go down on the sub-graph and gather the top-level mechanical objects (mechanical objects that
    //         aren't slaved of another mechanical object through a mechanical mapping). The matrix isn't built yet,
    //         we only keep a list of mechanical objects and mappings, and get the size of global system matrix that
    //         will be built in the next stage. Note, however, that the mapping matrices are accumulated (calling
    //         BaseMapping::getJ() on every mapping).
    Timer::stepBegin("PrepareMatrix");

    sofa::simulation::common::MechanicalOperations mops(mparams, this->getContext());
    p_accessor.clear();

    // Step 1.1 Get dimension of each top level mechanical states using BaseMechanicalState::getMatrixSize(),
    //          and accumulate mechanical objects and mapping matrices
    Timer::stepBegin("Dimension");
    mops.getMatrixDimension(nullptr, nullptr, &p_accessor);
    const auto n = static_cast<Eigen::Index>(p_accessor.getGlobalDimension());
    Timer::stepEnd("Dimension");

    // Step 2.2 Does nothing more than to accumulate from the previous step a list of
    //          "MatrixRef = <MechanicalState*, MatrixIndex>" where MatrixIndex is the
    //          (i,i) position of the given top level MechanicalState* inside the global
    //          system matrix. This global matrix hence contains one sub-matrix per top
    //          level mechanical state.
    Timer::stepBegin("SetupMatrixIndices");
    p_accessor.setupMatrices();
    Timer::stepEnd("SetupMatrixIndices");

    Timer::stepBegin("Clear");
    p_A.resize(n, n);
    Algebra::EigenMatrix<SparseMatrix &> wrapper (p_A);
    wrapper.set_symmetric(symmetric()); // Enables some optimization when the system matrix is symmetric
    p_accessor.setGlobalMatrix(&wrapper);
    Timer::stepEnd("Clear");

    Timer::stepEnd("PrepareMatrix");

    Timer::stepBegin("BuildMatrix");

    // Step 2. Building stage
    //         Here we go down on the current context sub-graph and call :
    //           1. ff->addMBKToMatrix(&K) for every force field "ff" found.
    //           2. pc->applyConstraint(&K) for every BaseProjectiveConstraintSet "pc" found.
    //         If a mechanical mapping "m" is found during the traversal, and m->areMatricesMapped() is false, the
    //         traversal stops in the subgraph of the mapping.

    sofa::simulation::common::VisitorExecuteFunc execute(*mops.ctx);
    Timer::stepBegin("AssembleGlobalMatrix");
    execute(visitor::AssembleGlobalMatrix(mparams, &p_accessor));
    Timer::stepEnd("AssembleGlobalMatrix");

    Timer::stepBegin("ConstrainGlobalMatrix");
    execute(visitor::ConstrainGlobalMatrix(mparams, &p_accessor));
    Timer::stepEnd("ConstrainGlobalMatrix");

    // Step 3. Mechanical mappings
    //         In case we have mapped matrices, which is, system matrix of a slave mechanical object, accumulate its
    //         contribution to the global system matrix with:
    //           [A]ij += Jt * [A']ij * J
    //         where A is the master mechanical object's matrix, A' is the slave mechanical object matrix and J=m.getJ()
    //         is the mapping relation between the slave and its master.
    Timer::stepBegin("MappedMatrices");
    p_accessor.computeGlobalMatrix();
    Timer::stepEnd("MappedMatrices");

    // Step 4. Convert the system matrix to a compressed sparse matrix
    Timer::stepBegin("ConvertToSparse");
    wrapper.compress();
    Timer::stepEnd("ConvertToSparse");
    Timer::stepEnd("BuildMatrix");

    // Remove the global matrix from the accessor since the wrapper was temporary
    p_accessor.setGlobalMatrix(nullptr);
}

template <class EigenSolver>
void EigenSparseSolver<EigenSolver>::setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) {
    using Timer = sofa::helper::AdvancedTimer;
    Timer::stepBegin("EigenSparseSolver::AssembleGlobalMatrix");

    // This will be set back to true once the system matrix has been successfully factorized
    // It is used to avoid trying to solve the system when something went wrong in this stage.
    p_A_is_factorized = false;

    // Save the current mechanical parameters (m, b and k factors of the mass (M), damping (B) and
    // stiffness (K) matrices)
    p_mechanical_params = *mparams;

    // Step 1. Assemble the system matrix
    auto previous_dimension = p_A.rows();
    assemble(mparams);
    auto current_dimension = p_A.rows();

    // Step 2. Let the solver analyse the matrix
    bool matrix_shape_has_changed = previous_dimension != current_dimension;
    if (matrix_shape_has_changed) {
        Timer::stepBegin("MatrixAnalysis");
        p_solver.analyzePattern(p_A);
        Timer::stepEnd("MatrixAnalysis");
        if (p_solver.info() != Eigen::Success) {
            msg_error() << "Failed to analyse the system matrix pattern (" +
            std::string((p_solver.info() == Eigen::NumericalIssue) ? "numerical issues" : "invalid input") + ")";
        }
    }

    // Step 5. Factorize the system matrix
    Timer::stepBegin("MatrixFactorization");
    p_solver.factorize(p_A);
    if (p_solver.info() != Eigen::Success) {
        msg_error() << "Failed to factorize the system matrix (" +
                       std::string((p_solver.info() == Eigen::NumericalIssue) ? "numerical issues" : "invalid input") + ")";
    } else {
        p_A_is_factorized = true;
    }
    Timer::stepEnd("MatrixFactorization");

    Timer::stepEnd("EigenSparseSolver::AssembleGlobalMatrix");
}

template <class EigenSolver>
void EigenSparseSolver<EigenSolver>::setSystemRHVector(sofa::core::MultiVecDerivId b_id) {
    using Timer = sofa::helper::AdvancedTimer;
    Timer::stepBegin("EigenSparseSolver::AssembleResidualVector");

    // Register the RHS to the mechanical parameters
    p_mechanical_params.setDf(b_id);

    sofa::simulation::common::MechanicalOperations mop(&p_mechanical_params, this->getContext());
    p_b_id = b_id;

    // Copy the vectors of the mechanical objects into a global eigen vector.
    p_b.resize(p_A.rows());
    sofa::component::linearsolver::EigenVectorWrapper<FLOATING_POINT_TYPE> b(p_b);
    mop.multiVector2BaseVector(p_b_id, &b, &p_accessor);

    Timer::stepEnd("EigenSparseSolver::AssembleResidualVector");
}

template <class EigenSolver>
void EigenSparseSolver<EigenSolver>::setSystemLHVector(sofa::core::MultiVecDerivId x_id) {
    using Timer = sofa::helper::AdvancedTimer;
    Timer::stepBegin("EigenSparseSolver::AssembleSolutionVector");

    // Register the dx vector to the mechanical parameters
    p_mechanical_params.setDx(x_id);

    sofa::simulation::common::MechanicalOperations mop(&p_mechanical_params, this->getContext());
    p_x_id = x_id;


    // Copy the vectors of the mechanical objects into a global eigen vector.
    p_x.resize(p_A.rows());
    sofa::component::linearsolver::EigenVectorWrapper<FLOATING_POINT_TYPE> x(p_x);
    mop.multiVector2BaseVector(p_x_id, &x, &p_accessor);

    Timer::stepEnd("EigenSparseSolver::AssembleSolutionVector");
}

template <class EigenSolver>
void EigenSparseSolver<EigenSolver>::solveSystem() {
    using Timer = sofa::helper::AdvancedTimer;
    using namespace sofa::component::linearsolver;
    using namespace sofa::core::behavior;

    sofa::simulation::common::MechanicalOperations mop( &p_mechanical_params, this->getContext() );

    Timer::stepBegin("EigenSparseSolver::solve");

    if (p_A_is_factorized) {
        // Solve the system
        p_x = p_solver.solve(p_b);
        if (p_solver.info() != Eigen::Success) {
            msg_error() << "Failed to solve the system (" +
                           std::string(
                               (p_solver.info() == Eigen::NumericalIssue) ? "numerical issues" : "invalid input") + ")";
        }

        // Copy the solution into the mechanical objects of the current context sub-graph.
        EigenVectorWrapper<FLOATING_POINT_TYPE> x_wrapper(p_x);
        mop.baseVector2MultiVector(&x_wrapper, p_x_id, &p_accessor);
    }

    Timer::stepEnd("EigenSparseSolver::solve");
}

template<typename EigenSolver>
std::string EigenSparseSolver<EigenSolver>::GetCustomTemplateName() {
    int status;
    char * name = abi::__cxa_demangle(typeid(EigenSolver).name(), 0, 0, &status);
    std::string namestring;
    std::string error;
    switch (status) {
        case 0: namestring = std::string(name); break;
        case -1: error = "A memory allocation failure occurred."; break;
        case -2: error = "The mangled name is not a valid name under the C++ ABI mangling rules."; break;
        case -3: error = "One of the arguments passed to the demangling is invalid."; break;
        default: error = "Unknown error."; break;
    }
    free(name);
    if (!error.empty()) {
        throw std::runtime_error("Demangling of '"+std::string(typeid(EigenSolver).name())+"' failed for the following reason: " + error);
    }
    return namestring;
}

template <class EigenSolver>
template<typename Derived>
auto EigenSparseSolver<EigenSolver>::canCreate(Derived*, sofa::core::objectmodel::BaseContext*, sofa::core::objectmodel::BaseObjectDescription* arg) -> bool {
    std::string requested_backend = arg->getAttribute( "backend", "");
    std::string current_backend = Derived::BackendName();

    // Let's check if the user has specified a backend
    if (not requested_backend.empty()) {
        if (requested_backend == current_backend) {
            return true;
        } else {
            arg->logError("The requested backend '" + requested_backend +
                          "' isn't compatible with the following template parameters: '" + GetCustomTemplateName() + "'.");
            return false;
        }
    }

    // No backend specified
    arg->setAttribute("backend", Derived::BackendName());
    return true;
}

template<class EigenSolver>
bool
EigenSparseSolver<EigenSolver>::solve(const sofa::defaulttype::BaseMatrix *A,
                                      const sofa::defaulttype::BaseVector *F,
                                      sofa::defaulttype::BaseVector *X) const
{
    using Timer = sofa::helper::AdvancedTimer;
    using namespace SofaCaribou::Algebra;

    auto eigen_A = dynamic_cast<const EigenMatrix<SparseMatrix> *>(A);
    auto eigen_F = dynamic_cast<const EigenVector<Vector> *>(F);
    auto eigen_X = dynamic_cast<EigenVector<Vector> *>(X);

    // Step 1. Let the solver analyse the matrix
    EigenSolver solver;
    Timer::stepBegin("MatrixAnalysis");
    solver.analyzePattern(eigen_A->matrix());
    Timer::stepEnd("MatrixAnalysis");
    if (solver.info() != Eigen::Success) {
        msg_error() << "Failed to analyse the system matrix pattern (" +
                       std::string((solver.info() == Eigen::NumericalIssue) ? "numerical issues" : "invalid input") + ")";

        return false;
    }

    // Step 2. Factorize the system matrix
    Timer::stepBegin("MatrixFactorization");
    solver.factorize(eigen_A->matrix());
    if (solver.info() != Eigen::Success) {
        msg_error() << "Failed to factorize the system matrix (" +
                       std::string((solver.info() == Eigen::NumericalIssue) ? "numerical issues" : "invalid input") + ")";
        return false;
    }
    Timer::stepEnd("MatrixFactorization");

    // Step3. Solve the system
    Timer::stepBegin("Solve");
    eigen_X->vector() = solver.solve(eigen_F->vector());
    Timer::stepEnd("Solve");
    if (solver.info() != Eigen::Success) {
        msg_error() << "Failed to solve the system (" +
                       std::string( (solver.info() == Eigen::NumericalIssue) ? "numerical issues" : "invalid input") + ")";
        return false;
    }

    return true;
}

} // namespace SofaCaribou::solver