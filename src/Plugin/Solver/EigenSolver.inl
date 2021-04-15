#pragma once

#include <SofaCaribou/Solver/EigenSolver.h>

#include <SofaCaribou/Solver/EigenSolver.h>
#include<SofaCaribou/Algebra/EigenMatrix.h>
#include <SofaCaribou/Visitor/AssembleGlobalMatrix.h>
#include <SofaCaribou/Visitor/ConstrainGlobalMatrix.h>

#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <SofaEigen2Solver/EigenVectorWrapper.h>

#ifndef _WIN32
#include <cxxabi.h>
#define CARIBOU_HAS_CXXABI_H
#endif

namespace SofaCaribou::solver {

template <class EigenMatrix_t>
void EigenSolver<EigenMatrix_t>::resetSystem() {
    p_A.resize(0, 0);
    p_x.resize(0);
    p_b.resize(0);
    p_accessor.clear();
}

template <class EigenMatrix_t>
auto EigenSolver<EigenMatrix_t>::assemble (const sofa::core::MechanicalParams* mparams, SofaCaribou::Algebra::EigenMatrix<Matrix> & A) const -> sofa::component::linearsolver::DefaultMultiMatrixAccessor
{
    using Timer = sofa::helper::AdvancedTimer;
    sofa::component::linearsolver::DefaultMultiMatrixAccessor accessor;

    // Step 1. Preparation stage
    //         This stage go down on the sub-graph and gather the top-level mechanical objects (mechanical objects that
    //         aren't slaved of another mechanical object through a mechanical mapping). The matrix isn't built yet,
    //         we only keep a list of mechanical objects and mappings, and get the size of global system matrix that
    //         will be built in the next stage. Note, however, that the mapping matrices are accumulated (calling
    //         BaseMapping::getJ() on every mapping).
    Timer::stepBegin("PrepareMatrix");

    auto context = const_cast<sofa::core::objectmodel::BaseContext *> (this->getContext());
    sofa::simulation::common::MechanicalOperations mops(mparams, context);

    // Step 1.1 Get dimension of each top level mechanical states using BaseMechanicalState::getMatrixSize(),
    //          and accumulate mechanical objects and mapping matrices
    Timer::stepBegin("Dimension");
    mops.getMatrixDimension(nullptr, nullptr, &accessor);
    const auto n = static_cast<Eigen::Index>(accessor.getGlobalDimension());
    Timer::stepEnd("Dimension");

    // Step 2.2 Does nothing more than to accumulate from the previous step a list of
    //          "MatrixRef = <MechanicalState*, MatrixIndex>" where MatrixIndex is the
    //          (i,i) position of the given top level MechanicalState* inside the global
    //          system matrix. This global matrix hence contains one sub-matrix per top
    //          level mechanical state.
    Timer::stepBegin("SetupMatrixIndices");
    accessor.setupMatrices();
    Timer::stepEnd("SetupMatrixIndices");

    Timer::stepBegin("Clear");
    A.resize(n, n);
    A.set_symmetric(symmetric()); // Enables some optimization when the system matrix is symmetric
    accessor.setGlobalMatrix(&A);
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
    execute(visitor::AssembleGlobalMatrix(mparams, &accessor));
    Timer::stepEnd("AssembleGlobalMatrix");

    Timer::stepBegin("ConstrainGlobalMatrix");
    execute(visitor::ConstrainGlobalMatrix(mparams, &accessor));
    Timer::stepEnd("ConstrainGlobalMatrix");

    // Step 3. Mechanical mappings
    //         In case we have mapped matrices, which is, system matrix of a slave mechanical object, accumulate its
    //         contribution to the global system matrix with:
    //           [A]ij += Jt * [A']ij * J
    //         where A is the master mechanical object's matrix, A' is the slave mechanical object matrix and J=m.getJ()
    //         is the mapping relation between the slave and its master.
    Timer::stepBegin("MappedMatrices");
    accessor.computeGlobalMatrix();
    Timer::stepEnd("MappedMatrices");

    // Step 4. Convert the system matrix to a compressed sparse matrix
    Timer::stepBegin("ConvertToSparse");
    A.compress();
    Timer::stepEnd("ConvertToSparse");
    Timer::stepEnd("BuildMatrix");

    return accessor;
}

template <class EigenMatrix_t>
void EigenSolver<EigenMatrix_t>::setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams) {
    using Timer = sofa::helper::AdvancedTimer;
    Timer::stepBegin("EigenSolver::AssembleGlobalMatrix");

    // This will be set back to true once the system matrix has been successfully factorized
    // It is used to avoid trying to solve the system when something went wrong in this stage.
    p_A_is_factorized = false;

    // Save the current mechanical parameters (m, b and k factors of the mass (M), damping (B) and
    // stiffness (K) matrices)
    p_mechanical_params = *mparams;

    // Step 1. Assemble the system matrix
    auto previous_dimension = p_A.rowSize();
    p_accessor = assemble(mparams, p_A);
    auto current_dimension = p_A.rowSize();

    // Step 2. Let the solver analyse the matrix
    bool matrix_shape_has_changed = previous_dimension != current_dimension;
    if (matrix_shape_has_changed) {
        Timer::stepBegin("MatrixAnalysis");
        if (not this->analyze_pattern(&p_A)) {
            msg_error() << "Failed to analyse the system matrix pattern";
        }
        Timer::stepEnd("MatrixAnalysis");

    }

    // Step 5. Factorize the system matrix
    Timer::stepBegin("MatrixFactorization");
    p_A_is_factorized = this->factorize(&p_A);
    if (not p_A_is_factorized) {
        msg_error() << "Failed to factorize the system matrix";
    }
    Timer::stepEnd("MatrixFactorization");

    Timer::stepEnd("EigenSolver::AssembleGlobalMatrix");
}

template <class EigenMatrix_t>
void EigenSolver<EigenMatrix_t>::setSystemRHVector(sofa::core::MultiVecDerivId b_id) {
    using Timer = sofa::helper::AdvancedTimer;
    Timer::stepBegin("EigenSolver::AssembleResidualVector");

    // Register the RHS to the mechanical parameters
    p_mechanical_params.setDf(b_id);

    sofa::simulation::common::MechanicalOperations mop(&p_mechanical_params, this->getContext());
    p_b_id = b_id;

    // Copy the vectors of the mechanical objects into a global eigen vector.
    p_b.resize(p_A.rowSize());
    mop.multiVector2BaseVector(p_b_id, &p_b, &p_accessor);

    Timer::stepEnd("EigenSolver::AssembleResidualVector");
}

template <class EigenMatrix_t>
void EigenSolver<EigenMatrix_t>::setSystemLHVector(sofa::core::MultiVecDerivId x_id) {
    using Timer = sofa::helper::AdvancedTimer;
    Timer::stepBegin("EigenSolver::AssembleSolutionVector");

    // Register the dx vector to the mechanical parameters
    p_mechanical_params.setDx(x_id);

    sofa::simulation::common::MechanicalOperations mop(&p_mechanical_params, this->getContext());
    p_x_id = x_id;


    // Copy the vectors of the mechanical objects into a global eigen vector.
    p_x.resize(p_A.rowSize());
    mop.multiVector2BaseVector(p_x_id, &p_x, &p_accessor);

    Timer::stepEnd("EigenSolver::AssembleSolutionVector");
}

template <class EigenMatrix_t>
void EigenSolver<EigenMatrix_t>::solveSystem() {
    using Timer = sofa::helper::AdvancedTimer;
    using namespace sofa::component::linearsolver;
    using namespace sofa::core::behavior;

    sofa::simulation::common::MechanicalOperations mop( &p_mechanical_params, this->getContext() );

    Timer::stepBegin("EigenSolver::solve");
    bool success = this->solve(&p_b, &p_x);
    if (success) {
        // Copy the solution into the mechanical objects of the current context sub-graph.
        mop.baseVector2MultiVector(&p_x, p_x_id, &p_accessor);
    }

    Timer::stepEnd("EigenSolver::solve");
}

template<typename EigenMatrix_t>
std::string EigenSolver<EigenMatrix_t>::GetCustomTemplateName() {
    std::string namestring;

#if defined( CARIBOU_HAS_CXXABI_H )
    int status;
    char* name = abi::__cxa_demangle(typeid(EigenSolver).name(), 0, 0, &status);
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
        throw std::runtime_error("Demangling of '"+std::string(typeid(EigenMatrix_t).name())+"' failed for the following reason: " + error);
    }
#else
    const char* name = typeid(EigenSolver).name();
#endif

    return namestring;
}

template <class EigenMatrix_t>
template<typename Derived>
auto EigenSolver<EigenMatrix_t>::canCreate(Derived*, sofa::core::objectmodel::BaseContext*, sofa::core::objectmodel::BaseObjectDescription* arg) -> bool {
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

} // namespace SofaCaribou::solver