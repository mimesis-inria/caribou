/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2016 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
// Author: Hadrien Courtecuisse
//
// Copyright: See COPYING file that comes with this distribution
#include <fstream>
#include <iostream>

#include "PardisoSolver.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <cmath>
#include <sofa/helper/system/thread/CTime.h>
#include <sofa/helper/AdvancedTimer.h>
#include <SofaBaseLinearSolver/CompressedRowSparseMatrix.inl>
#include <sofa/simulation/AnimateBeginEvent.h>

#include <sys/types.h>
#include <sys/stat.h>

#ifndef WIN32
#include <unistd.h>
#else
#include <windows.h>
#endif
#include <cstdlib>
#include "PardisoSolver.h"

/* Change this if your Fortran compiler does not append underscores. */
/* e.g. the AIX compiler:  #define F77_FUNC(func) func                */

#ifdef AIX
#define F77_FUNC(func)  func
#else
#define F77_FUNC(func)  func ## _
#endif
extern "C" {

/* PARDISO prototype. */
extern  int F77_FUNC(pardisoinit)
    (void *, int *, int *, int *, double *, int *);

extern  int F77_FUNC(pardiso)
    (void *, int *, int *, int *, int *, int *,
     double *, int *, int *, int *, int *, int *,
     int *, double *, double *, int *, double *);

extern int F77_FUNC(pardiso_get_schur)
    (void *, int*, int*, int*, double*, int*, int*);

} // "C"


namespace sofa
{

namespace component
{

namespace linearsolver
{

template<class TMatrix, class TVector>
SparsePARDISOSolver<TMatrix,TVector>::SparsePARDISOSolverInvertData::SparsePARDISOSolverInvertData(int d_symmetric)
    : solver(nullptr)
    , pardiso_initerr(1)
    , pardiso_mtype(0)
    , factorized(false)
{
    factorized = false;
    pardiso_initerr = 0;

    msg_info("SparsePARDISOSolverInvertData") << "FSYM: " << d_symmetric;

    switch(d_symmetric)
    {
        case  0: pardiso_mtype = 11; break; // real and nonsymmetric
        case  1: pardiso_mtype = -2; break; // real and symmetric indefinite
        case  2: pardiso_mtype =  2; break; // real and symmetric positive definite
        case -1: pardiso_mtype =  1; break; // real and structurally symmetric
        default:
            pardiso_mtype = 11; break; // real and nonsymmetric
    }

    for (size_t i = 0; i < 64; i++)
        pardiso_iparm[i] = 0;

    pardiso_iparm[1-1] = -1;   // not default
    pardiso_iparm[2-1] = 2;   // use METIS
    const char* var = getenv("OMP_NUM_THREADS");
    pardiso_iparm[3-1] = (var == NULL) ? 1 : atoi(var);
    msg_info("SparsePARDISOSolverInvertData") << "Using " << pardiso_iparm[2] << " thread(s), set OMP_NUM_THREADS environment variable to change.";
    pardiso_iparm[18-1] = -1;
    pardiso_iparm[24-1] = 1;
    pardiso_iparm[25-1] = 1;
    pardiso_iparm[38-1] = 0;
    pardiso_iparm[52-1] = 1;
    int solver = 0; /* use sparse direct solver */
    /* Numbers of processors, value of OMP_NUM_THREADS */

    F77_FUNC(pardisoinit) (pardiso_pt,  &pardiso_mtype, &solver, &pardiso_iparm[0], pardiso_dparm, &pardiso_initerr);

    switch(pardiso_initerr)
    {
        case 0:   msg_info("SparsePARDISOSolverInvertData")  << "License check was successful"; break;
        case -10: msg_error("SparsePARDISOSolverInvertData") << "No license file found"; break;
        case -11: msg_error("SparsePARDISOSolverInvertData") << "License is expired"; break;
        case -12: msg_error("SparsePARDISOSolverInvertData") << "Wrong username or hostname"; break;
        default:  msg_error("SparsePARDISOSolverInvertData") << "Unknown error " << pardiso_initerr; break;
    }
}


template<class TMatrix, class TVector>
SparsePARDISOSolver<TMatrix,TVector>::SparsePARDISOSolver()
    : d_symmetric( initData(&d_symmetric,1,"symmetric","type of the matrix: 0 = nonsymmetric arbitrary matrix, 1 = symmetric matrix, 2 = symmetric positive definite, -1 = structurally symmetric matrix") )
    , d_verbose( initData(&d_verbose,false,"verbose","print detailed information about the system being solved") )
    , d_exportDataToDir( initData(&d_exportDataToDir, std::string(""), "exportDataToFolder", "export data (matrix, RHS, solution) to files in directory specified by the parameter"))
    , d_iterativeSolverNumbering( initData(&d_iterativeSolverNumbering,false,"iterativeSolverNumbering","if true, the naming convention is incN_itM where N is the time step and M is the iteration within the time step (e.g. if used with the Newton-Raphson)") )
    , d_pardisoSchurComplement( initData(&d_pardisoSchurComplement,0,"pardisoSchurComplement","use Pardiso to compute Schur complement employing the partial factorization") )
    , d_trustRegionCoefficient( initData(&d_trustRegionCoefficient,0.0,"trustRegionCoefficient","add the value specified by the parameter on the diagonal of the matrix before inversion to perform the regularization") )
    , d_showTiming( initData(&d_showTiming,false,"showTiming","if true, show length of computation in microseconds") )
{
}

template<class TMatrix, class TVector>
void SparsePARDISOSolver<TMatrix,TVector>::init()
{
    this->f_listening.setValue(true);
    numIterationInStep = 0;
    numPrevNZ = 0;
    numActNZ = 0;
    Inherit::init();

    std::string exportDir=d_exportDataToDir.getValue();

    doExportData = false;
    if (!exportDir.empty()) {
        int status;
        status = mkdir(exportDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (status == 0 || errno == EEXIST)
            doExportData = true;
        else
            msg_error() << " Cannot create directory " << exportDir << " for exporting the solver data";

    }
    timeStep = -1;

    defaulttype::BaseMatrix* dummyM = NULL;
    pardisoInvertData = (SparsePARDISOSolverInvertData*) this->getMatrixInvertData(dummyM);

    patJL.clear();
    patJC.clear();
    solveNum=0;
}

template<class TMatrix, class TVector>
int SparsePARDISOSolver<TMatrix,TVector>::callPardiso(SparsePARDISOSolverInvertData* data, int phase, Vector* vx, Vector* vb)
{
    int maxfct = 1; // Maximum number of numerical factorizations
    int mnum = 1; // Which factorization to use
    int n = data->Mfiltered.rowSize();
    double* a = NULL;
    int* ia = NULL;
    int* ja = NULL;
    int* perm = NULL; // User-specified permutation vector
    int nrhs = 0; // Number of right hand sides
    int msglvl = (d_verbose.getValue())?1:0; // Print statistical information
    double* b = NULL;
    double* x = NULL;
    int error = 0;

    if (phase > 0)
    {
        ia = (int *) &(data->Mfiltered.getRowBegin()[0]);
        ja = (int *) &(data->Mfiltered.getColsIndex()[0]);
        a  = (double*) &(data->Mfiltered.getColsValue()[0]);

        if (vx)
        {
            nrhs = 1;
            x = vx->ptr();
            b = vb->ptr();
        }
    }
    //sout << "Solver phase " << phase << "..." << sendl;
    sofa::helper::AdvancedTimer::stepBegin("PardisoRealSolving");
    F77_FUNC(pardiso)(data->pardiso_pt, &maxfct, &mnum, &data->pardiso_mtype, &phase,
                      &n, a, ia, ja, perm, &nrhs,
                      data->pardiso_iparm, &msglvl, b, x, &error,  data->pardiso_dparm);
    sofa::helper::AdvancedTimer::stepEnd("PardisoRealSolving");
    printError(error, phase);
    return error;
}

template<class TMatrix, class TVector>
void SparsePARDISOSolver<TMatrix,TVector>::invert(Matrix& M)
{
    SparsePARDISOSolverInvertData * data = pardisoInvertData;
    if (data->pardiso_initerr)
        return;

    sofa::helper::AdvancedTimer::stepBegin("PardisoInvert");

    M.compress();

    data->Mfiltered.clear();
    if (d_symmetric.getValue() > 0)
    {
        data->Mfiltered.copyUpperNonZeros(M);
        double coeff = d_trustRegionCoefficient.getValue();
        if (fabs(coeff) > 0.0) {
            msg_info() << "Trust region coeff: " << coeff;
            for (size_t i = 0; i < (size_t)data->Mfiltered.colSize(); i++)
                data->Mfiltered.add(i,i,coeff);
        }
        data->Mfiltered.fullDiagonal();
        msg_info() << "Filtered upper part of M, nnz = " << data->Mfiltered.getRowBegin().back();
    }
    else if (d_symmetric.getValue() < 0)
    {
        data->Mfiltered.copyNonZeros(M);
        data->Mfiltered.fullDiagonal();
        msg_info() << "Filtered M, nnz = " << data->Mfiltered.getRowBegin().back();
    }
    else
    {
        data->Mfiltered.copyNonZeros(M);
        data->Mfiltered.fullRows();
        msg_info() << "Filtered M, nnz = " << data->Mfiltered.getRowBegin().back();
    }

    CompressedRowSparseMatrix<double>& refM = data->Mfiltered;
    numActNZ = refM.getRowBegin().back();

    //  Convert matrix from 0-based C-notation to Fortran 1-based notation (only for the purpose of solving)
    refM.shiftIndices(1);
    //data->pardiso_iparm[37] = d_schurComplementSize.getValue();
    if (numPrevNZ != numActNZ) {
        if (this->f_printLog.getValue())
            msg_info() << "Analyze the matrix, nnz = " << numActNZ;
        TIC
        if (callPardiso(data, 11))
            goto END;
        data->factorized = true;
        TOC("Invert: analyze")
    }
    numPrevNZ = numActNZ;

    if (this->f_printLog.getValue())
        msg_info() << "Factorize the matrix";

    if (doExportData) {
        char nm[100];
        numIterationInStep++;
        sprintf(nm, "step%04ld_iter%04d", timeStep, numIterationInStep);
        suffix = nm;
        saveSparseMatrix(refM, "M");
        solveNum = -1;
    }

    TIC
    if (callPardiso(data, 22)) {
        data->factorized = false;
        goto END;
    }
    TOC("Invert: factorize");

    refM.shiftIndices(-1);

    END:
    sofa::helper::AdvancedTimer::stepEnd("PardisoInvert");



}

template<class TMatrix, class TVector>
void SparsePARDISOSolver<TMatrix,TVector>::solve (Matrix& M, Vector& z, Vector& r)
{
    if (doExportData){
        solveNum++;
        std::string exportDir=d_exportDataToDir.getValue();
        std::ofstream f;
        char name[100];
        sprintf(name, "%s/rhs_PARD_%s_solve%03d.txt", exportDir.c_str(), suffix.c_str(), solveNum);
        f.open(name);
        f << r;
        f.close();
    }

    SparsePARDISOSolverInvertData * data = (SparsePARDISOSolverInvertData *) this->getMatrixInvertData(&M);

    if (data->pardiso_initerr) return;
    if (!data->factorized) return;

    CompressedRowSparseMatrix<double>& refM = data->Mfiltered;
    //CompressedRowSparseMatrix<double>& refJMJ = data->JMJfiltered;

    data->pardiso_iparm[7] = 0;       /* Max numbers of iterative refinement steps. */

    msg_info() << "Solve the system";
    bool success = false;
    sofa::helper::AdvancedTimer::stepBegin("PardisoSolve");
    refM.shiftIndices(1);
    if (!callPardiso(data, 33, &z, &r)) {
        success = true;
        refM.shiftIndices(-1);
    }
    sofa::helper::AdvancedTimer::stepEnd("PardisoSolve");

    if (success && doExportData){
        std::string exportDir=d_exportDataToDir.getValue();
        std::ofstream f;
        char name[100];
        sprintf(name, "%s/solution_PARD_%s_solve%03d.txt", exportDir.c_str(), suffix.c_str(), solveNum);
        f.open(name);
        f << z;
        f.close();
    }
}

template<class TMatrix, class TVector>
bool SparsePARDISOSolver<TMatrix,TVector>::addJMInvJtLocal(Matrix * M,ResMatrixType * result,const JMatrixType * J, double fact) {
    if (!d_pardisoSchurComplement.getValue()) {
        if (this->f_printLog.getValue())
            msg_info() << "Compute the Schur complement using iterated solve ";
        /*if (masterPardiso) {
            std::cerr << "[" << this->getName() << "] cannot use classical Schur computation in slave mode" << std::endl;
            return false;
        }*/
        TIC;
        bool res = Inherit::addJMInvJtLocal(M, result, J, fact);
        TOC("addJMInvJtLocal");
        //std::cout << "Schur complement using iterated solve: = " << *result << std::endl;
        return  res;
    }

    msg_info() << "Compute the Schur complement using incomplete factorization";

    SparsePARDISOSolverInvertData * data = pardisoInvertData;
    CompressedRowSparseMatrix<double>* sysM;
    //if (masterPardiso) {
    //    sysM = &masterPardiso->getPardisoInvertData()->Mfiltered;
    //} else {
    sysM = &pardisoInvertData->Mfiltered;
    //}

    CompressedRowSparseMatrix<double> refJMJ;

    size_t sizeM = sysM->colSize();
    size_t sizeSchur = J->rowSize();
    size_t sizeJMJ = sizeM + sizeSchur;

    msg_info() << "Size of the system matrix: " << sizeM;
    msg_info() << "Number of rows in J: " << sizeSchur;

    TIC;
    refJMJ.clear();
    refJMJ.resize(sizeJMJ, sizeJMJ);
    refJMJ.rowIndex = sysM->rowIndex;
    refJMJ.rowBegin = sysM->rowBegin;
    refJMJ.colsIndex = sysM->colsIndex;
    refJMJ.colsValue = sysM->colsValue;

    //std::cout << "J = " << *J << std::endl;

    helper::vector<int> patJLnew, patJCnew;
    for (typename sofa::component::linearsolver::SparseMatrix<double>::LineConstIterator jit1 = J->begin(); jit1 != J->end(); jit1++)
    {
        int l = jit1->first;
        patJLnew.push_back(l);
        for (typename sofa::component::linearsolver::SparseMatrix<double>::LElementConstIterator i1 = jit1->second.begin(); i1 != jit1->second.end(); i1++)
        {
            int c = i1->first;
            patJCnew.push_back(c);
            double val = i1->second;
            if (fabs(val) > 1e-12)
                refJMJ.set(c,l+sizeM,val);
        }
    }
    refJMJ.fullDiagonal();
    TOC("initialize JMJ");

    TIC
    bool patJremains = (patJLnew.size() == patJL.size()) && (patJCnew.size() == patJC.size());
    //std::cout << patJLnew.size() << " " <<  patJL.size() << " " << patJCnew.size() << " " << patJC.size() << std::endl;

    if (patJremains)
        for (size_t i = 0; i < patJLnew.size(); i++)
            patJremains = (patJLnew[i] == patJL[i]);

    if (patJremains)
        for (size_t i = 0; i < patJCnew.size(); i++)
            patJremains = (patJCnew[i] == patJC[i]);

    if (!patJremains) {
        patJL = patJLnew;
        patJC = patJCnew;
        msg_info() << "Pattern of J has changed!";
    }

    TOC("J pattern check");

    refJMJ.shiftIndices(1);
    int maxfct = 1; // Maximum number of numerical factorizations
    int mnum = 1; // Which factorization to use
    int n = refJMJ.rowSize();
    int* ia = (int *) &(refJMJ.getRowBegin()[0]);
    int* ja = (int *) &(refJMJ.getColsIndex()[0]);
    double* a  = (double*) &(refJMJ.getColsValue()[0]);
    int* perm = NULL; // User-specified permutation vector
    int nrhs = 0; // Number of right hand sides
    int msglvl = (d_verbose.getValue())?10:0; // Print statistical information
    double* b = NULL;
    double* x = NULL;
    int error = 0;
    int phase;

    if (doExportData) {
        saveSparseMatrix(refJMJ, "JMJ");
    }


    sofa::helper::AdvancedTimer::stepBegin("PardisoPartialFactorization");
    phase=11;
    //if (!patJremains) {
    TIC;
    data->pardiso_iparm[38-1] = sizeSchur;
    msg_info() << "Size schur: " << sizeSchur;
    F77_FUNC(pardiso)(data->pardiso_pt, &maxfct, &mnum, &data->pardiso_mtype, &phase,
                      &n, a, ia, ja, perm, &nrhs,
                      data->pardiso_iparm, &msglvl, b, x, &error,  data->pardiso_dparm);
    TOC("Analyze with Schur");
    //}

    phase=22;
    TIC
        F77_FUNC(pardiso)(data->pardiso_pt, &maxfct, &mnum, &data->pardiso_mtype, &phase,
                          &n, a, ia, ja, perm, &nrhs,
                          data->pardiso_iparm, &msglvl, b, x, &error,  data->pardiso_dparm);
    TOC("Factorize with Schur");
    sofa::helper::AdvancedTimer::stepEnd("PardisoPartialFactorization");

    TIC
    int schurNNZ = data->pardiso_iparm[39-1];
    int *iS = new int[sizeSchur+1];
    int *jS = new int[schurNNZ];
    double *S = new double[schurNNZ];
    F77_FUNC(pardiso_get_schur)(data->pardiso_pt, &maxfct, &mnum, &data->pardiso_mtype, S, iS, jS);
    TOC("Compute Schur");
    if (error < 0)
        printError(error, phase);

    /*std::cout << "Schur newway" << std::endl;
    for (int i = 0; i < schurNNZ; i++)
        std::cout << S[i] << " ";
    std::cout << std::endl;*/

    if (doExportData){
        std::string exportDir=d_exportDataToDir.getValue();
        std::ofstream f;
        char name[100];
        sprintf(name, "%s/schurIn_PARD_%s.txt", exportDir.c_str(), suffix.c_str());
        f.open(name);
        f << *result;
        f.close();
    }

    TIC
    int rw = 0;
    for (int i = 0; i < iS[sizeSchur]-1; i++) {
        if (iS[rw] == i+1)
            rw++;
        int l = rw-1;
        int c = jS[i]-1;
        double v = -fact*S[i];
        result->add(l, c, v);
    }

    for (size_t i = 0; i < (size_t)result->colSize(); i++)
        for (size_t j = 0; j < i; j++)
            result->add(i,j,result->element(j,i));
    delete[] iS;
    delete[] jS;
    delete[] S;
    TOC("Convert Shur to res");

    //if (!masterPardiso) {
    TIC;
    data->pardiso_iparm[38-1] = 0;
    sysM->shiftIndices(1);
    callPardiso(data, 12);
    sysM->shiftIndices(-1);
    TOC("Refactorize without Schur");
    //}

    if (doExportData){
        std::string exportDir=d_exportDataToDir.getValue();
        std::ofstream f;
        char name[100];
        sprintf(name, "%s/schurOut_PARD_%s.txt", exportDir.c_str(), suffix.c_str());
        f.open(name);
        f << *result;
        f.close();
    }

    //std::cout << "Schur size: " << sizeSchur << std::endl;

    /*double sm1 = 0;
    for (size_t i = 0; i < sizeSchur; i++)
        for (size_t j = 0; j < sizeSchur; j++)
            sm1 += fabs(result->element(i,j));*/

    //std::cout << "SM1 = " <<sm1 << " SM2= " << sm2 << std::endl;

    return true;
}

template<class TMatrix, class TVector>
void SparsePARDISOSolver<TMatrix,TVector>::printError(const int error, const int phase) {
    const char* msg = NULL;
    switch(error)
    {
        case 0: break;
        case -1: msg="Input inconsistent"; break;
        case -2: msg="Not enough memory"; break;
        case -3: msg="Reordering problem"; break;
        case -4: msg="Zero pivot, numerical fact. or iterative refinement problem"; break;
        case -5: msg="Unclassified (internal) error"; break;
        case -6: msg="Preordering failed (matrix types 11, 13 only)"; break;
        case -7: msg="Diagonal matrix problem"; break;
        case -8: msg="32-bit integer overflow problem"; break;
        case -10: msg="No license file pardiso.lic found"; break;
        case -11: msg="License is expired"; break;
        case -12: msg="Wrong username or hostname"; break;
        case -100: msg="Reached maximum number of Krylov-subspace iteration in iterative solver"; break;
        case -101: msg="No sufficient convergence in Krylov-subspace iteration within 25 iterations"; break;
        case -102: msg="Error in Krylov-subspace iteration"; break;
        case -103: msg="Break-Down in Krylov-subspace iteration"; break;
        default: msg="Unknown error"; break;
    }
    if (msg)
        msg_error() << "Solver phase " << phase << ": ERROR " << error << " : " << msg;
}

template<class TMatrix, class TVector>
void SparsePARDISOSolver<TMatrix,TVector>::saveSparseMatrix(CompressedRowSparseMatrix<double>& _M, std::string _name) {
    int n = _M.rowSize();
    int* ia = (int *) &(_M.getRowBegin()[0]);
    int* ja = (int *) &(_M.getColsIndex()[0]);
    double* a  = (double*) &(_M.getColsValue()[0]);

    std::string exportDir=d_exportDataToDir.getValue();
    std::ofstream f;
    char name[100];
    sprintf(name, "%s/spmat%s_PARD_%s.txt", exportDir.c_str(), _name.c_str(), suffix.c_str());
    f.open(name);
    int rw = 0;
    for (int i = 0; i < ia[n]-1; i++) {
        if (ia[rw] == i+1)
            rw++;
        f << rw << " " << ja[i] << " " << a[i] << std::endl;
    }
    f.close();
}

template<class TMatrix, class TVector>
void SparsePARDISOSolver<TMatrix,TVector>::printSparseMatrix(CompressedRowSparseMatrix<double>& _M, std::string _name) {
    int n = _M.rowSize();
    int* ia = (int *) &(_M.getRowBegin()[0]);
    int* ja = (int *) &(_M.getColsIndex()[0]);
    double* a  = (double*) &(_M.getColsValue()[0]);

    int rw = 0;
    for (int i = 0; i < ia[n]-1; i++) {
        if (ia[rw] == i+1)
            rw++;
        std::cout << _name << " = \n" << rw << " " << ja[i] << " " << a[i] << std::endl;
    }
}

template<class TMatrix, class TVector>
void SparsePARDISOSolver<TMatrix,TVector>::handleEvent(core::objectmodel::Event *event)
{
    if (dynamic_cast<sofa::simulation::AnimateBeginEvent *>(event)) {
        if (doExportData) {
            double dt = Inherit::getContext()->getDt();
            double actualTime = Inherit::getContext()->getTime();
            timeStep = lround(actualTime/dt);
            numIterationInStep = -1;
        }
    }
}


SOFA_DECL_CLASS(SparsePARDISOSolver)

int SparsePARDISOSolverClass = core::RegisterObject("Direct linear solvers implemented with the PARDISO library")
                                   .add< SparsePARDISOSolver< CompressedRowSparseMatrix<double>,FullVector<double> > >()
.add< SparsePARDISOSolver< CompressedRowSparseMatrix< defaulttype::Mat<3,3,double> >,FullVector<double> > >(true)
.addAlias("PardisoSolver")
;

} // namespace linearsolver

} // namespace component

} // namespace sofa

