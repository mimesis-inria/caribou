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
// Authors: Hadrien Courtecuisse, Igor Peterlik

#ifndef SOFA_COMPONENT_LINEARSOLVER_SparsePARDISOSolver_H
#define SOFA_COMPONENT_LINEARSOLVER_SparsePARDISOSolver_H

#include <SofaBaseLinearSolver/MatrixLinearSolver.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <SofaBaseLinearSolver/SparseMatrix.h>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaBaseLinearSolver/CompressedRowSparseMatrix.h>
#include <sofa/helper/map.h>

#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <fstream>

#define TIC  if (this->d_showTiming.getValue()) this->startTime = double(this->timer->getTime()) ;
#define TOC(arg) if (this->d_showTiming.getValue()) { this->stopTime = double(this->timer->getTime()); \
                 std::cout << "[" << this->getName() << "] WTIME: " << arg << " " << this->stopTime - this->startTime << std::endl; }

#define TOCTIC(arg) if (this->d_showTiming.getValue()) { this->stopTime = double(this->timer->getTime()); \
    std::cout << "[" << this->getName() << "] WTIME: " << arg << " " << this->stopTime - this->startTime << std::endl; \
                 this->startTime = double(this->timer->getTime()); }


namespace sofa
{

namespace component
{

namespace linearsolver
{

/// Direct linear solvers implemented with the PARDISO library
template<class TMatrix, class TVector>
class SparsePARDISOSolver : public sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(SparsePARDISOSolver,TMatrix,TVector),SOFA_TEMPLATE2(sofa::component::linearsolver::MatrixLinearSolver,TMatrix,TVector));

    typedef TMatrix Matrix;
    typedef TVector Vector;
    typedef typename Matrix::Real Real;
    typedef sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector> Inherit;
    typedef typename Inherit::ResMatrixType ResMatrixType;
    typedef typename Inherit::JMatrixType JMatrixType;

    //Data< helper::vector<std::string> > d_options;
    Data<int> d_symmetric;
    Data<bool> d_verbose;
    Data<std::string> d_exportDataToDir;
    Data<bool> d_iterativeSolverNumbering;
    Data<int> d_pardisoSchurComplement;
    Data<double> d_trustRegionCoefficient;
    Data<bool> d_showTiming;

    //SingleLink<MyType, MyType, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> masterPardisoLink;
    //MyType* masterPardiso;

    SparsePARDISOSolver();
    ~SparsePARDISOSolver();
    virtual void init();

    void solve (Matrix& M, Vector& x, Vector& b);
    void invert(Matrix& M);
    virtual bool addJMInvJtLocal(Matrix * M,ResMatrixType * result,const JMatrixType * J, double fact);
    void saveSparseMatrix(CompressedRowSparseMatrix<double>& M, std::string _name);
    void printSparseMatrix(CompressedRowSparseMatrix<double>& M, std::string _name = "");

    void handleEvent(core::objectmodel::Event *event);

    MatrixInvertData * createInvertData()
    {
        return new SparsePARDISOSolverInvertData(d_symmetric.getValue(),std::cout,std::cerr);
    }

private:
    int solveNum;
    void printError(const int error, const int phase);

    sofa::helper::system::thread::CTime *timer;
    double startTime, stopTime;

    bool initialized;

    helper::vector<int> patJL;
    helper::vector<int> patJC;

protected:
    bool doExportData;
    std::string suffix;

    int numIterationInStep;
    long int timeStep;
    int numPrevNZ, numActNZ;
    class SparsePARDISOSolverInvertData : public MatrixInvertData
    {
    public :
        CompressedRowSparseMatrix<double> Mfiltered;
        SparsePARDISOSolver<Matrix,Vector>* solver;
        void*  pardiso_pt[64];
        int    pardiso_iparm[64];
        double pardiso_dparm[64];
        int pardiso_initerr;
        int pardiso_mtype;
        bool factorized;

        SparsePARDISOSolverInvertData(int d_symmetric,std::ostream & sout,std::ostream & serr);

        ~SparsePARDISOSolverInvertData()
        {
            if (solver && pardiso_initerr == 0)
            {
                solver->callPardiso(this, -1);  // Release internal memory.
            }
        }

    };

    SparsePARDISOSolverInvertData* pardisoInvertData;


    int callPardiso(SparsePARDISOSolverInvertData* data, int phase, Vector* vx = NULL, Vector* vb = NULL);
public:
    SparsePARDISOSolverInvertData* getPardisoInvertData() {
        return this->pardisoInvertData;
    }
};

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
