/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_LINEARSOLVER_SparseLDLSolver_INL
#define SOFA_COMPONENT_LINEARSOLVER_SparseLDLSolver_INL

#include <SofaSparseSolver/SparseLDLSolver.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include "sofa/helper/system/thread/CTime.h"
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <cmath>
#include <sofa/helper/system/thread/CTime.h>
#include <SofaBaseLinearSolver/CompressedRowSparseMatrix.inl>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <string>

#include "CudaMatrix.h"
#include <ctime>

namespace sofa {

namespace component {

namespace linearsolver {

template<class TMatrix, class TVector, class TThreadManager>
SparseLDLSolver<TMatrix,TVector,TThreadManager>::SparseLDLSolver()
    : numStep(0)
    , f_saveMatrixToFile( initData(&f_saveMatrixToFile, false, "savingMatrixToFile", "save matrix to a text file (can be very slow, as full matrix is stored"))
    , d_filename( initData(&d_filename, std::string("MatrixInLDL_%04d.txt"),"savingFilename", "Name of file where system matrix (mass, stiffness and damping) will be stored."))
    , d_precision( initData(&d_precision, 6, "savingPrecision", "Number of digits used to store system's matrix. Default is 6."))
{}

template<class TMatrix, class TVector, class TThreadManager>
void SparseLDLSolver<TMatrix,TVector,TThreadManager>::solve (Matrix& M, Vector& z, Vector& r) {
    Inherit::solve_cpu(&z[0],&r[0],(InvertData *) this->getMatrixInvertData(&M));
}

template<class TMatrix, class TVector, class TThreadManager>
void SparseLDLSolver<TMatrix,TVector,TThreadManager>::invert(Matrix& M) {
    if (f_saveMatrixToFile.getValue())
    {
        std::ofstream f;
        std::string name=d_filename.getValue().c_str();
        if (d_filename.getValue().find("%d"))
        {
            char bname[100];
            snprintf(bname, 100, d_filename.getValue().c_str(), numStep);
            name=bname;
        }
        f.open(name);
        f << std::scientific << std::setprecision(d_precision.getValue()) << M;
        f.close();
    }

    Mfiltered.copyNonZeros(M);
    Mfiltered.compress();

    int n = M.colSize();

    int * M_colptr = (int *) &Mfiltered.getRowBegin()[0];
    int * M_rowind = (int *) &Mfiltered.getColsIndex()[0];
    Real * M_values = (Real *) &Mfiltered.getColsValue()[0];

    if(M_colptr==nullptr || M_rowind==nullptr || M_values==nullptr || Mfiltered.getRowBegin().size() < (size_t)n )
    {
        msg_warning() << "Invalid Linear System to solve. Please insure that there is enough constraints (not rank deficient)." ;
        return ;
    }

    Inherit::factorize(n,M_colptr,M_rowind,M_values,(InvertData *) this->getMatrixInvertData(&M));

    numStep++;
}


void printMatrix(unsigned N, unsigned M, double* array) {
    for (int i = 0; i < N * M; i++) {
        std::cout << "\t" << array[i];

        if (!((i + 1) % M))
            std::cout << std::endl;
    }
    std::cout << std::endl;
}

/// Default implementation of Multiply the inverse of the system matrix by the transpose of the given matrix, and multiply the result with the given matrix J
template<class TMatrix, class TVector, class TThreadManager>
bool SparseLDLSolver<TMatrix, TVector, TThreadManager>::addJMInvJtLocal(TMatrix* M, ResMatrixType* result, const JMatrixType* J, double fact) {

    if (J->rowSize() == 0) return true;

    InvertData* data = (InvertData*)this->getMatrixInvertData(M);

    Jdense.clear();
    Jdense.resize(J->rowSize(), data->n);
    Jminv.resize(J->rowSize(), data->n);

    for (typename SparseMatrix<Real>::LineConstIterator jit = J->begin(), jitend = J->end(); jit != jitend; ++jit) {

        Real* line = Jdense[jit->first];

        for (typename SparseMatrix<Real>::LElementConstIterator it = jit->second.begin(), itend = jit->second.end(); it != itend; ++it)
            line[data->invperm[it->first]] = it->second;
    }

    //Solve the lower triangular system
    for (unsigned c = 0; c < (unsigned)J->rowSize(); c++) {

        Real* line = Jdense[c];

        for (int j = 0; j < data->n; j++) {
            for (int p = data->LT_colptr[j]; p < data->LT_colptr[j + 1]; p++)
                line[j] -= (data->LT_values[p]) * line[data->LT_rowind[p]];
        }
    }

    //Apply diagonal
    for (unsigned j = 0; j < (unsigned)J->rowSize(); j++) {
        Real* lineD = Jdense[j];
        Real* lineM = Jminv[j];
        for (unsigned i = 0; i < (unsigned)J->colSize(); i++) {
            lineM[i] = lineD[i] * data->invD[i];
        }
    }

    clock_t begin_time = clock();
    std::cout << begin_time << " , SparseLDLSolver - addJMInvJtLocal - multiply the result with the given matrix J, " << std::endl;

    const unsigned n = (unsigned)J->rowSize();
    const unsigned m = (unsigned)J->colSize();
    const unsigned r = (unsigned)J->rowSize();

    double* response = new double[n * r];
    double* Jminv_T  = new double[m * r];
    
    clock_t begin_time0 = clock();
    std::cout << begin_time0 << " , Transpose, " << std::endl;

    Transpose(Jminv[0] , Jminv_T , n, m);

    std::cout << clock() << " , Transpose, " << double(clock() - begin_time0) / CLOCKS_PER_SEC << " seconds " << std::endl;

    clock_t begin_time1 = clock();
    std::cout << begin_time1 << " , MultiplyCUBLAS, " << std::endl;

    MultiplyCUBLAS(Jdense[0], Jminv_T, response, m, n, r);

    std::cout << clock() << " , MultiplyCUBLAS, " << double(clock() - begin_time1) / CLOCKS_PER_SEC << " seconds " << std::endl;

    clock_t begin_time2 = clock();
    std::cout << begin_time2 << " , Copy result, " << std::endl;

    for (unsigned j = 0; j < n ; j++) {
        for (unsigned i = 0; i < n ; i++)
            result->add(j, i, response[j * n + i] * fact);
    }

    std::cout << clock() << " , Copy result, " << double(clock() - begin_time2) / CLOCKS_PER_SEC << " seconds " << std::endl;

    delete[] response;
    delete[] Jminv_T;
    
    std::cout << clock() << " , SparseLDLSolver - addJMInvJtLocal - multiply the result with the given matrix J , " << " Matrix "<< n << "x" << m << " time: " << double(clock() - begin_time) / CLOCKS_PER_SEC << " seconds "<<std::endl;

    return true;
}

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
