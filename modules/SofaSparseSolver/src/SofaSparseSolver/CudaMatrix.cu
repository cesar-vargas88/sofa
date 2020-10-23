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
#ifndef SOFA_COMPONENT_LINEARSOLVER_CUDAMATRIX_H
#define SOFA_COMPONENT_LINEARSOLVER_CUDAMATRIX_H

#include "CudaMatrix.h"

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/inner_product.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{
    template<class T>
    struct dp
    {
        T* A, * B;
        int m, n, r;

        dp(T* _A, T* _B, int _m, int _n, int _r) : A(_A), B(_B), m(_m), n(_n), r(_r) {};

        __host__ __device__
            T operator()(size_t idx) {

            T sum = 0.0f;
            int row = idx / r;
            int col = idx - (row * r); // cheaper modulo

            for (int i = 0; i < m; i++)
                sum += A[row * m + i] * B[col * m + i];

            return sum;
        }
    };

    void MultiplyMattix(double* m1, double* m2, double* result, unsigned m, unsigned n, unsigned r) {

        thrust::device_vector<double> data_     (m1 , m1+n*m);
        thrust::device_vector<double> other_    (m2 , m2+m*r);
        thrust::device_vector<double> result_   (result,result  + n * r);

        //dp<double> _dp(thrust::raw_pointer_cast(data_.data()), thrust::raw_pointer_cast(other_.data()), m, n, r);
        
        thrust::transform(thrust::counting_iterator<unsigned>(0),
            thrust::counting_iterator<unsigned>(n* r),
            result_.begin(),
            dp<double>(thrust::raw_pointer_cast(data_.data()), thrust::raw_pointer_cast(other_.data()), m, n, r));
    
        cudaDeviceSynchronize();

        for (unsigned i = 0; i < n * r; i++)
            result[i] = result_[i];
    }

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif