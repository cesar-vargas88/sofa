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

#include <cuda_runtime.h>
#include "cublas_v2.h"

#include <ctime>
#include <iostream>

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
    
    void Transpose(double* src, double* dst, unsigned n, unsigned m) {

        // Allocate device memory
        double* d_src;
        double* d_dst;

        // Allocate device memory
        if (cudaMalloc(&d_src, sizeof(double) * n * m) != cudaSuccess) std::cout << "cudaMalloc failed!" << std::endl;
        if (cudaMalloc(&d_dst, sizeof(double) * m * n) != cudaSuccess) std::cout << "cudaMalloc failed!" << std::endl;

        if (cudaMemcpy(d_src, src, m * n * sizeof(double), cudaMemcpyHostToDevice)) std::cout << "cudaMemcpy failed!" << std::endl;
        //cudaDeviceSynchronize();

        // cuBLAS handle
        cublasHandle_t handle;

        if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS)
            std::cout << "CUBLAS initialization failed" << std::endl;

        // Scalaing factors
        double alpha = 1.0;
        double beta = 0.0;

        // Tranpose d_matrix2
        cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N, n, m, &alpha, d_src, m, &beta, d_src, n, d_dst, n);
        //cudaDeviceSynchronize();

        // Copy back the three matrices
        cudaMemcpy(dst, d_dst, sizeof(double) * m * n, cudaMemcpyDeviceToHost);
        //cudaDeviceSynchronize();

        // Free our memory
        cudaFree(d_src);
        cudaFree(d_dst);

        cublasDestroy(handle);
    }

    void MultiplyThrust(double* m1, double* m2, double* result, unsigned m, unsigned n, unsigned r) {

        thrust::device_vector<double> matrix1(m1, m1 + n * m);
        thrust::device_vector<double> matrix2(m2, m2 + m * r);
        thrust::device_vector<double> matrix_result(n * r, 0);

        thrust::transform(thrust::counting_iterator<unsigned>(0),
            thrust::counting_iterator<unsigned>(n * r),
            matrix_result.begin(),
            dp<double>(thrust::raw_pointer_cast(matrix1.data()), thrust::raw_pointer_cast(matrix2.data()), m, n, r));

        cudaDeviceSynchronize();

        thrust::copy(matrix_result.begin(), matrix_result.end(), result);
    }

    void MultiplyCUBLAS(double* m1, double* m2, double* result, unsigned m, unsigned n, unsigned r) {

        // Allocate device memory
        double* d_matrix1;
        double* d_matrix2;
        double* d_result;

        // Allocate device memory
        if (cudaMalloc(&d_matrix1, sizeof(double) * n * m) != cudaSuccess) std::cout << "cudaMalloc failed!" << std::endl;
        if (cudaMalloc(&d_matrix2, sizeof(double) * m * r) != cudaSuccess) std::cout << "cudaMalloc failed!" << std::endl;
        if (cudaMalloc(&d_result, sizeof(double) * n * r) != cudaSuccess) std::cout << "cudaMalloc failed!" << std::endl;

        // Copy host to device memory
        if (cudaMemcpy(d_matrix1, m1, n * m * sizeof(double), cudaMemcpyHostToDevice)) std::cout << "cudaMemcpy failed!" << std::endl;
        if (cudaMemcpy(d_matrix2, m2, m * r * sizeof(double), cudaMemcpyHostToDevice)) std::cout << "cudaMemcpy failed!" << std::endl;
        //cudaDeviceSynchronize();

        // cuBLAS handle
        cublasHandle_t handle;

        if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS)
            std::cout << "CUBLAS initialization failed" << std::endl;

        // Scalaing factors
        double alpha = 1.0;
        double beta = 0.0;

        // Calculate: c = (alpha*a) * b + (beta*c)
        // nxr = nxm * mxr
        // Signature: handle, operation, operation, n, r, m, alpha, A, lda, B, ldb,
        // beta, C, ldc    
        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, r, m, &alpha, d_matrix2, n, d_matrix1, m, &beta, d_result, r);
        //cudaDeviceSynchronize();

        // Copy back the three matrices
        cudaMemcpy(result, d_result, sizeof(double) * n * r, cudaMemcpyDeviceToHost);
        //cudaDeviceSynchronize();

        // Free our memory
        cudaFree(d_matrix1);
        cudaFree(d_matrix2);
        cudaFree(d_result);

        cublasDestroy(handle);
    }

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif



