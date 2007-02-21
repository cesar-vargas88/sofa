#ifndef SOFA_GPU_CUDA_CUDAUNIFORMMASS_H
#define SOFA_GPU_CUDA_CUDAUNIFORMMASS_H

#include "CudaTypes.h"
#include <sofa/component/mass/UniformMass.h>

namespace sofa
{

namespace component
{

namespace mass
{

// -- Mass interface
template <>
void UniformMass<gpu::cuda::CudaVec3fTypes, float>::addMDx(VecDeriv& res, const VecDeriv& dx);

template <>
void UniformMass<gpu::cuda::CudaVec3fTypes, float>::accFromF(VecDeriv& a, const VecDeriv& f);

template <>
void UniformMass<gpu::cuda::CudaVec3fTypes, float>::addForce(VecDeriv& f, const VecCoord&, const VecDeriv&);

} // namespace mass

} // namespace component

} // namespace sofa

#endif
