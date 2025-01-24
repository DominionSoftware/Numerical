#include "MemoryManagerFiniteDifference2D.cuh"

#include <stdexcept>


using namespace numerical;

MemoryManagerFiniteDifference2D::MemoryManagerFiniteDifference2D(size_t sz)
{
	inputData_.allocate(sz);
	outputData_.allocate(sz);

}

void MemoryManagerFiniteDifference2D::copyHostInputDataToDevice(float* data, size_t sz)
{
	cudaError_t error = cudaMemcpy(getInputData(), data, sz, cudaMemcpyHostToDevice);

	if (error != cudaSuccess)
	{
		throw std::runtime_error(cudaGetErrorString(error));
	}
}

void MemoryManagerFiniteDifference2D::copyDeviceOutputDataToHost(float* data, size_t sz)
{
	cudaError_t error = cudaMemcpy(data,getOutputData(), sz, cudaMemcpyDeviceToHost);

	if (error != cudaSuccess)
	{
		throw std::runtime_error(cudaGetErrorString(error));
	}
}
