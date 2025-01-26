#ifndef MEMORY_MANAGER_FINITE_DIFFERENCE2D_
#define MEMORY_MANAGER_FINITE_DIFFERENCE2D_
#include "DeviceMemory.cuh"


namespace numerical
{
	class MemoryManagerFiniteDifference2D

	{
	public:

		explicit MemoryManagerFiniteDifference2D(size_t sz);

		void copyHostInputDataToDevice(const float* data, size_t sz);
		void copyDeviceOutputDataToHost(float* data, size_t sz);

		float* getInputData()
		{
			return inputData_.d_memory;
		}

		float* getOutputData()
		{
			return outputData_.d_memory;
		}
		MemoryManagerFiniteDifference2D() = delete;


		MemoryManagerFiniteDifference2D(const MemoryManagerFiniteDifference2D&) = delete;
		MemoryManagerFiniteDifference2D& operator=(const MemoryManagerFiniteDifference2D&) = delete;
		MemoryManagerFiniteDifference2D(MemoryManagerFiniteDifference2D&& other) = delete;
		MemoryManagerFiniteDifference2D& operator=(MemoryManagerFiniteDifference2D&& other) = delete;

	private:
		DeviceMemory<float> inputData_;
		DeviceMemory<float> outputData_;

	};
}


#endif