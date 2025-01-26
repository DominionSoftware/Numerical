#include "RunFiniteDifference2D.h"

#include <iostream>

#include "MemoryManagerFiniteDifference2D.cuh"
#include "FiniteDifference2D.cuh"

 

namespace numerical
{


	void RunFiniteDifference2D(const float* inputData, int width, int height, float* outputData)
	{
		try
		{
			size_t sz = width * height * sizeof(float);
			MemoryManagerFiniteDifference2D memory(width * height * sizeof(float));

			memory.copyHostInputDataToDevice(inputData, sz);

			dim3 blockDim(16, 16);
			dim3 gridDim((width + blockDim.x - 1) / blockDim.x,

				(height + blockDim.y - 1) / blockDim.y);
			finiteDiffKernel <<<gridDim, blockDim >>> (memory.getInputData(), memory.getOutputData(), width, height);

			memory.copyDeviceOutputDataToHost(outputData, sz);



		}
		catch (std::exception& ex)
		{
			std::cout << ex.what() << std::endl;
		}

	}
}