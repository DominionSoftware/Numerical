#ifndef DEVICE_MEMORY_
#define DEVICE_MEMORY_

#include "cuda_runtime.h"
#include <exception>

namespace numerical
{
	template<typename T>
	class DeviceMemory
	{
	public:
		void allocate(size_t sz)
		{
			cudaError_t e = cudaMalloc(&d_memory, sz);

			if (e != cudaSuccess)
			{
				throw std::bad_alloc();
			}
		}


		DeviceMemory() : d_memory{}
		{
		}

		~DeviceMemory()
		{
			if (d_memory != nullptr)
			{
				cudaFree(d_memory);
			}
		}


	 
		T* d_memory;


	};
}



#endif