/*

                      Copyright (c) 2025 Dominion Software, Inc.
                      Copyright (c) 2025 Rick Frank rfrank@dominionsw.com.
                      
                                 All Rights Reserved.
                      
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify,   
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to  
permit persons to whom the Software is furnished to do so, subject to the following 
conditions:     
                                 
The above copyright notice and this permission notice shall be included in ALL copies
or ANY portions of the Software.  
                         
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,  
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A        
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                        

*/

#include "cuda_runtime.h"
#include "cuda.h"

__device__ const float hDiff[3][3] =
{
    {1, 0, -1},
    {1, 0, -1},
    {1, 0, -1}
};

__device__ float getDataSafe(
    const float* data,
    int x,
    int y,
    int width,
    int height)
{
    // Clamp coordinates to boundaries
    x = max(0, min(x, width - 1));
    y = max(0, min(y, height - 1));

    // Return  value at clamped coordinates
    return data[y * width + x];
}

__global__ void finiteDiffKernel(
    const float* input,
    float* output,
    int width,
    int height)
{
   


    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x >= width || y >= height)
    {
        return;
    }
    float d{ };
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            float v = getDataSafe(input, x + j, y + i, width, height);
            d += v * hDiff[i + 2][j + 2];
        }
    }

    output[y * width + x] = d;

}
