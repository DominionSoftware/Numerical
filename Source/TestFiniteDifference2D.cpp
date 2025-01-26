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

 

 
#include "RunFiniteDifference2D.h"
#include "linspace.h"
 
#include "WriteVector.h"
#include "gtest/gtest.h"
#include <numbers>

using namespace numerical;




TEST(NumericalTestSuite, TestFiniteDifference2D1)
{
    // Function to test
    auto func = [](float x)->float
        {
           return std::cos(x);
        };

    namespace ranges = std::ranges;
    std::vector<float> xData = numerical::Linspace<float>(-std::numbers::pi, 2 * std::numbers::pi, 256);
    std::vector<float> inputData;
    std::vector<float> outputData;
    inputData.resize(256 * 256);
    outputData.resize(256 * 256);
    auto iter = inputData.begin();

    for (size_t i = 0; i < 256; i++)
    {
        ranges::transform(xData.begin(), xData.end(), iter, func);
        iter += 256;

    };
    numerical::RunFiniteDifference2D(&inputData[0], 256, 256, &outputData[0]);



}

