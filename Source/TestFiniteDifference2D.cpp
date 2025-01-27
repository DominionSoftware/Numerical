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
#include "Linspace.h"
 
#include "WriteVector.h"
#include "gtest/gtest.h"
#include <numbers>
#include <cmath>
#include <filesystem>

using namespace numerical;




TEST(NumericalTestSuite, TestFiniteDifference2D1)
{
    // Function to test
    auto func = [](float x)->float
        {
           return std::cos(x);
        };

    auto test_path = []() ->std::filesystem::path
        {
            return std::filesystem::current_path() / "finitedDifference2D_outout" / "";
        };


    namespace ranges = std::ranges;
    std::vector<float> xData = numerical::Linspace<float>(-std::numbers::pi_v<float>, 2 * std::numbers::pi_v<float>, 256);
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


    std::filesystem::create_directories(test_path());


    auto saveInputFileName = test_path().replace_filename("inputData.csv");

   
    std::ofstream inputDataFile;

    inputDataFile.open(saveInputFileName.string());

    if (!inputDataFile.is_open())
    {
        throw std::runtime_error("Unable to open file.");
    }

    std::vector<float>::const_iterator inputIter = inputData.begin();
    
    for (size_t i = 0; i < 256; i++)

    {
        numerical::FileIO<float>::FileWriteResult result = numerical::FileIO<float>::WriteVector(inputDataFile, inputIter, inputIter + 256);
        if (result == numerical::FileIO<float>::FileWriteResult::FileWriteErr)
        {
            throw std::runtime_error("Error Writing File.");
        }
        inputIter += 256;
    }
    inputDataFile.close();

    auto saveOutputFileName = test_path().replace_filename("outputData.csv");

    std::ofstream outputDataFile;

    outputDataFile.open(saveOutputFileName.string());

    if (!outputDataFile.is_open())
    {
        throw std::runtime_error("Unable to open file.");
    }

    std::vector<float>::const_iterator outputIter = outputData.begin();

    for (size_t i = 0; i < 256; i++)

    {
        numerical::FileIO<float>::FileWriteResult result = numerical::FileIO<float>::WriteVector(outputDataFile, outputIter, outputIter + 256);
        if (result == numerical::FileIO<float>::FileWriteResult::FileWriteErr)
        {
            throw std::runtime_error("Error Writing File.");
        }
        outputIter += 256;
    }
    outputDataFile.close();
   

}

