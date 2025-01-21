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


#include <iostream>
#include <numbers>
#include <ostream>
#include <gtest/gtest.h>
#include "Roots.h"


 


TEST(NumericalTestSuite, TestRootFinding)
{
	
    const std::function<double(double)>  f = [](double x)->double
        {
            double x2 = (x * x);

            return (x2 * x) - x2 - x - 1;
        };

	
    std::tuple<numerical::FZeroErrors, double> result = numerical::BrentZero<double>(f, 1.0, 2.0);
    std::cout << static_cast<int>(std::get<0>(result)) << " " << std::get<1>(result) << std::endl;



    const std::function<float(float)>  ff = [](float x)->float
        {
            float x2 = (x * x);

            return (x2 * x) - x2 - x - 1;
        };


    std::tuple<numerical::FZeroErrors, float> resultf = numerical::BrentZero<float>(ff, 1.0f, 2.0f);
    std::cout << static_cast<int>(std::get<0>(resultf)) << " " << std::get<1>(resultf) << std::endl;

}

TEST(NumericalTestSuite, TestRootFindingNewton)
{

    const std::function<double(double)>  ff = [](double x)->double
        {
            double x2 = (x * x);

            return (x2 * x) - x2 - x - 1;
        };

    const std::function<double(double)>  fPrime = [](double x)->double
        {
            

            return (3 * x * x) - (x + x) - 1;
        };


    double root = numerical::NewtonRhapson(400,1.4, ff, fPrime);
    std::cout << root << std::endl;

}