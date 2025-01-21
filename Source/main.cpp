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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <iso646.h>

#include "GenerateLegendreCoeffs.h"


int main()
{
	Eigen::VectorXd x;
	Eigen::VectorXd w;

	ComputeAbscissasAndWeights(6, x, w);
	std::ofstream outputFileE;
	outputFileE.open("lg_eigen.txt");
	outputFileE << "// abscissae for lg_eigen" << std::endl;
	for(auto xx : x)
	{
		outputFileE << std::setw(2) << std::setprecision(25) << xx << std::endl;
	}

	outputFileE << "// weights for lg_eigen" << std::endl;

	for (auto ww : w)
	{
		outputFileE << std::setw(2) << std::setprecision(25) <<  ww << std::endl;
	}



	std::ofstream outputFile;
	outputFile.open("lg.txt");

	if (not outputFile.is_open())
	{
		std::cerr << "unable to open file." << std::endl;
		return EXIT_FAILURE;
	}


	for (unsigned int i = 2; i < 100; i += 2)
	{
		auto coeffs = GenerateLegendreAbscissae(i);

		auto weights = GenerateLegendreWeights(i, coeffs);

 		outputFile << "// coefficients for " << i << std::endl;
		for (auto& c : coeffs)
		{
			 
			
			if (c > 0)
			{
				outputFile << std::setw(2) << std::setprecision(25) << c << std::endl;
			}
		}
		outputFile << "// weights for n = " << i << std::endl;
	
		for (auto& w : weights)
		{


			if (w > 0)
			{
				outputFile << std::setw(2) << std::setprecision(25) << w << std::endl;
			}
		}
	}

}



