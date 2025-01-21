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

#include <iomanip>

#include <cstdlib>


#include "GenerateLegendreCoeffs.h"
#include "WriteVector.h"


int GenerateEigenLegendreConstants()
{
	try
	{
		std::ofstream outputFile;
		outputFile.open("EigenLegendreConstants.h");

		if (not outputFile.is_open())
		{
			std::cerr << "unable to open file." << std::endl;
			return EXIT_FAILURE;
		}


		outputFile << "#ifndef EIGEN_LEGENDRE_CONSTANTS_\n";
		outputFile << "#define EIGEN_LEGENDRE_CONSTANTS_\n";
		outputFile << "#include <map>\n";
		outputFile << "#include <array>\n";
		outputFile << "#include <vector>\n";
		outputFile << "namespace numerical {\n";


		for (unsigned int i = 2; i <= 64; i++)
		{

			Eigen::VectorXd x;
			Eigen::VectorXd w;

			ComputeAbscissasAndWeights(i, x, w);

			std::vector<double> stdW(w.data(), w.data() + w.size());
			std::vector<double> stdX(x.data(), x.data() + x.size());


			outputFile << "// abscissae for " << i << std::endl;

			outputFile << "const std::vector<double> " << " x" << i << " = {";
			for (size_t j = 0; j < stdX.size() - 1; ++j)
			{
				{
					outputFile << std::setw(2) << std::setprecision(25) << stdX[j] << ",";
				}
			}
			outputFile << std::setw(2) << std::setprecision(25) << stdX[stdX.size() - 1] << "};" << std::endl;


			outputFile << "// weights for " << i << std::endl;

			outputFile << "const std::vector<double> " << " w" << i << " = {";
			for (size_t j = 0; j < stdW.size() - 1; ++j)
			{
				{
					outputFile << std::setw(2) << std::setprecision(25) << stdW[j] << ",";
				}
			}
			outputFile << std::setw(2) << std::setprecision(25) << stdW[stdW.size() - 1] << "};" << std::endl <<
				std::endl;
		}


		outputFile << "const std::map<unsigned int, std::vector<double>> legendreAbscissae = {\n";

		for (unsigned int i = 2; i <= 64; i++)
		{
			outputFile << "    { " << i << ", x" << i << "},\n";
		}
		outputFile << "};" << std::endl << std::endl;


		outputFile << "const std::map<unsigned int, std::vector<double>> legendreWeights = {\n";
		for (unsigned int i = 2; i <= 64; i++)
		{
			outputFile << "    { " << i << ", w" << i << "},\n";
		}
		outputFile << "};" << std::endl;
		outputFile << "}" << std::endl;
		outputFile << "#endif" << std::endl;

		return EXIT_SUCCESS;
	}
	catch (std::exception& ex)
	{
		std::cerr << ex.what();
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}


int GenerateLegendreConstants()
{
	try
	{
		std::ofstream outputFile;
		outputFile.open("LegendreConstants.h");

		if (not outputFile.is_open())
		{
			std::cerr << "unable to open file." << std::endl;
			return EXIT_FAILURE;
		}


		outputFile << "#ifndef LEGENDRE_CONSTANTS_\n";
		outputFile << "#define LEGENDRE_CONSTANTS_\n";
		outputFile << "#include <map>\n";
		outputFile << "#include <array>\n";
		outputFile << "#include <vector>\n";
		outputFile << "namespace numerical {\n";


		for (unsigned int i = 2; i <= 64; i++)
		{
			std::vector<double> stdX = GenerateLegendreAbscissae(i);

			std::vector<double> stdW = GenerateLegendreWeights(i, stdX);


			outputFile << "// abscissae for " << i << std::endl;

			outputFile << "const std::vector<double> " << " x" << i << " = {";
			for (size_t j = 0; j < stdX.size() - 1; ++j)
			{
				{
					outputFile << std::setw(2) << std::setprecision(25) << stdX[j] << ",";
				}
			}
			outputFile << std::setw(2) << std::setprecision(25) << stdX[stdX.size() - 1] << "};" << std::endl;


			outputFile << "// weights for " << i << std::endl;

			outputFile << "const std::vector<double> " << " w" << i << " = {";
			for (size_t j = 0; j < stdW.size() - 1; ++j)
			{
				{
					outputFile << std::setw(2) << std::setprecision(25) << stdW[j] << ",";
				}
			}
			outputFile << std::setw(2) << std::setprecision(25) << stdW[stdW.size() - 1] << "};" << std::endl <<
				std::endl;
		}


		outputFile << "const std::map<unsigned int, std::vector<double>> legendreAbscissae = {\n";

		for (unsigned int i = 2; i <= 64; i++)
		{
			outputFile << "    { " << i << ", x" << i << "},\n";
		}
		outputFile << "};" << std::endl << std::endl;


		outputFile << "const std::map<unsigned int, std::vector<double>> legendreWeights = {\n";
		for (unsigned int i = 2; i <= 64; i++)
		{
			outputFile << "    { " << i << ", w" << i << "},\n";
		}
		outputFile << "};" << std::endl;
		outputFile << "}" << std::endl;
		outputFile << "#endif" << std::endl;

		return EXIT_SUCCESS;
	}
	catch (std::exception& ex)
	{
		std::cerr << ex.what();
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int main()
{
	int r = GenerateEigenLegendreConstants();

	if (r == EXIT_SUCCESS)
	{
		r = GenerateLegendreConstants();
	}

	return r;
}
