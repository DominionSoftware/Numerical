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

#ifndef GENERATE_LEGENDRE_COEFFS_
#define GENERATE_LEGENDRE_COEFFS_


// Warnings triggered by Eigen
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4127) 
#pragma warning(disable: 5054) 
#endif

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wtype-limits"  // MSVC 4127 (conditional expression is constant)
#pragma GCC diagnostic ignored "-Wenum-compare" // MSVC 5054 (enum type comparison)
#endif




#include <algorithm>
#include <numbers>
#include <vector>
#include <Eigen/Dense>
#include "Integral.h"
#include "Roots.h"




template<std::floating_point T = double>
std::vector<T> GenerateLegendreWeights(unsigned int n, const std::vector<T>& roots)
{

	std::vector<T> weights(n);

	for (unsigned int i = 0; i < n; ++i)
	{
		weights[i] = (2.0 / ((1.0 - roots[i] * roots[i]) * numerical::legendrePolynomialDerivative(n, roots[i]) * numerical::legendrePolynomialDerivative(n, roots[i])));
	}

	return weights;
}



template<typename T = double>
requires std::is_floating_point_v<T>
std::vector<T> GenerateLegendreAbscissae(unsigned int n)
{


	auto l = [](const unsigned int nn, const T x)->T
		{
			return std::legendre(nn, x);
		};


	auto ld = [](const unsigned int nn, const T x)->T
		{
			return numerical::legendrePolynomialDerivative(nn, x);
		};


	auto guess = [](int i, int n)-> T 
	{
		return cos(std::numbers::pi * (i + 0.75) / (n + 0.5));
	};


	return numerical::NewtonRhapson<T>(n, guess, l, ld);

	

}



// Ported from https://www.mathworks.com/matlabcentral/fileexchange/26737-legendre-laguerre-and-hermite-gauss-quadrature
// Version 1.2.0.0 (3.25 KB) by Geert Van Damme

void ComputeAbscissasAndWeights(unsigned int n, Eigen::VectorXd& x, Eigen::VectorXd& w)
{
	// Create vector i = 1:n-1

	 auto i = Eigen::VectorXd::LinSpaced(n - 1, 1, n - 1);

	// Calculate vector a = i./sqrt(4*i.^2-1)
	Eigen::VectorXd a = i.array() / ((4 * i.array().square() - 1).sqrt());

	// Create matrix CM = diag(a, 1) + diag(a, -1)
	Eigen::MatrixXd CM = Eigen::MatrixXd::Zero(n, n);
	CM.diagonal(1) = a;
	CM.diagonal(-1) = a;

	// Calculate eigenvalues (L) and eigenvectors (V) of CM
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(CM);
	Eigen::VectorXd L = solver.eigenvalues();
	Eigen::MatrixXd V = solver.eigenvectors();

	// Sort the eigenvalues and rearrange the corresponding eigenvectors
	std::vector<int> indices(n);
	std::iota(indices.begin(), indices.end(), 0);
	std::sort(indices.begin(), indices.end(), [&L](int a, int b) { return L(a) < L(b); });

	x.resize(n);
	w.resize(n);

	for (unsigned int j = 0; j < n; ++j) {
		x(j) = L(indices[j]);
		auto temp = V.col(indices[j])(0);
		w(j) = 2 * (temp * temp);
	}
}


#ifdef _MSC_VER
#pragma warning(pop)
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif


#endif
