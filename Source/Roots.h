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

#ifndef ROOTS_
#define ROOTS_
#include <functional>
#include <vector>
#include <concepts>
#include <fstream>
#include <optional>
#include <cmath>

namespace numerical
{



	// initial guess, e.g.,  cos(std::numbers::pi * (i + 0.75) / (n + 0.5));
	template<std::floating_point T = double>
	std::vector<T> NewtonRhapson(unsigned int n, std::function<T(unsigned int, unsigned int)> initialGuess, std::function<T(const unsigned int, const T)> f, std::function<T(const unsigned int, const T)> fPrime, double tol = 1e-10)
	{
		std::vector<double> roots;

		double x;
		double x1;

		for (unsigned int i = 0; i < n; ++i) {


			x = initialGuess(i, n);

			do {
				x1 = x;
				x = x1 - f(n, x1) / fPrime(n, x1);
			} while (std::fabs(x - x1) > tol);

			roots.push_back(x);
		}

		return roots;

	}

	template<std::floating_point T = double>
	T NewtonRhapson(unsigned int n,T initialGuess, std::function<T(const T)> f, std::function<T( const T)> fPrime, double tol = 1e-10)
	{

		double x = initialGuess;
		double x1;
		unsigned int limit = n;
		do {
			x1 = x;
			x = x1 - f(x1) / fPrime (x1);
		} while (fabs(x - x1) > tol && --limit > 0);


		return x;

	}


	enum class FZeroErrors
	{
		NoError,
		NotBracketed,
		FileError
	};

	template<typename T>
	bool equalZero(T n1, T eps = std::numeric_limits<T>::epsilon())
	{
		return std::fabs(n1 - 0.0) < eps;
	};

	template<typename T>
	T safeDivide(T numerator, T denominator, T eps = std::numeric_limits<T>::epsilon())
	{
		if (equalZero(denominator, eps))
		{
			return 0.0;
		}

		return numerator / denominator;
	}


	template<typename T>
	std::tuple<FZeroErrors, T> BrentZero(const std::function<T(T)>& f, T a, T b)
	{
		T Ya = f(a);
		T Yb = f(b);
		T Yc = Ya;

		if (std::signbit(Ya) == std::signbit(Yb))
		{
			return std::make_tuple(FZeroErrors::NotBracketed, std::nanf("NaN"));
		}

		T c = a;
		T d = b - c;
		T e = d;

		while (!equalZero<T>(Yb))
		{
			bool aEQc{ false };

			if (std::signbit(Ya) == std::signbit(Yb))
			{
				a = c;
				Ya = Yc;
				d = b - c;
				e = d;
				aEQc = true;
			}
			if (std::fabs(Ya) < std::fabs(Yb))
			{

				c = b;
				b = a;
				a = c;

				Yc = Yb;
				Yb = Ya;
				Ya = Yc;
				aEQc = true;

			}

			T m = static_cast<T>(0.5) * (a - b);
			T tol = static_cast<T>(2.0) * std::numeric_limits<T>::epsilon() * std::max<T>(std::fabs(b), static_cast<T>(1.0));

			if (std::fabs(m) <= tol || equalZero<T>(Yb))
			{
				break;
			}

			if (std::fabs(e) < tol || std::fabs(Yc) <= std::fabs(Yb))
			{
				// Bisection 
				d = m;
				e = m;
			}
			else
			{

				T s = Yb / Yc;
				T p{ 0.0 };
				T q{ 0.0 };
				T r{ 0.0 };

				if (aEQc)
				{
					// Linear, Secant
					p = static_cast<T>(2.0) * m * s;
					q = static_cast<T>(1.0) - s;

				}
				else
				{
					// Inverse Quadratic = 

					q = Yc / Ya;
					r = Yb / Ya;
					p = s * static_cast<T>(2.0) * m * q * (q - r) - (b - c) * (r - static_cast<T>(1.0));
					q = (q - static_cast<T>(1.0)) * (r - static_cast<T>(1.0)) * (s - static_cast<T>(1.0));

				}
				if (p > 0.0)
				{
					q = -q;
				}
				else
				{
					p = -p;
				}

				if (static_cast<T>(2.0) * p < static_cast<T>(3.0) * m * q - std::fabs(tol * q) && p < std::fabs(static_cast<T>(0.5) * e * q))
				{
					// Interpolated Point acceptable
					e = d;
					d = p / q;
				}
				else
				{
					// Interpolated Point Not acceptable 
					d = m;
					e = m;
				}
			}

			c = b;
			Yc = Yb;
			if (std::fabs(d) > tol)
			{
				b = b + d;
			}
			else
			{
				T sg = std::signbit(b - a) ? static_cast<T>(1.0) : static_cast<T>(-1.0);
				b = b - sg * tol;
			}
			Yb = f(b);
		}

		return std::make_tuple(FZeroErrors::NoError, b);

	}

}


#endif


