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

#ifndef INTEGRAL_
#define INTEGRAL_
#include <cassert>
#include <numeric>

#include "Abscissae.h"
#include "FiniteDifference.h"
#ifdef EIGENLEGENDRE
#include "EigenLegendreConstants.h"
#else
#include "LegendreConstants.h"
#endif

#include "GaussKronrodRules.h"
#include "ClenshawCurtis.h"

namespace numerical
{

	/**
	 * @brief Performs numerical integration using the trapezoidal rule.
	 *
	 * @tparam T The data type of the abscissae.
	 * @tparam R The return type of the function to be integrated.
	 * @param abscissae A pointer to the abscissae.
	 * @param f The function to be integrated.
	 * @return The result of the integration.
	 *
	 * This function calculates the definite integral of a function `f` over a specified interval.
	 * The interval is defined by the `abscissae` parameter. The function uses the trapezoidal rule for numerical integration.
	 * The trapezoidal rule works by approximating the region under the graph of the function as a trapezoid and calculating its area.
	 */

	template<typename T, typename R>
	requires std::floating_point<T>

	T IntegrateTrapezoid(Abscissae<T>* abscissae, std::function<R(T)> f)
	{
		auto lower = abscissae->at(0);
		auto upper = abscissae->at(abscissae->size() - 1);

		double step = (upper - lower) / static_cast<T>(abscissae->size() - 1);

		auto sum = 0.5 * (f(upper) - f(lower));
		for (size_t i = 1; i < abscissae->size(); ++i)
		{
			double x = abscissae->at(i);
			sum += f(x);
		}

		return  step * sum;
	}

	/**
	 * @brief Computes the definite integral of a function using Simpson's rule.
	 *
	 * @tparam T The type of the abscissae, must be a floating point type.
	 * @tparam R The type of the return value, must be a floating point type.
	 * @param abscissae Pointer to the abscissae, which must have an even count and more than 3 points.
	 * @param f The function to integrate.
	 *
	 * @return The definite integral of the function from the first to the last point in the abscissae.
	 *
	 * @throws std::runtime_error If the abscissae has an odd count or less than 4 points.
	 */

	template<typename T, typename R>
	requires std::floating_point<T>
	T IntegrateSimpsons(Abscissae<T>* abscissae, std::function<R(T)> f)
	{


		if (((abscissae->size() & 1) == 1) || (abscissae->size() <= 3))
		{
			throw std::runtime_error("abscissae must have even count and more than 3 points");
		}


		auto lower = abscissae->at(0);
		auto upper = abscissae->at(abscissae->size() - 1);

		const double sum = f(lower) + f(upper);
		double odd{ 0 };
		double even{ 0 };


		for (size_t i = 1; i < abscissae->size() - 1; ++i)
		{
			double x = abscissae->at(i);
			if (i & 0x01)
			{
				odd += f(x);
			}
			else
			{
				even += f(x);
			}
		}

		const double step = (upper - lower) / (abscissae->size() - 1);

		return ((step / 3) * (sum + (4 * odd) + (2 * even)));
	}


	/**
	 * @brief Performs numerical integration using Simpson's rule.
	 *
	 * This function calculates the definite integral of a function over a specified interval using Simpson's rule.
	 * The function to be integrated is passed as a parameter.
	 *
	 * @tparam T The type of the lower and upper bounds of the interval. Must be a floating point type.
	 * @tparam R The type of the return value. Must be a floating point type.
	 *
	 * @param lower The lower bound of the interval over which to integrate.
	 * @param upper The upper bound of the interval over which to integrate.
	 * @param count The number of subintervals to use in the approximation. Must be an even number greater than 3.
	 * @param f The function to integrate. This is a function that takes a single argument of type T and returns a value of type R.
	 *
	 * @return The approximate value of the definite integral of the function over the specified interval.
	 *
	 * @throws std::runtime_error If count is not an even number or is less than or equal to 3.
	 *
	 * @note The accuracy of the approximation increases with the number of subintervals (count).
	 */

	template<typename T, typename R>
	requires std::floating_point<T> || std::floating_point<R>
	T IntegrateSimpsons(T lower, T upper, size_t count, std::function<R(T)> f)
	{
		static_assert(std::is_floating_point_v<T>, "Type must be real value");
		static_assert(std::is_floating_point_v<R>, "Type must be real value");


		if (((count & 1) == 1) || (count <= 3))
		{
			throw std::runtime_error("abscissae must have even count and more than 3 points");
		}

		const double sum = f(lower) + f(upper);
		double odd{ 0 };
		double even{ 0 };
		double h = (upper - lower) / (count - 1);

		for (size_t i = 1; i < count - 1; ++i)
		{
			double x = lower + (i * h);

			if (i & 0x01)
			{
				odd += f(x);

			}
			else
			{
				even += f(x);
			}
		}


		return { (h / 3) * (sum + (4 * odd) + (2 * even)) };
	}


	// Utilities for Legendre Polynomials.

	// Function to calculate the derivative for Gauss-Legendre quadrature
	template<typename T>
	requires std::floating_point<T>
	T legendrePolynomialDerivative(unsigned int n, T x)
	{
		return n * (x * std::legendre(n, x) - std::legendre(n - 1, x)) / (x * x - 1.0);
	}


	// Function to calculate the weights for Gauss-Legendre quadrature
	template<typename T = double>
	requires std::floating_point<T>
	std::vector<T> calculateWeights(unsigned int n, const std::vector<T>& roots)
	{
		std::vector<T> weights(n);

		for (int i = 0; i < n; ++i) {
			weights[i] = 2.0 / ((1.0 - roots[i] * roots[i]) * legendrePolynomialDerivative(n, roots[i]) * legendrePolynomialDerivative(n, roots[i]));
		}

		return weights;
	}

	/**
	 * @brief Performs numerical integration using the Gauss-Legendre method.
	 *
	 * This function calculates the definite integral of a given function within a specified interval [a, b] using the Gauss-Legendre method.
	 * The Gauss-Legendre method approximates the definite integral by summing the weighted function values at specific points within the interval.
	 *
	 * @param a The lower limit of the interval.
	 * @param b The upper limit of the interval.
	 * @param n The number of points to be used in the Gauss-Legendre method. Must be less than or equal to 64.
	 * @param func The function to be integrated. This function should take a double as input and return a double.
	 *
	 * @return The approximate value of the definite integral of the function within the interval [a, b].
	 *
	 * @note The function `func` should be continuous and differentiable within the interval [a, b] for accurate results.
	 * @note The accuracy of the approximation increases with the number of points `n`.
	 */

	inline double IntegrateGaussLegendre(double a, double b, unsigned int n, const std::function<double(double)>& func)
	{
		assert(n <= 64);

		const double c = (b - a);
		const double d = (b + a);

		double sum{ 0 };

		for (unsigned i = 0; i < n; ++i)
		{
			const double w = numerical::legendreWeights.at(n)[i];
			const double x = numerical::legendreAbscissae.at(n)[i];
			sum += w * func((c * x + d) / 2.0) * c / 2.0;
		}

		return sum;
	}


	inline std::tuple<numerical::GuassKronrodErrorCode, const char *, double>   IntegrateGaussKronrod(numerical::GuassKronrodParameters<double>& parameters,
		const std::function<double(double)>& func,
		void (*rule)(const std::function<double(double)>&, double, double, double&, double&, double&, double&), const double epsilonAbs = 0.0,const double epsRelative = 1.0e-7)
	{
		int roundoffType1 = 0;
		int roundoffType2 = 0;
		int errorType = 0;

		double result{ 0 };
		double absoluteError{ 0 };
		double defabs{ 0 };
		double resultAbsolute{ 0 };

		rule(func, parameters.getInitialLowerBounds(), parameters.getInitialUpperBounds(), result, absoluteError, defabs, resultAbsolute);

		parameters.setFirstValues(result, absoluteError);

		double errorBound = std::max(epsilonAbs, epsRelative * std::fabs(result));

		const double roundOff = (50 * std::numeric_limits<double>::epsilon() * resultAbsolute);

		if (absoluteError <= roundOff && absoluteError > errorBound)
		{
			return std::make_tuple(numerical::GuassKronrodErrorCode::RoundoffError,"Cannot reach tolerance because of round-off error.", result);  
		}

		if ((absoluteError <= errorBound && absoluteError != resultAbsolute) || absoluteError == 0.0)
		{

			return std::make_tuple(numerical::GuassKronrodErrorCode::NoError,"Normal Result.", result);  
		}

		if (parameters.limit() == 1)
		{
			return std::make_tuple(numerical::GuassKronrodErrorCode::MaxSubIntervalsReached, "Max sub-intervals reached.",result);  
		}


		double area = result;
		double errorSum = absoluteError;

		size_t iteration = 1;
		do
		{
			double aCurrent = 0;
			double bCurrent = 0;
			double resultCurrent = 0;
			double errorCurrent = 0;
			double area1 = 0;
			double area2 = 0;
			double error1 = 0;
			double error2 = 0;
			double resasc1 = 0;
			double resasc2 = 0;
			double resabs1 = 0;
			double resabs2 = 0;

			parameters.getValuesAtMaxError(aCurrent, bCurrent, resultCurrent, errorCurrent);

			const double a1 = aCurrent;
			const double b1 = 0.5 * (aCurrent + bCurrent);
			const double a2 = b1;
			const double b2 = bCurrent;

			rule(func, a1, b1, area1, error1, resabs1, resasc1);

			rule(func, a2, b2, area2, error2, resabs2, resasc2);

			const double area12 = area1 + area2;
			const double error12 = error1 + error2;
			errorSum += (error12 - errorCurrent);
			area += area12 - resultCurrent;

			if (resasc1 != error1 && resasc2 != error2)
			{
				const double deltaArea = resultCurrent - area12;

				if (fabs(deltaArea) <= 1.0e-5 * fabs(area12) && error12 >= 0.99 * errorCurrent)
				{
					roundoffType1++;
				}
				if (iteration >= 10 && error12 > errorCurrent)
				{
					roundoffType2++;
				}
			}
			errorBound = std::max(epsilonAbs, epsRelative * fabs(area));
			if (errorSum > errorBound)
			{
				if (roundoffType1 >= 6 || roundoffType2 >= 20)
				{
					errorType = 2;   /* round off error */
				}

				/* set error flag in the case of bad integrand behaviour at
				   a point of the integration range */
				double tmp = (1 + 100 * std::numeric_limits<double>::epsilon()) * (fabs(a2) + 1000 * std::numeric_limits<double>::min());

				if (const bool tooSmall = fabs(a1) <= tmp && fabs(b2) <= tmp)
				{
					errorType = 3;
				}
			}
			parameters.appendToLists(a1, b1, area1, error1, a2, b2, area2, error2);

			parameters.getValuesAtMaxError(aCurrent, bCurrent, resultCurrent, errorCurrent);

			iteration++;

		} while (iteration < parameters.limit() && !errorType && errorSum > errorBound);

		double resultSum = parameters.sumResults();

		if (errorSum <= errorBound)
		{
			return std::make_tuple(numerical::GuassKronrodErrorCode::NoError, "Normal Result",resultSum);  
			
		}
		else if (errorType == 2)
		{

			return std::make_tuple(numerical::GuassKronrodErrorCode::RoundoffError, "Round-off error prevents tolerance from being achieved",resultSum);  
		}
		else if (errorType == 3)
		{

			return std::make_tuple(numerical::GuassKronrodErrorCode::BadIntegrandBehavior, "Bad integrand behavior found in the integration interval.",resultSum); 
		}
		else if (iteration == parameters.limit())
		{
			return std::make_tuple(numerical::GuassKronrodErrorCode::MaxSubIntervalsReached, "Max subintevals reached.", result);  
		}
		
		return std::make_tuple(numerical::GuassKronrodErrorCode::UnableToIntegrate,"Unable to Integrate.", resultSum); 

	}


	inline std::tuple<numerical::ClenshawCurtis::ErrorCode, const char*, double>   IntegrateClenshawCurtis(numerical::ClenshawCurtis::Parameters<double>& parameters)
	{
		return std::make_tuple(numerical::ClenshawCurtis::ErrorCode::NoError, "NoError.", 0.0);

	}


}
#endif