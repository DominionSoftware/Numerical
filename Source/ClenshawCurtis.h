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

#ifndef CLENSHAW_CURTIS_
#define CLENSHAW_CURTIS_

#include <ranges>
#include <algorithm>
#include <numeric>

namespace numerical::ClenshawCurtis
{
	enum class ErrorCode
	{
		NoError,
		MaxSubIntervalsReached,
		RoundoffError,
		ConvergenceNotAchieved,
		InvalidInput,
		BadIntegrandBehavior,
		UnableToIntegrate
	};

	template <std::floating_point T>
	struct Parameters
	{
		T absoluteTolerance = 1e-10;
		T relativeTolerance = 1e-10;
		size_t maxSubIntervals = 1000;
		T roundoffTolerance = 1e-30;
	};

	template <std::floating_point T>
	class Function
	{
	public:
		std::function<T(T)> integrand;
		void* additionalParameters;

		T operator()(T input) const
		{
			return integrand(input);
		}
	};

	template <std::floating_point T>
	struct IntegrationResult
	{
		T integral;
		T absoluteError;
		size_t evaluations;
		ErrorCode errorCode;
	};

	template <std::floating_point T>
	std::array<T, 12> ComputeChebyshevCoefficients12(const std::array<T, 25>& functionValues,
	                                                 const std::array<T, 12> differences)
	{
		std::array<T, 12> coefficients;

		coefficients[0] = std::reduce(std::next(functionValues.begin()), functionValues.begin() + 12,
		                              functionValues[0]) / 12.0;
		for (auto j : std::views::iota(1u, 12u))
		{
			coefficients[j] = 0.5 * functionValues[0];

			auto evenView = std::views::iota(1u, 12u) |
				std::views::transform([&](size_t i)
				{
					return functionValues[i] * std::cos(j * i * std::numbers::pi / 12);
				});
			coefficients[j] += std::accumulate(evenView.begin(), evenView.end(), T{0});

			if (j & 0x01)
			{
				auto oddView = std::views::iota(0u, 12u) |
					std::views::transform([&](size_t i)
					{
						return differences[i] * std::sin(j * (i + 0.5) * std::numbers::pi / 12);
					});
				coefficients[j] += std::accumulate(oddView.begin(), oddView.end(), T{0});
			}
			coefficients[j] *= 1.0 / 6.0;
		}
		return coefficients;
	}

	template <std::floating_point T>
	std::array<T, 24> ComputeChebyshevCoefficients24(
		const std::array<T, 25>& functionValues,
		const std::array<T, 12>& differences)
	{
		std::array<T, 24> coefficients{};

		coefficients[0] = std::reduce(
			std::next(functionValues.begin()),
			functionValues.begin() + 24,
			functionValues[0]
		) / 24.0;

		for (auto j : std::views::iota(1u, 24u))
		{
			coefficients[j] = 0.5 * functionValues[0];

			auto evenView = std::views::iota(1u, 24u) |
				std::views::transform([&](size_t i)
				{
					return functionValues[i] * std::cos(j * i * std::numbers::pi / 24);
				});
			coefficients[j] += std::accumulate(evenView.begin(), evenView.end(), T{0});

			if (j % 2 == 1)
			{
				// Odd coefficients
				auto oddView = std::views::iota(0u, 12u) |
					std::views::transform([&](size_t i)
					{
						return differences[i] * std::sin(j * (i + 0.5) * std::numbers::pi / 24);
					});
				coefficients[j] += std::accumulate(oddView.begin(), oddView.end(), T{0});
			}
			coefficients[j] *= 1.0 / 12.0;
		}
		return coefficients;
	}


	template <std::floating_point T>
	std::tuple<std::array<T, 12>, std::array<T, 24>> GenerateChebyshev(const Function<T>& function, T lowerBound,
	                                                                   T upperBound)
	{
		constexpr std::array<T, 11> chebyshevNodes = {
			0.9914448613738104, 0.9659258262890683, 0.9238795325112868,
			0.8660254037844386, 0.7933533402912352, 0.7071067811865475,
			0.6087614290087206, 0.5000000000000000, 0.3826834323650898,
			0.2588190451025208, 0.1305261922200516
		};

		std::array<T, 25> functionValues = {};
		std::array<T, 12> differences = {};
		std::array<T, 25> sumValues{};
		const T intervalCenter = 0.5 * (upperBound + lowerBound);
		const T intervalHalfLength = 0.5 * (upperBound - lowerBound);
		functionValues[0] =  function(upperBound);
		functionValues[12] = function(intervalCenter);
		functionValues[24] = function(lowerBound);


		for (auto i : std::views::iota(1u, 12u)) {
			const size_t oppositeIndex = 24 - i;
			const T nodeOffset = intervalHalfLength * chebyshevNodes[i - 1];
			const T xPlus = intervalCenter + nodeOffset;
			const T xMinus = intervalCenter - nodeOffset;
			functionValues[i] = function(xPlus);
			functionValues[oppositeIndex] = function(xMinus);
		}

		for (auto i : std::views::iota(0u, 12u)) {
			const size_t oppositeIndex = 24 - i;
			differences[i] = functionValues[i] - functionValues[oppositeIndex];
			functionValues[i] = functionValues[i] + functionValues[oppositeIndex];
		}
		return {
			ComputeChebyshevCoefficients12(functionValues, differences),
			ComputeChebyshevCoefficients24(functionValues, differences)
		};
	}


	template <std::floating_point T>
	IntegrationResult<T> Integrate(const Function<T>& function, T lowerBound, T upperBound,
	                               const Parameters<T>& parameters = Parameters<T>())
	{
		using SubInterval = std::vector<std::tuple<T, T>>;

		if (lowerBound > upperBound)
		{
			return IntegrationResult<T>{0, 0, 0, ErrorCode::InvalidInput};
		}


		SubInterval subIntervals;
		size_t numEvaluations{};
		T totalIntegral{};
		T totalError{};
		subIntervals.push_back({lowerBound, upperBound});
		while (!subIntervals.empty() && numEvaluations < parameters.maxSubIntervals)
		{
			auto [lower,upper] = subIntervals.back();
			subIntervals.pop_back();

			auto [chebychev12, chebychev24] = GenerateChebyshev<T>(function, lower, upper);

			numEvaluations += 25;

			T integral24 = chebychev24[0];
			T integral12 = chebychev12[0];

			for (size_t k = 2; k < 24; k += 2) {
				integral24 -= chebychev24[k] / (k * k - 1);
			}
			integral24 *= 2;
			for (size_t k = 2; k < 12; k += 2) {
				integral12 -= chebychev12[k] / (k * k - 1);
			}
			integral12 *= 2;
			integral24 *= (upper - lower);
			integral12 *= (upper - lower);

			T error = std::abs(integral24 - integral12);


			T tolerance = std::max(parameters.absoluteTolerance, parameters.relativeTolerance * std::abs(integral24));


			if (error < tolerance)
			{
				totalIntegral += integral24;
				totalError += error;
			}
			else if (numEvaluations >= parameters.maxSubIntervals)
			{
				return IntegrationResult<T>(totalIntegral, totalError, numEvaluations,
				                            ErrorCode::MaxSubIntervalsReached);
			}
			else
			{
				T mid = (lower + upper) * 0.5;
				subIntervals.push_back({mid, upper});
				subIntervals.push_back({lower, mid});
			}
		}
		if (std::abs(totalIntegral) < parameters.roundoffTolerance)
		{
			return IntegrationResult<T>(0, totalError, numEvaluations, ErrorCode::RoundoffError);
		}

		return IntegrationResult<T>(totalIntegral, totalError, numEvaluations, ErrorCode::NoError);
	}
}
#endif // !CHEBYSHEV12_24_
