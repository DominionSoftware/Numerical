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

namespace numerical
{
	
	template<std::floating_point T>
	class ClenshawCurtis
	{
	public:
		explicit ClenshawCurtis(int numberOfPoints) : numberOfPoints_(numberOfPoints)
		{
			generateNodes(numberOfPoints);
			generateWeights(numberOfPoints);

		}

		T chebyshev(int degree, T x)
		{
			if (degree == 0)
			{
				return 1.0;
			}
			if (degree == 1)
			{
				return x;
			}

			T previous = 1.0;
			T current = x;
			T next{};

			for (int i = 2; i <= degree; i++) {
				next = 2.0 * x * current - previous;
				previous = current;
				current = next;
			}
			return current;
		}

		void generateNodes(int numberOfPoints)
		{
			nodes_.resize(numberOfPoints);
			for (int nodeIndex = 0; nodeIndex < numberOfPoints; nodeIndex++)
			{
				nodes_[nodeIndex] = -std::cos(nodeIndex * std::numbers::pi / (numberOfPoints - 1));
			}

		}

		void generateWeights(int numberOfPoints)
		{

			weights_.resize(numberOfPoints);
			for (int j = 0; j < numberOfPoints; j++) 
			{
				weights_[j] = 1.0;
				for (int k = 1; k <= (numberOfPoints - 1) / 2; k++) 
				{
					weights_[j] -= 2.0 * std::cos(2.0 * k * j * std::numbers::pi / (numberOfPoints - 1))
						/ (4 * k * k - 1);
				}
			}

			weights_[0] /= (numberOfPoints - 1);
			weights_[numberOfPoints - 1] /= (numberOfPoints - 1);

			for (int j = 1; j < numberOfPoints - 1; j++) 
			{
				weights_[j] *= 2.0 / (numberOfPoints - 1);
			}
		}


		T Integrate(std::function<T(T)>function, T lowerBound, T upperBound) const

		{

			T integralSum = 0.0;
			T intervalScale = (upperBound - lowerBound) / 2.0;
			T intervalMidpoint = (upperBound + lowerBound) / 2.0;

			for (int nodeIndex = 0; nodeIndex < numberOfPoints_; nodeIndex++) {
				double transformedPoint = intervalScale * nodes_[nodeIndex] + intervalMidpoint;
				integralSum += weights_[nodeIndex] * function(transformedPoint);
			}

			return intervalScale * integralSum;
		}

	private:
		std::vector<T> nodes_;
		std::vector<T> weights_;
		int numberOfPoints_;
	};
}
#endif 
