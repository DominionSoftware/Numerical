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

#pragma once
#include <functional>
#include <iostream>
#include <vector>
#include "Abscissae.h"



#pragma optimize("",off)

namespace numerical
{


	/**
	 * @brief A class for computing finite difference approximations of derivatives.
	 *
	 * This class provides methods for computing first-order derivative approximations
	 * using finite difference methods.
	 *
	 * @tparam T The type of values for the function and abscissae.
	 */
	template<std::floating_point T, std::floating_point D>
	class FiniteDifference
	{
	public:
		/**
		 * @brief Compute the first-order derivatives of a function over abscissae using finite differences.
		 *
		 * This method computes the first-order derivatives of a given function over a abscissae
		 * using finite difference methods.
		 *
		 * @param f The function for which derivatives are computed.
		 * @param abscissae The abscissae over which the function is defined.
		 * @return A vector containing the first-order derivatives of the function at each point in the abscissae.
		 */
		std::vector<T> gradient(std::function<T(D)> f, Abscissae<D>* abscissae)
		{
			std::vector<T> result;
			for (size_t index = 0; index < abscissae->size() - 1; index++)
			{
				if (index == 0)
				{
					result.push_back(forwardDifferenceAtFirstPoint(f, abscissae));
				}
				else
				{
					result.push_back(centeredDifference(f, abscissae, index));
				}
			}

			result.push_back(backwardDifferenceAtLastPoint(f, abscissae));

			return result;
		}

	private:

		T forwardDifferenceAtFirstPoint(std::function<T(D)> f, Abscissae<D>* abscissae)
		{
			T h = (abscissae->at(1) - abscissae->at(0));
			T x = abscissae->at(0);

			T n = -3 * f(x) + (4 *  f(x + h)) - f(x + 2 * h);
		
			return n / (2 * h);
		}


		T backwardDifferenceAtLastPoint(std::function<T(D)> f, Abscissae<D>* abscissae)
		{
			size_t sz = abscissae->size();
			T h = (abscissae->at(sz - 1) - abscissae->at(sz - 2));
			T x = abscissae->at(sz - 1);
			T n = 3 * f(x) - (4 * f(x-h)) + f(x -2 * h);
			return n / (2 * h);
		}



		/**
		* @brief Compute the forward difference of a function at a given index.
		*
		* @param f The function for which the forward difference is computed.
		* @param abscissae The abscissae over which the function is defined.
		* @param index The index at which to compute the forward difference.
		* @return The forward difference of the function at the specified index.
		*/
		T forwardDifference(std::function<T(D)> f, Abscissae<D>* abscissae, size_t index)
		{
			T n = f(abscissae->at(index + 1)) - f(abscissae->at(index));
			T d = (abscissae->at(index + 1) - abscissae->at(index));
			return n / d;
		}

		/**
		* @brief Compute the backward difference of a function at a given index.
		*
		* @param f The function for which the backward difference is computed.
		* @param abscissae The abscissae over which the function is defined.
		* @param index The index at which to compute the backward difference.
		* @return The backward difference of the function at the specified index.
		*/
		T backwardDifference(std::function<T(D)> f, Abscissae<D>* abscissae, size_t index)
		{
			T n = f(abscissae->at(index)) - f(abscissae->at(index - 1));
			T d = (abscissae->at(index) - abscissae->at(index - 1));
			return n / d;
		}
		/**
		* @brief Compute the centered difference of a function at a given index.
		*
		* @param f The function for which the centered difference is computed.
		* @param abscissae The abscissae over which the function is defined.
		* @param index The index at which to compute the centered difference.W
		* @return The centered difference of the function at the specified index.
		*/
		T centeredDifference(std::function<T(D)> f, Abscissae<D>* abscissae, size_t index)
		{
			T n = f(abscissae->at(index + 1)) - f(abscissae->at(index - 1));
			T h = (abscissae->at(index + 1) - abscissae->at(index - 1));
			return n / h;
		}

	};
}

