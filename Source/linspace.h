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

#ifndef LINSPACE_
#define LINSPACE_

#include <vector>
#include <algorithm>
#include <numeric>


template<std::floating_point T>
bool isEqualFloat(T a, T b)
{
	return std::abs(b - a) < 0.001;
}

namespace numerical
{


	/**
	 * @brief Generate a vector containing evenly spaced values within a specified range.
	 *
	 * This function generates `n` evenly spaced values between `a` and `b`, inclusive.
	 *
	 * @tparam T      The type of the values in the vector.
	 * @param a       The starting value of the range.
	 * @param b       The ending value of the range.
	 * @param n       The number of values to generate.
	 * @return        A vector containing `n` evenly spaced values between `a` and `b`, inclusive.
	 *
	 * @note          If `n` is less than 2, the returned vector will contain only `b`.
	 */
	template<typename T>
	requires std::floating_point<T> || std::integral<T>
	std::vector<T> Linspace(T a, T b, size_t n)
	{
		bool endsEqual = false;

		std::vector<T> temp1;

		if (std::is_floating_point<T>())
		{
			if (isEqualFloat(a, -b))
			{
				endsEqual = true;
			}
		}

		if (std::is_integral<T>())
		{
			if (a == -b)
			{
				endsEqual = true;
			}
		}
		if (endsEqual)
		{
			size_t l = n - 1;
			temp1.resize(n - 1);


			std::iota(temp1.begin(), temp1.end(), 0);

			T div = (b - a) / static_cast<T>(l);

			std::for_each(std::begin(temp1), std::end(temp1), [=](T& v)
				{
					v = a + v * div;
				});
			temp1.push_back(b);
			temp1[0] = a;
			temp1[temp1.size() - 1] = b;

			if (n & 1)
			{
				temp1[n / 2 + 1] = 0;
			}

			return temp1;
		}
		else
		{
			size_t l = n - 1;
			T v = a;
			for (size_t i = 0; i <= l; ++i)
			{
				v = a + static_cast<T>(i) * (b - a) / l;
				temp1.push_back(v);

			}
		}
		return temp1;
	}
}

#endif

