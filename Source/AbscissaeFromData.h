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

#ifndef ABSCISSAE_FROMDATA_
#define ABSCISSAE_FROMDATA_
#include "Abscissae.h"
#include <vector>

namespace numerical
{


	/**
	 * @brief Abscissae representing values obtained from existing data.
	 *
	 * This class represents the abscissae of values derived from existing data,
	 * where each value in the abscissae corresponds to a value from the provided data.
	 *
	 * @tparam T The type of values in the abscissae.
	 */

	template <typename T>
	requires std::floating_point<T>
	class AbscissaeFromData : public Abscissae<T>
	{
	public:

		AbscissaeFromData() = default;

		 
		/**
		 * @brief Constructs the abscissae from existing data.
		 *
		 * Constructs the abscissae using the provided vector of data.
		 *
		 * @param data The vector containing the data to be used for the abscissae.
		 */
		explicit AbscissaeFromData(const std::vector<T>& data) : data_(data)
		{

		}

		/**
		* @brief Get the value at the specified index in the abscissae (const version).
		*
		* @param i The index of the value to retrieve.
		* @return The value at the specified index.
		*/
		T operator [](size_t i) const override
		{
			return data_[i];
		}

		/**
		* @brief Get a reference to the value at the specified index in the abscissae.
		*
		* @param i The index of the value to retrieve.
		* @return A reference to the value at the specified index.
		*/
		T& operator [](size_t i) override { return data_[i]; }

		/**
		* @brief Get the value at the specified index in the abscissae.
		*
		* This method retrieves the value at the specified index in the abscissae.
		*
		* @param i The index of the value to retrieve.
		* @return The value at the specified index.
		*/
		double at(size_t i) override
		{
			return data_.at(i);
		}

		/**
		* @brief Get the size of the abscissae.
		*
		* This method returns the number of values in the abscissae.
		*
		* @return The size of the abscissae.
		*/
		size_t size() override
		{
			return data_.size();
		}

		std::vector<T> data()
		{
			return data_;
		}
	private:

		std::vector<T> data_;


	};
}


#endif
