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

#ifndef ABSCISSAE_
#define ABSCISSAE_
/**
 * @brief Abstract base class for representing a abscissae of values.
 *
 * This class serves as an abstract base class for representing a abscissae of values,
 * typically used in numerical computations or mathematical modeling.
 *
 * @tparam T The type of values in the abscissae.
 */
template<typename T>
requires std::floating_point<T>
class Abscissae
{
public:
	virtual ~Abscissae() {}

	/**
	* @brief Get the value at the specified index in the abscissae.
	*
	* This method retrieves the value at the specified index in the abscissae.
	*
	* @param i The index of the value to retrieve.
	* @return The value at the specified index.
	*/
	virtual T at(size_t i) = 0;

	/**
	 * @brief Get the size of the abscissae.
	 *
	 * This method returns the number of values in the abscissae.
	 *
	 * @return The size of the abscissae.
	 */
	virtual size_t size() = 0;

	/**
	 * @brief Get the value at the specified index in the abscissae (const version).
	 *
	 * This operator allows accessing the value at the specified index in the abscissae
	 * in a read-only manner.
	 *
	 * @param i The index of the value to retrieve.
	 * @return The value at the specified index.
	 */
	virtual T operator [](size_t i) const = 0;
	 
	/**
	* @brief Get a reference to the value at the specified index in the abscissae.
	*
	* This operator allows accessing and modifying the value at the specified index
	* in the abscissae.
	*
	* @param i The index of the value to retrieve.
	* @return A reference to the value at the specified index.
	*/
	virtual T& operator [](size_t i) = 0;

};


#endif


