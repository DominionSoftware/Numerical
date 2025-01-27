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

#ifndef WRITEVECTOR_
#define WRITEVECTOR_
#include <fstream>
#include <string>

namespace numerical
{

	template<typename T>
	class FileIO
	{
	public:
		enum FileWriteResult
		{
			FileWriteOK,
			FileWriteErr
		};


		[[nodiscard]] static FileWriteResult WriteVector(std::vector<T>& v, const std::string& path)
		{
			std::ofstream output(path);




			if (!output.is_open())
			{
				return FileWriteResult::FileWriteErr;
			}

			for (size_t i = 0; i < v.size() - 1; i++)
			{
				output << v[i] << ",";
			}
			output << v[v.size() - 1] << std::endl;

			return FileWriteResult::FileWriteOK;

		}


		[[nodiscard]] static FileWriteResult WriteVector(std::ofstream& ofStream, std::vector<T>::const_iterator start, std::vector<T>::const_iterator end)
		{

			if (!ofStream.is_open())
			{
				return FileWriteResult::FileWriteErr;
			}

			for (auto iter = start; iter != end; iter++)
			{
				ofStream << *iter << ",";
			}
			ofStream << std::endl;

			return FileWriteResult::FileWriteOK;

		}


	};
}


#endif
