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
namespace numerical
{
    namespace ClenshawCurtis
    {
        enum class ErrorCode {
            NoError,
            MaxSubIntervalsReached,
            RoundoffError,
            ConvergenceNotAchieved,
            InvalidInput,
            BadIntegrandBehavior,
            UnableToIntegrate
        };

        template<std::floating_point T>
        struct Parameters
        {
        };

        template<std::floating_point T>
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

        template<std::floating_point T>
        void GenerateChebyshev(const ClenshawCurtis::Function<T>& function, T lowerBound, T upperBound,
            std::array<T, 13>& chebyshevCoefficients12, std::array<T, 25>& chebyshevCoefficients24) {
            constexpr std::array<T, 11> chebyshevNodes = {
                0.9914448613738104, 0.9659258262890683, 0.9238795325112868,
                0.8660254037844386, 0.7933533402912352, 0.7071067811865475,
                0.6087614290087206, 0.5000000000000000, 0.3826834323650898,
                0.2588190451025208, 0.1305261922200516
            };
            std::array<T, 25> functionValues;
            std::array<T, 12> differences;
            const T intervalCenter = 0.5 * (upperBound + lowerBound);
            const T intervalHalfLength = 0.5 * (upperBound - lowerBound);
            functionValues[0] = 0.5 * function(upperBound);
            functionValues[12] = function(intervalCenter);
            functionValues[24] = 0.5 * function(lowerBound);
            for (size_t i = 1; i < 12; ++i)
            {
                const size_t oppositeIndex = 24 - i;
                const T nodeOffset = intervalHalfLength * chebyshevNodes[i - 1];
                functionValues[i] = function(intervalCenter + nodeOffset);
                functionValues[oppositeIndex] = function(intervalCenter - nodeOffset);
            }
            for (size_t i = 0; i < 12; ++i)
            {
                const size_t oppositeIndex = 24 - i;
                differences[i] = functionValues[i] - functionValues[oppositeIndex];
                functionValues[i] = functionValues[i] + functionValues[oppositeIndex];
            }
            // Compute cheb12 and cheb24 coefficients
            auto computeChebyshevCoefficients = [&]<size_t N>(std::array<T, N>&coefficients)
            {
                const T lambda1 = differences[0] - differences[8];
                const T lambda2 = chebyshevNodes[5] * (differences[2] - differences[6] - differences[10]);
                coefficients[3] = lambda1 + lambda2;
                coefficients[9] = lambda1 - lambda2;
            };
            computeChebyshevCoefficients(chebyshevCoefficients12);
            computeChebyshevCoefficients(chebyshevCoefficients24);
            // Scale coefficients
            for (size_t i = 1; i < 12; ++i)
            {
                chebyshevCoefficients12[i] *= 1.0 / 6.0;
            }
            chebyshevCoefficients12[0] *= 1.0 / 12.0;
            chebyshevCoefficients12[12] *= 1.0 / 12.0;
            for (size_t i = 1; i < 24; ++i)
            {
                chebyshevCoefficients24[i] *= 1.0 / 12.0;
            }
            chebyshevCoefficients24[0] *= 1.0 / 24.0;
            chebyshevCoefficients24[24] *= 1.0 / 24.0;
        }
    }
}
#endif // !CHEBYSHEV12_24_