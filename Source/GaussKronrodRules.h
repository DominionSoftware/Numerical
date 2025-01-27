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

#ifndef GAUSS_KRONROD_RULES_
#define GAUSS_KRONROD_RULES_
#include <functional>
#include <vector>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <span>
#include <numeric>
#include <cmath>
#include "GaussKronrodConstants.h"

namespace numerical
{
	
	enum class GuassKronrodErrorCode {
		NoError,
		MaxSubIntervalsReached,
		RoundoffError,
		ConvergenceNotAchieved,
		InvalidInput,
		BadIntegrandBehavior,
		UnableToIntegrate
	};

	/**
	 * @brief class GuassKronrodParameters 
	 * Based on from Fortran code from netlib.org.
	 *
	 * This class and methods is based on the algorithm from the book:
	 * R. Piessens, E. de Doncker-Kapenga, C. W. ?berhuber, and D. K. Kahaner,
	 * "QUADPACK: A Subroutine Package for Automatic Integration",
	 * Springer-Verlag, New York, NY, 1983.
	 *
	 * @note The original Fortran source can be found at netlib.org.
	 *
	 * @see QUADPACK: A Subroutine Package for Automatic Integration.
	 *
	 * The algorithm uses a vector - index for ordering the errors in descending
	 * order. This requires pre-allocating vectors in constructor.
	 */

	template<std::floating_point T>
	class GuassKronrodParameters
	{
	public:
		explicit GuassKronrodParameters(T a,T b,size_t limit) : initialLowerBounds_(a),initialUpperBounds_(b), LIMIT_(limit), SIZE_(0), MAXNR_(0), MAXERR_(0),
		                                                        MAXLEVEL_(0), ALIST_(limit, 0.0), BLIST_(limit, 0.0),
		                                                        RLIST_(limit, 0.0), ELIST_(limit, 0.0),
		                                                        IORD_(limit, 0), LEVEL_(limit, 0)
		{
			ALIST_[0] = a;
			BLIST_[0] = b;
		}

		GuassKronrodParameters() = delete;

		T getInitialLowerBounds()
		{
			return initialLowerBounds_;
		}

		T getInitialUpperBounds()
		{
			return initialUpperBounds_;
		}

		void getValuesAtMaxError(T& a,T& b,T& r,T& e)
		{
			a = ALIST_[MAXERR_];
			b = BLIST_[MAXERR_];
			r = RLIST_[MAXERR_];
			e = ELIST_[MAXERR_];
		}

		void setFirstValues(T result,T abserr)
		{
			SIZE_ = 1;
			RLIST_[0] = result;
			ELIST_[0] = abserr;
			IORD_[0] = 0;
		}

		size_t limit() const
		{
			return LIMIT_;
		}

		T sumResults()
		{
			auto vi = std::span(RLIST_.begin(), RLIST_.begin() + SIZE_);
			return std::accumulate(vi.begin(),vi.end(),0.0);
		}

		void appendToLists(double a1, double b1, double area1, double error1,
			double a2, double b2, double area2, double error2)
		{

			const size_t indexMax = MAXERR_; // set value
			const size_t indexAdd = SIZE_;	 // add value.

			const size_t newLevel = LEVEL_[indexMax] + 1;


			if (error2 > error1)
			{
				ALIST_[indexMax] = a2;        
				RLIST_[indexMax] = area2;
				ELIST_[indexMax] = error2;
				LEVEL_[indexMax] = newLevel;

				ALIST_[indexAdd] = a1;
				BLIST_[indexAdd] = b1;
				RLIST_[indexAdd] = area1;
				ELIST_[indexAdd] = error1;
				LEVEL_[indexAdd] = newLevel;
			}
			else
			{
				BLIST_[indexMax] = b1;       
				RLIST_[indexMax] = area1;
				ELIST_[indexMax] = error1;
				LEVEL_[indexMax] = newLevel;

				ALIST_[indexAdd] = a2;
				BLIST_[indexAdd] = b2;
				RLIST_[indexAdd] = area2;
				ELIST_[indexAdd] = error2;
				LEVEL_[indexAdd] = newLevel;
			}

			SIZE_++;

			if (newLevel > MAXLEVEL_)
			{
				MAXLEVEL_ = newLevel;
			}

			// resort not that they are added....
			qpsrt();
		}

	private:

		/**
		 * @brief Ported from Fortran code qpsort.f from netlib.org.
		 * and from the book:
		 * 
		 * R. Piessens, E. de Doncker-Kapenga, C. W. ?berhuber, and D. K. Kahaner,
		 * "QUADPACK: A Subroutine Package for Automatic Integration",
		 * Springer-Verlag, New York, NY, 1983.
		 *
		 * see also:
		 * D. Kahaner, C. B. Moler, and S. Nash, 
		 * Numerical Methods and Software. 
		 * Englewood Cliffs, NJ: Prentice Hall, 1989.

		 * @note Fortran source can be found at netlib.org.
		 *
		 * @see QUADPACK: A Subroutine Package for Automatic Integration.
		 *
		 * note: all caps used to better map to fortan code in book.
		 */
		 
		void qpsrt()
		{

			size_t NRMAX = MAXNR_;
			size_t MAXERR = IORD_[NRMAX];
			const size_t LAST = SIZE_ - 1;

			/* CHECK WHETHER THE LIST CONTAINS MORE THAN TWO ERROR ESTIMATES */

			if (LAST < 2)
			{
				IORD_[0] = 0;
				IORD_[1] = 1;
				MAXERR_ = MAXERR;
				return;
			}

			double ERRMAX = ELIST_[MAXERR];

			/* THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A DIFFICULT
			   INTEGRAND, SUBDIVISION INCREASED THE ERROR ESTIMATE. IN THE NORMAL
			   CASE THE INSERT PROCEDURE SHOULD START AFTER THE NRMAX-TH LARGEST
			   ERROR ESTIMATE.
			*/

			while (NRMAX > 0 && ERRMAX > ELIST_[IORD_[NRMAX - 1]])
			{
				IORD_[NRMAX] = IORD_[NRMAX - 1];
				NRMAX--;
			}

			/* COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED IN
			   DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
			   SUBDIVISIONS STILL ALLOWED.
			*/
			const size_t LIMIT = LIMIT_;
			size_t top;

			if (LAST < (LIMIT / 2 + 2))
			{
				top = LAST;
			}
			else
			{
				top = LIMIT - LAST + 1;
			}

			/* INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN, STARTING
			   COMPARISON FROM THE ELEMENT ELIST(ORDER(I_NRMAX+1)).
			*/

			size_t I = NRMAX + 1;

			while (I < top && ERRMAX < ELIST_[IORD_[I]])
			{
				IORD_[I - 1] = IORD_[I];
				I++;
			}

			IORD_[I - 1] = MAXERR;

			/* INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP */

			double ERRMIN = ELIST_[LAST];

			size_t k = top - 1;

			while (k > I - 2 && ERRMIN >= ELIST_[IORD_[k]])
			{
				IORD_[k + 1] = IORD_[k];
				k--;
			}

			IORD_[k + 1] = LAST;


			MAXERR = IORD_[NRMAX];

			MAXERR_ = MAXERR;
			MAXNR_ = NRMAX;
		}


		T initialLowerBounds_;
		T initialUpperBounds_;

		size_t LIMIT_;
 		size_t SIZE_;
		size_t MAXNR_;
		size_t MAXERR_;
		size_t MAXLEVEL_;

		std::vector<T> ALIST_;
		std::vector<T> BLIST_;
		std::vector<T> RLIST_;
		std::vector<T> ELIST_;
		std::vector<size_t> IORD_;
		std::vector<size_t> LEVEL_;
	};

	/**
	 * @brief This function implements the dkNN algorithm, ported from Fortran code at netlib.org.
	 *
	 * It calculates numerical integration using Kronrod quadrature, within the limits of integration.
	 *
	 * @param a The lower limit of integration.
	 * @param b The upper limit of integration.
	 *
	 * The result is calculated as:
	 * \f[
	 * result = \text{Kronrod} \times \frac{b - a}{2}
	 * \f]
	 *
	 * \f$ resabs \f$ is the approximation of the integral as:
	 * \f[
	 * resabs = \sum_{a}^{b} |f(x)|
	 * \f]
	 *
	 * \f$ resasc \f$ is the approximation of the absolute error as:
	 * \f[
	 * resasc = \sum_{a}^{b} |f(x) - \frac{result}{b - a}|
	 * \f]
	 *
	 *
	 *
	 * result = konrond x (b - a) / 2
	 *
	 *				    b
	 * resabs  approx = sigma abs(f(x))
	 *				    a
	 *
	 *
	 *				    b
	 * resasc  approx = sigma abs(f(x) - result/(b-a)
	 *                  a
	 */
	template<typename T>
		requires std::floating_point<T>
	void dk(const std::function<T(T)>& func, T a, T b, const std::vector<T>& xgk, const std::vector<T>& wg, const std::vector<T>& wgk, T& result, T& abserr, T& resabs, T& resasc)
	{
		auto epsilonMachine = std::numeric_limits<T>::epsilon();
		auto underFlow = std::numeric_limits<T>::min();

		T resultAbs{ 0 };
		T resAsc{ 0 };
		T center = 0.5 * (a + b);
		T halfLength = 0.5 * (b - a);
		T absHalfLength = fabs(halfLength);

		T midPoint = func(center);
		T sumGauss = midPoint * wg[wg.size() - 1];
		T sumGaussKronrod = midPoint * wgk[wgk.size() - 1];

		resultAbs = std::fabs(sumGaussKronrod);

		std::vector<T> functionValuesLeft(wgk.size(), T());
		std::vector<T> functionValuesRight(wgk.size(), T());

		// Compute Kronrod approximation and estimate absolute error
		const int loopEnd1 = static_cast<int>(wgk.size() - 1) / 2;

		// fortran J=1,3  JTW = J*2
		// i.e., 2,4,6
		// in C++ 1,3,5 so
		// j = 0 => 1
		// j = 1 => 3
		// j = 2 => 5
		//
		for (int j = 0; j < loopEnd1; j++) {
			int JTW = j * 2 + 1;
			const T abscissaOffset = halfLength * xgk[JTW];
			const T functionValueLeft = func(center - abscissaOffset);
			const T functionValueRight = func(center + abscissaOffset);
			functionValuesLeft[JTW] = functionValueLeft;
			functionValuesRight[JTW] = functionValueRight;
			const T sumLeftRight = functionValueLeft + functionValueRight;
			sumGauss += wg[j] * sumLeftRight;
			sumGaussKronrod += wgk[JTW] * sumLeftRight;
			resultAbs += wgk[JTW] * (fabs(functionValueLeft) + fabs(functionValueRight));
		}
		// in fortran, j=1,4  jtwm1 = j*2-1
		// i.e. 1,3,5
		// in C++ 0,2,4 so
		// j = 0 => 0
		// j = 1 => 2
		// j = 2 => 4
		for (int j = 0; j < wgk.size() / 2; j++) {
			int JTWM1 = j * 2;
			T abscissaOffset = halfLength * xgk[JTWM1];
			T functionValueLeft = func(center - abscissaOffset);
			T functionValueRight = func(center + abscissaOffset);
			functionValuesLeft[JTWM1] = functionValueLeft;
			functionValuesRight[JTWM1] = functionValueRight;
			T sumLeftRight = functionValueLeft + functionValueRight;
			sumGaussKronrod += wgk[JTWM1] * sumLeftRight;
			resultAbs += wgk[JTWM1] * (fabs(functionValueLeft) + fabs(functionValueRight));
		}

		const T mean = sumGaussKronrod * 0.5;
		resAsc = wgk[wgk.size() - 1] * fabs(midPoint - mean);

		for (int j = 0; j < wgk.size() - 1; j++) {
			resAsc += wgk[j] * (fabs(functionValuesLeft[j] - mean) + fabs(functionValuesRight[j] - mean));
		}

		// set ouputs
		T err = (sumGaussKronrod - sumGauss) * halfLength;
		result = sumGaussKronrod * halfLength;
		resabs = resultAbs *= absHalfLength;
		resasc = resAsc *= absHalfLength;

		abserr = std::fabs((sumGaussKronrod - sumGauss) * halfLength);
		/*
		 * Experimentation has suggested that
		 *		Integral Estimate = Konrod 2n + 1
		 *		Error Estimate = (200 * abs(Gn - K2n + 1)) ^ 1.5
		 *
		 *		D. Kahaner, C. B. Moler, and S. Nash,
		 *		Numerical Methods and Software.
		 *		Englewood Cliffs, NJ: Prentice Hall, 1989, p.154.
		 *
		 */

		if (resAsc != 0.0 && err != 0.0)
		{
			T scale =  std::pow(200.0 * err / resAsc, 1.5);
			if (scale < 1)
			{
				abserr = resAsc * scale;
			}
			else
			{
				abserr = resAsc;
			}
		}

		if (resultAbs > underFlow / (50.0 * epsilonMachine))
		{
			T minErr =  50.0 * epsilonMachine * resultAbs;
			if (minErr > err)
			{
				abserr = minErr;
			}
		}
	}
	


	template<typename T>
	requires std::floating_point<T>
	void dk15(const std::function<T(T)>& func, T a, T b,T& result,T& abserr,T& resabs,T& resasc)
	{
		dk<T>(func, a, b, numerical::dqk15::xgk, numerical::dqk15::wg, numerical::dqk15::wgk, result, abserr, resabs, resasc);
	}

	template<typename T>
		requires std::floating_point<T>
	void dk21(const std::function<T(T)>& func, T a, T b, T& result, T& abserr, T& resabs, T& resasc)
	{
		dk<T>(func, a, b, numerical::dqk21::xgk, numerical::dqk21::wg, numerical::dqk21::wgk, result, abserr, resabs, resasc);
	}

	template<typename T>
		requires std::floating_point<T>
	void dk31(const std::function<T(T)>& func, T a, T b, T& result, T& abserr, T& resabs, T& resasc)
	{
		dk<T>(func, a, b, numerical::dqk31::xgk, numerical::dqk31::wg, numerical::dqk31::wgk, result, abserr, resabs, resasc);
	}
	template<typename T>
		requires std::floating_point<T>
	void dk41(const std::function<T(T)>& func, T a, T b, T& result, T& abserr, T& resabs, T& resasc)
	{
		dk<T>(func, a, b, numerical::dqk41::xgk, numerical::dqk41::wg, numerical::dqk41::wgk, result, abserr, resabs, resasc);
	}

	template<typename T>
		requires std::floating_point<T>
	void dk51(const std::function<T(T)>& func, T a, T b, T& result, T& abserr, T& resabs, T& resasc)
	{
		dk<T>(func, a, b, numerical::dqk51::xgk, numerical::dqk51::wg, numerical::dqk51::wgk, result, abserr, resabs, resasc);
	}

	template<typename T>
		requires std::floating_point<T>
	void dk61(const std::function<T(T)>& func, T a, T b, T& result, T& abserr, T& resabs, T& resasc)
	{
		dk<T>(func, a, b, numerical::dqk61::xgk, numerical::dqk61::wg, numerical::dqk61::wgk, result, abserr, resabs, resasc);
	}
}

#endif