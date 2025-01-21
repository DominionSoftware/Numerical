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


#include <iomanip>
#include <iostream>
#include <numbers>
#include <ostream>
#include <gtest/gtest.h>
#include <sciplot/Canvas.hpp>

#include "AbscissaeFromData.h"
#include "GaussKronrodRules.h"
#include "Integral.h"
#include "linspace.h"
#include "valarray2vector.h"


TEST(NumericalTestSuite, TestGaussKonrod1)
{
    auto f = [](double x)->double
        {
            return sin(x);
        };

	numerical::GuassKronrodParameters<double> params(0, std::numbers::pi,100);

	std::tuple<numerical::GuassKronrodErrorCode, const char *, double> result =   IntegrateGaussKronrod( params, f, numerical::dk15);

    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::NoError);

    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::get<2>(result) << std::endl;


    numerical::GuassKronrodParameters<double> params2(0, std::numbers::pi, 1000);

    std::tuple<numerical::GuassKronrodErrorCode, const char*, double> result2 = IntegrateGaussKronrod(params2, f, numerical::dk21,  0.0,  1.0e-8);
    std::cout << std::get<1>(result2) << std::endl;
    std::cout << std::get<2>(result2) << std::endl;

    numerical::GuassKronrodErrorCode ec = static_cast<numerical::GuassKronrodErrorCode>(std::get<0>(result2));

    bool ok = (std::get<0>(result2) == numerical::GuassKronrodErrorCode::NoError);

    ASSERT_TRUE(ok);

}

TEST(NumericalTestSuite, TestGaussKonrod2)
{
    auto f = [](double x)->double
        {
            return std::exp(x) * std::log(x);
        };

    numerical::GuassKronrodParameters<double> params(0,1,2000);

    const std::tuple<numerical::GuassKronrodErrorCode,const char *, double> result = IntegrateGaussKronrod( params, f, numerical::dk15);

    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::NoError);

    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::get<2>(result) << std::endl;

}

TEST(NumericalTestSuite, TestGaussKonrod3)
{

    std::vector<std::function<double(double)>> functions = {
        {[](double x)->double
                {
                    return sin(20000 * std::numbers::pi * x);
                }
            }

    };


    numerical::GuassKronrodParameters<double> params(0, 1, 1000);

    std::tuple<numerical::GuassKronrodErrorCode, const char *, double> result = IntegrateGaussKronrod(params, functions[0], numerical::dk21);
    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::get<2>(result) << std::endl;
    
    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::RoundoffError);

}


TEST(NumericalTestSuite, TestGaussKonrod4)
{
    auto f = [](double x)->double
        {
            return (x * x * x) + (x * x) + x + 1;
        };

    numerical::GuassKronrodParameters<double> params(-1, 1, 2000);

    const std::tuple<numerical::GuassKronrodErrorCode, const char*, double> result = IntegrateGaussKronrod( params, f, numerical::dk15);


    std::cout << "I(x^3 + x^2 + x + 1) = " <<  std::get<2>(result) << std::endl;
    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::NoError);
    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::get<2>(result) << std::endl;


}
//(lambda x : x * *2, 0, 1, 1 / 3), # Integral of x ^ 2 from 0 to 1
//(lambda x : np.sin(x), 0, np.pi, 2), # Integral of sin(x) from 0 to pi
//(lambda x : np.exp(x), 0, 1, np.e - 1), # Integral of e^ x from 0 to 1
//(lambda x : 1 / x, 1, 2, np.log(2))  # Integral of 1 / x from 1 to 2

TEST(NumericalTestSuite, TestGaussKonrod5)
{
	  std::function<double(double)> function = 
      {[](double x)->double
              {
                  return exp(x);
              }
      };
    numerical::GuassKronrodParameters<double> params(0, 1, 1000);

    std::tuple<numerical::GuassKronrodErrorCode, const char*, double> result = IntegrateGaussKronrod(params, function, numerical::dk21);
    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::get<2>(result) << std::endl;

    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::NoError);
    ASSERT_NEAR(1.71828182845905, std::get<2>(result), 1.0e-8);

}

TEST(NumericalTestSuite, TestGaussKonrod6)
{
    std::function<double(double)> function =
    { [](double x)->double
            {
                return 1/x;
            }
    };
    numerical::GuassKronrodParameters<double> params(1, 2, 1000);

    std::tuple<numerical::GuassKronrodErrorCode, const char*, double> result = IntegrateGaussKronrod(params, function, numerical::dk21);
    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::get<2>(result) << std::endl;

    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::NoError);
   
	ASSERT_NEAR(0.693147180559945, std::get<2>(result), 1.0e-8);
}

TEST(NumericalTestSuite, TestGaussKonrod7)
{
    std::function<double(double)> function =
    { [](double x)->double
            {
                return exp(x);
            }
    };
    numerical::GuassKronrodParameters<double> params(0, 1, 1000);

    std::tuple<numerical::GuassKronrodErrorCode, const char*, double> result = IntegrateGaussKronrod(params, function, numerical::dk31);
    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::setprecision(14) << 1.71828182845905<< " "<<  std::get<2>(result) << std::endl;

    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::NoError);
    ASSERT_NEAR(1.71828182845905, std::get<2>(result), 1.0e-14);
}

TEST(NumericalTestSuite, TestGaussKonrod8)
{
    std::function<double(double)> function =
    { [](double x)->double
            {
                return sin(x);
            }
    };
    numerical::GuassKronrodParameters<double> params(0, std::numbers::pi, 1000);

    std::tuple<numerical::GuassKronrodErrorCode, const char*, double> result = IntegrateGaussKronrod(params, function, numerical::dk41);
    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::setprecision(14) << 2.0 << " " << std::get<2>(result) << std::endl;

    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::NoError);
    ASSERT_NEAR(2.0, std::get<2>(result), 1.0e-14);
}

TEST(NumericalTestSuite, TestGaussKonrod9)
{
    std::function<double(double)> function =
    { [](double x)->double
            {
                return sin(x);
            }
    };
    numerical::GuassKronrodParameters<double> params(0, std::numbers::pi, 1000);

    std::tuple<numerical::GuassKronrodErrorCode, const char*, double> result = IntegrateGaussKronrod(params, function, numerical::dk51);
    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::setprecision(14) << 2.0 << " " << std::get<2>(result) << std::endl;

    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::NoError);
    EXPECT_DOUBLE_EQ(2.0, std::get<2>(result));
}

TEST(NumericalTestSuite, TestGaussKonrod10)
{
    std::function<double(double)> function =
    { [](double x)->double
            {
                return sin(x);
            }
    };
    numerical::GuassKronrodParameters<double> params(0, std::numbers::pi, 1000);

    std::tuple<numerical::GuassKronrodErrorCode, const char*, double> result = IntegrateGaussKronrod(params, function, numerical::dk61);
    std::cout << std::get<1>(result) << std::endl;
    std::cout << std::setprecision(14) << 2.0 << " " << std::get<2>(result) << std::endl;

    ASSERT_EQ(std::get<0>(result), numerical::GuassKronrodErrorCode::NoError);
    EXPECT_DOUBLE_EQ(2.0, std::get<2>(result));
}

TEST(NumericalTestSuite, TestSimpsons)
{

    auto f = [](double x)->double
        {
            return (x * x * x) + (x * x) + x + 1;
        };


    auto x = numerical::Linspace<double>(-1, 1, 1000000);

    std::vector<double> diffs;

    std::adjacent_difference(x.begin(), x.end(), std::back_inserter(diffs));

    numerical::AbscissaeFromData<double> abscissae(x);

    auto s = numerical::IntegrateSimpsons<double, double>(&abscissae, f);

    EXPECT_NEAR(s, 2.66667, 1e-4);
    auto s2 = numerical::IntegrateSimpsons<double, double>(-1.0, 1.0, 1000000, f);
    EXPECT_NEAR(s2, 2.66667, 1e-4);
    std::cout << s2 << std::endl;


    sciplot::Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");





    plot.yrange(-1, 8);
    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2)
        .displayExpandHeightBy(2);


    std::vector<double> d = numerical::Linspace<double>(-1.0, 1.0, 100);
    std::vector<double> curve;
    std::transform(std::begin(d), std::end(d), std::back_inserter(curve), f);
    plot.drawCurve(numerical::vector2valarray(d), numerical::vector2valarray(curve)).label("function");


    // Create figure to hold plot
    sciplot::Figure fig = { {plot} };
    // Create canvas to hold figure
    sciplot::Canvas canvas = { {fig} };
    canvas.size(1200, 1800);
    // Show the plot in a pop-up window
    canvas.show();


}

TEST(NumericalTestSuite, TestTrapezoid)
{
    auto x = numerical::Linspace<double>(-1, 1, 1000);

    numerical::AbscissaeFromData<double> abscissae(x);

    auto f = [](double x)->double
        {
            return (x * x * x) + (x * x) + x + 1;
        };

    auto s = numerical::IntegrateTrapezoid<double, double>(&abscissae, f);
    std::cout << s << std::endl;  // NOLINT(performance-avoid-endl)
    sciplot::Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");

    plot.xrange(x[0], x[x.size() - 1]);



    plot.yrange(-1, 8);
    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2)
        .displayExpandHeightBy(2);


    std::vector<double> xs = abscissae.data();
    std::vector<double> r;

    std::transform(std::begin(xs), std::end(xs), std::back_inserter(r), [&](double x)
        {
            return f(x);
        });

    plot.drawCurveWithPoints(numerical::vector2valarray(x), numerical::vector2valarray(r)).label("trapezoid");

    // Create figure to hold plot
    sciplot::Figure fig = { {plot} };
    // Create canvas to hold figure
    sciplot::Canvas canvas = { {fig} };
    canvas.size(1200, 1800);
    // Show the plot in a pop-up window
    canvas.show();


}

TEST(NumericalTestSuite, TestGaussianLegendre)
{


    auto func = [](double x)
        {
            return (x * x * x) + (x * x) + x + 1;

        };


    double result = numerical::IntegrateGaussLegendre(-1, 1, 5, func);


    std::cout << result << std::endl;

}


