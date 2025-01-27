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
#ifdef USE_SCIPLOT
#include <sciplot/Canvas.hpp>
#include <sciplot/Constants.hpp>
#include <sciplot/Plot2D.hpp>
#endif
#include "AbscissaeFromData.h"
#include "FiniteDifference.h"
#include "Linspace.h"
#include "valarray2vector.h"
#include "WriteVector.h"
#include "gtest/gtest.h"
#include <numbers>
using namespace numerical;

#define VISUAL_TEST


TEST(NumericalTestSuite, TestFiniteDifference1)
{
    // Function to test
    std::function<double(double)> func = [](double x)
        {
            return std::cos(x);
        };

    auto x = numerical::Linspace<double>(-std::numbers::pi, std::numbers::pi, 100);
    std::vector<double> csx(x.size());
    std::transform(std::begin(x), std::end(x), std::begin(csx), func);



    numerical::AbscissaeFromData<double> abscissaeFromData(x);

    // Calculate derivative using finite difference
    FiniteDifference<double, double> fd;
    auto fPrime = fd.gradient(func, &abscissaeFromData);

    ASSERT_EQ(fPrime.size(), x.size());  // Ensure the gradient is the same size as x

    // Optional: Check some known values (e.g., derivative of cos(x) is -sin(x))
    EXPECT_NEAR(fPrime[0], -std::sin(x[0]), 1e-4);
    EXPECT_NEAR(fPrime[x.size() / 2], -std::sin(x[x.size() / 2]), 1e-4);

#ifdef USE_SCIPLOT
    sciplot::Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");
    plot.xrange(x[0], x[x.size() - 1]);
    plot.yrange(-2, 2);
    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2)
        .displayExpandHeightBy(2);

    plot.drawCurve(numerical::vector2valarray(x), numerical::vector2valarray(csx)).label("cos x");
    plot.drawCurve(numerical::vector2valarray(x), numerical::vector2valarray(fPrime)).label("fPrime");

    auto err = FileIO<double>::WriteVector(fPrime, R"(D:\Projects\numerical\cosx_deriv.csv)");
    EXPECT_EQ(err, 0);

    // Create figure to hold plot
    sciplot::Figure fig = { {plot} };
    // Create canvas to hold figure
    sciplot::Canvas canvas = { {fig} };
    canvas.size(1200, 2000);

    canvas.show();
#endif
}


TEST(NumericalTestSuite, TestFiniteDifference2)
{
    auto abscissae = numerical::Linspace<double>(-4.0, 4, 50);


    std::vector<double> xSquared;

    std::for_each(std::begin(abscissae), std::end(abscissae), [&](double x)
        {
            xSquared.push_back(x * x);
        });


    numerical::FiniteDifference<double, double> fd;
    std::function<double(double)> func = [](double x)
        {
            return x * x;
        };

    numerical::AbscissaeFromData<double> abscissaeFromData(abscissae);

    auto fPrime = fd.gradient(func, &abscissaeFromData);


#ifdef USE_SCIPLOT

    sciplot::Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");

    plot.xrange(abscissae[0], abscissae[abscissae.size() - 1]);

    auto min_it = std::min_element(xSquared.begin(), xSquared.end());
    auto max_it = std::max_element(xSquared.begin(), xSquared.end());

    double min = *max_it * -1;

    plot.yrange(min, *max_it);
    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2)
        .displayExpandHeightBy(2);


    plot.drawCurve(numerical::vector2valarray(abscissae), numerical::vector2valarray(xSquared)).label("x squared");
    plot.drawCurve(numerical::vector2valarray(abscissae), numerical::vector2valarray(fPrime)).label("fPrime");

    // Create figure to hold plot
    sciplot::Figure fig = { {plot} };
    // Create canvas to hold figure
    sciplot::Canvas canvas = { {fig} };
    canvas.size(1200, 1800);
    // Show the plot in a pop-up window
    canvas.show();
#endif


}


TEST(NumericalTestSuite, TestFiniteDifference3)
{
    std::vector<double> yData = { -1,-0.839234918784,-0.6834639200335,-0.5326946364015,-0.3869231903375,-0.246146194086,-0.110370545606,0.0204072653060001,0.1461902590585,0.2669722723865,0.3827564481465,0.4935454394,0.599333817576,0.700124358184,0.7959193469385,0.8867140899625,0.9725109954185,1.053311981674,1.129113089546,1.19991635985,1.2657233436065,1.3265308163265,1.3823404514785,1.43315284089963,1.47896705805914,1.51978345805914,1.55560184089962,1.5864224514785,1.6122448163265,1.6330693436065,1.64889635985,1.659725089546,1.665555981674,1.6663889954185,1.6622240899625,1.6530613469385,1.638900358184,1.619741817576,1.5955854394,1.5664304481465,1.5322782723865,1.4931282590585,1.448979265306,1.399833454394,1.345689805914,1.2865468096625,1.2224073635985,1.1532700799665,1.079133081216,1
    };
    std::function<double(double)> func = [](double x)
        {
            constexpr double p1 = -0.0025;
            constexpr double p2 = 0.1683;
            constexpr double p3 = -1.1568;

            return (p1 * x * x) + (p2 * x) + p3;

        };


    auto x = numerical::Linspace<double>(0, 50, 100);
    std::vector<double> curve(x.size());

    std::transform(std::begin(x), std::end(x), std::begin(curve), func);


    numerical::FiniteDifference<double, double> fd;
    numerical::AbscissaeFromData<double> fromData(x);

    auto g = fd.gradient(func, &fromData);
    ASSERT_EQ(curve.size(), x.size());
    ASSERT_EQ(g.size(), x.size());
#ifdef USE_SCIPLOT
    sciplot::Plot2D plot;
    plot.xlabel("x");
    plot.ylabel("y");

    plot.xrange(x[0], x[x.size() - 1]);



    plot.yrange(-1, 2);
    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2)
        .displayExpandHeightBy(2);


    plot.drawCurve(numerical::vector2valarray(x), numerical::vector2valarray(curve)).label("x squared");
    plot.drawCurve(numerical::vector2valarray(x), numerical::vector2valarray(g)).label("fPrime");

    // Create figure to hold plot
    sciplot::Figure fig = { {plot} };
    // Create canvas to hold figure
    sciplot::Canvas canvas = { {fig} };
    canvas.size(1200, 1800);
    // Show the plot in a pop-up window
    canvas.show();

    #endif
}

