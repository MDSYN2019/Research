// Standard libraries

/**
- *  This is where we put in the text 
 *  e^{} = e^{}
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <set>

// OpenCV libaries

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "MD_alg.hpp"
#include "MD_boostlib.hpp"
#include "MD_friend_function.hpp"
#include "MD_gsl.hpp"
#include "MD_gslmat.hpp"
#include "MD_object.hpp"

// --- GSL libraries --- //

// --- Boost libraries --- //

#include <boost/array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/distributions.hpp>
#include <boost/current_function.hpp>
#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <boost/detail/lightweight_test.hpp>
#include <boost/lexical_cast.hpp>

// --- Algorithms libraries --- //

using namespace boost::numeric::odeint;
using namespace boost;
using namespace cv;
using namespace std;

int main(int argc, char *argv[]) {
    if ( argc != 2 )
    {
        printf("usage: DisplayImage.out <Image_Path>\n");
        return -1;
    }

    Mat image;

    image = imread( argv[1], 1 );

    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }
    namedWindow("Display Image", WINDOW_AUTOSIZE );
    imshow("Display Image", image);

    waitKey(0);

    
return 0;    
}