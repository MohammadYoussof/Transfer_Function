# Transfer_Function
MTF determines how much contrast in the original object is preserved by the imager. In other words, it describes how accurately the object's spatial frequency content is translated to the image.

#######
This is a readme file explaining the C++ functions that have been converted from the original MATLAB codes.

Translating MATLAB to C++ is a non-trivial task that requires a good deal of time. This is due to the different levels of complexity and the different structures used in MATLAB and C++.

However, if I can provide the definitions and implementations of all the functions used in the code (such as light_source_setting, paraxial_focal_length, data_reshape, view_lens, spot_diagram, transmission_plane, line_spread_function, MTF, point_spread_function, etc.), I can then target a more acceptable translation into C++.

The current MATLAB and Python codes perform a series of operations related to laser ray tracing and optical computations, including:

= Defining various optical parameters.
= Setting up a light source.
= Calculating paraxial focal length.
= Performing ray tracing.
= Displaying various visualizations such as lens view, spot diagram, optical path difference at the transmission plane, etc.
= Calculating and displaying line spread function (LSF), modulation transfer function (MTF), and point spread function (PSF).

There is a significant usage of MATLAB-specific functionality such as matrix operations, built-in functions, and cell arrays, which do not have a very direct equivalents in C++. Therefore, the conversion will require a good flexibility and understanding of both language and may still need significant changes to make the code more suitable for C++ work flow.

I provide here an outline of how to approach this conversion. A full conversion of this code would be quite extensive and require a considerable amount of time. I guess looking on the existing and functioning C++ code that is working with LSF will resolve the purpose of grafting a new MTF routine. 

Here is the general process for converting this MATLAB code into C++:

1-	Data types: 
MATLAB uses dynamic typing, while C++ is a statically typed language. I need to specify data types for all your variables. For example, the MATLAB variable lambda (λ) would be defined in C++ as double .
2-	Arrays and matrices: 
MATLAB has built-in support for matrix operations, but C++ does not. I use libraries like Eigen, Armadillo, or raw arrays with manual implementation of matrix operations.
3-	Functions: 
MATLAB functions can return multiple values, but in C++, functions can only return one value. I use pass-by-reference parameters or std::tuple to return multiple values in C++.
4-	Libraries: 
Some MATLAB functions might have no direct equivalent in C++, so I may need to use external libraries or write own functions. For example, num2str can be replaced with std::to_string.
5-	Loops and conditionals: 
These can be translated relatively directly between MATLAB and C++, but C++ still uses zero-based indexing, not one-based indexing like MATLAB.
6-	Cell arrays: 
MATLAB cell arrays do not have a direct equivalent in C++. I will need to decide on a suitable alternative, such as std::vector or std::array.
7-	Memory management: 
MATLAB manages memory automatically, but in C++, I will manually manage memory, especially for dynamic data structures.
8-	Input/Output: 
Replace MATLAB disp function with C++ std::cout.

To translate this code into C++, each of these functions would need to be rewritten in C++.
This requires not just the knowledge of what these functions do, but also how they are implemented.  

I can provide a structure on how to convert the loops and operations in the MATLAB code into C++.  

Given the main MATLAB code, a first rough skeleton will be as bellow:

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

using namespace std;

int main() {
    // Parameter Setting
    double lambda = 546.1; // Unit : nm
    double sk16_schott = 1.62286;
    int surface_num = 2;
    vector<double> distance = {0.1, 0.02, 0.15}; // Unit : mm
    vector<double> material = {1, sk16_schott, 1}; // Unit : mm
    vector<double> y_radius = {numeric_limits<double>::infinity(), -0.1}; // Unit : mm
    double aperture = 0.05; // Unit : mm (Pupil Diameter)

    // Additional parameters
    double ang_x = 0;
    double ang_y = 0;
    int cross_diameter_num = 201;

    // Switches
    int Use_Paraxial_Solve = 1;
    int View_Lens = 1;
    int viewplane = 1;
    int display_line = 11;
    int Spot_Diagram = 1;
    int Transmission_Plane = 1;
    int Line_Spread_Function = 1;
    int MTF = 1;
    vector<int> Point_Spread_Function = {1, 1};

    // Source Setting
    lambda = lambda * 1e-6;   // nm -> mm

    // Here, we should call light_source_setting function with appropriate arguments
    // For example:
    // tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> result = light_source_setting(aperture, distance, cross_diameter_num, ang_x, ang_y);
    // And then, we need to unpack the returned tuple into variables s_x, s_y, s_z, L, M, N

    // And so on...

    return 0;
}

==========================================

The attached additional Python and MATLAB codes contain functions for calculating the MTF of simplified figures.

The code starts by loading an image from a specified file path and returning the image as grayscale data.

The LoadImageAsArray function converts it into a numpy array, with pixel values ranging between 0.0 and 1.0.
The ImageToArray function converts a PIL image to a numpy array. The ArrayToImage function converts a numpy array back to a PIL image. The code adjusts the orientation of an image for analysis purposes.

The GetEdgeSpreadFunctionCrop function calculates the Edge Spread Function (ESF) for a slanted edge target. It detects the edge in the image, calculates the ESF, and crops it around the center of the transition. The SimplifyEdgeSpreadFunction function removes duplicate distance occurrences in the ESF to simplify the data.

The GetLineSpreadFunction function calculates the Line Spread Function (LSF) from the ESF by taking the derivative.
The GetMTF function calculates MTF from the LSF by performing the Fourier transform and normalizing the values.

The CalculateMtf function combines these steps to calculate the MTF of an image array, including the necessary preprocessing steps.

Overall, these functions provide a set of tools for analyzing simple figures and calculating their MTF as a measure of image quality. The next step involves more complicated tests.
I still need to develop a control loop or automation routine to iterate through the focus settings, acquire images, calculate the MTF, and adjust the focus until the best focus resolution is achieved. This could be done in a feedback loop, continuously refining the focus based on MTF measurements.
