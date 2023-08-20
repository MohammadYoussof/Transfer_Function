 // I will give here a high-level structure of the C++ main code based on the translated functions and the structure of the original MATLAB main code.
 // Therefore, this may give a good starting point to further refine, and debug the C++ code in the local environment.

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <numeric>
#include "Eigen/Dense" // for matrix operations, you need to install and include the Eigen library

using namespace std;
using namespace Eigen;

// Here, you should define or replace the `Lens`, `Data`, `TransPlaneData`, `LSFData`, and `PSFData` types
// according to the actual C++ program structure.

// Here, you should include all the previously translated functions. They are just prototypes and
// might not work exactly like the MATLAB functions. You need to adjust and debug them according to your needs and the structure of your codes.

// light_source_setting function
vector<double> light_source_setting(double aperture, vector<double> distance, int cross_diameter_num, double ang_x, double ang_y) {
    // function body.
}

// paraxial_focal_length function
pair<double, double> paraxial_focal_length(int surface_num, vector<double> distance, vector<double> material, vector<double> sur_radius) {
    // function body..
}

// data_reshape function
Data data_reshape(vector<vector<double>> x_all, vector<vector<double>> y_all, vector<vector<double>> z_all, int cross_diameter_num) {
    // function body...
}

// trans_plane_position_and_optical_path function
TransPlaneData trans_plane_position_and_optical_path(int surface_num, vector<double> distance, vector<double> material, Data data, vector<vector<double>> L, vector<vector<double>> M, vector<vector<double>> N) {
    // function body....
}

// diffraction_limit function
double diffraction_limit(double lambda, double aperture, double BFL, double EFL, double focal_plane_position) {
    // function body.....
}

// line_spread_function function
LSFData line_spread_function(double lambda, double aperture, TransPlaneData trans_plane_data, double focal_plane_position) {
    // function body......
}

// point_spread_function function
PSFData point_spread_function(double lambda, double aperture, TransPlaneData trans_plane_data, double propagate_distance, vector<int> perspective_switch) {
    // function body.......
}

// display_tools class  [I am not sure]
class display_tools {
public:
    static void view_lens(Lens lens, Data data, int display_line, std::string viewplane) {
        // Visualization code here.
    }

    static void spot_diagram(Data data) {
        // Visualization code here..
    }

    static void transmission_plane(TransPlaneData trans_plane_data) {
        // Visualization code here...
    }

    static void line_spread_function(LSFData LSF_data, double diffra_limit) {
        // Visualization code here....
    }

    static void MTF(Lens lens, LSFData LSF_data, double diffra_limit, double ang_y) {
        // Visualization code here.....
    }

    static void point_spread_function(bool Switch, PSFData PSF_data, TransPlaneData trans_plane_data, double focal_plane_position) {
        // Visualization code here......
    }
};

// main function
int main() {
    // Here you can follow the structure of the MATLAB main script, calling the above defined functions at appropriate places.
    // Remember, there are several things that needed to be adjusted and debugged according to the needs, and also you need to install
    // and include necessary libraries for matrix operations and data visualization.

    return 0;
}
