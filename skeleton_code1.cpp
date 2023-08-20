#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <numeric>
#include <limits>

using namespace std;

struct TransPlaneData {
    vector<double> y, OP;
};

struct LSFData {
    vector<double> Monitor_y, Intensity, Intensity_normalize, Power, Power_normalize;
    double Strehl_ratio;
};

LSFData line_spread_function(double lambda, double aperture, TransPlaneData trans_plane_data, double focal_plane_position) {
    // Body of the line_spread_function translated earlier
}

pair<double, double> paraxial_focal_length(int surface_num, vector<double> distance, vector<double> material, vector<double> sur_radius) {
    // Body of the paraxial_focal_length function translated earlier
}

struct Data {
    vector<vector<vector<double>>> X_2, Y_2, Z_2;
};

Data data_reshape(vector<vector<double>> x_all, vector<vector<double>> y_all, vector<vector<double>> z_all, int cross_diameter_num) {
    // Body of the data_reshape function translated earlier
}

struct PSFData {
    vector<double> Monitor1_y, Monitor_z;
    vector<vector<double>> Intensity, Intensity_normalize, Power, Power_normalize;
    double Strehl_ratio;
};

PSFData point_spread_function(double lambda, double aperture, TransPlaneData trans_plane_data, double propagate_distance, vector<int> perspective_switch) {
    // Body of the point_spread_function function translated earlier
}

TransPlaneData trans_plane_position_and_optical_path(int surface_num, vector<double> distance, vector<double> material, Data data, vector<vector<double>> L, vector<vector<double>> M, vector<vector<double>> N) {
    // Body of the trans_plane_position_and_optical_path function translated earlier
}

vector<double> light_source_setting(double aperture, vector<double> distance, int cross_diameter_num, double ang_x, double ang_y) {
    // Body of the light_source_setting function translated earlier
}

double diffraction_limit(double lambda, double aperture, double BFL, double EFL, double focal_plane_position) {
    // Body of the diffraction_limit function translated earlier
}

// main function
int main() {
    // Add the main logic of your program here, using the functions defined above.
    // This is where you would call the above functions in the appropriate order, pass the necessary arguments to them,
    // and handle their return values, similar to how it is done in the MATLAB main script. However, we would need to adjust the logic to fit C++ syntax and semantics.

    return 0;
}
