// This function calculate PSF given the parameters lambda, aperture, trans_plane_data, propagate_distance, (and maybe perspective_switch).
// It returns a structure containing the results.
// This translation includes only the first part of the function, corresponding to when perspective_switch[0] == 1.

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;

struct TransPlaneData {
    vector<double> x, y, OP, dz;
};

struct PSFData {
    vector<double> Monitor1_y, Monitor_z;
    vector<vector<double>> Intensity, Intensity_normalize, Power, Power_normalize;
    double Strehl_ratio;
};

PSFData point_spread_function(double lambda, double aperture, TransPlaneData trans_plane_data, double propagate_distance, vector<int> perspective_switch) {
    double k0 = 2 * M_PI / lambda;

    vector<double> amp = trans_plane_data.OP;
    for (double& val : amp)
        val = std::isnan(val) ? 0 : 1;

    for (double& val : trans_plane_data.x)
        if (std::isnan(val))
            val = 0;

    for (double& val : trans_plane_data.y)
        if (std::isnan(val))
            val = 0;

    for (double& val : trans_plane_data.OP)
        if (std::isnan(val))
            val = 0;

    PSFData result;

    if (perspective_switch[0] == 1) {
        vector<double> y_bound = {-aperture / 2, aperture / 2};
        vector<vector<double>> Monitor_Boundary = {
            {y_bound[0], y_bound[1], (double)trans_plane_data.y.size()},
            {(propagate_distance + propagate_distance * 0.3) / 200, propagate_distance + propagate_distance * 0.3, (propagate_distance + propagate_distance * 0.3) / 200}
        };
        vector<double> Monitor1_y = linspace(Monitor_Boundary[0][0], Monitor_Boundary[0][1], Monitor_Boundary[0][2]);
        vector<double> Monitor1_z = linspace(Monitor_Boundary[1][0], Monitor_Boundary[1][1], Monitor_Boundary[1][2]);

        vector<complex<double>> phase_mask;
        for (double val : amp)
            phase_mask.push_back(val * std::exp(complex<double>(0, 1) * k0 * val));

        vector<vector<double>> Monitor(Monitor1_y.size(), vector<double>(Monitor1_z.size()));
        for (int i = 0; i < Monitor1_z.size(); i++) {
            vector<double> Monitor_xy(phase_mask.size());
            for (int j = 0; j < Monitor1_y.size(); j++) {
                vector<double> R(phase_mask.size());
                for (int k = 0; k < R.size(); k++)
                    R[k] = sqrt(pow(Monitor1_z[i], 2) + pow(trans_plane_data.y[k] - Monitor1_y[j], 2) + pow(trans_plane_data.x[k], 2));
                for (int k = 0; k < R.size(); k++)
                    Monitor_xy[j] += phase_mask[k] * std::exp(complex<double>(0, 1) * k0 * R[k]);
            }
            for (int j = 0; j < Monitor.size(); j++)
                Monitor[j][i] = abs(Monitor_xy[j]);
        }

        int index_z = 0;
        double min_val = std::numeric_limits<double>::max();
        for (int i = 0; i < Monitor1_z.size(); i++) {
            double val = abs(Monitor1_z[i] - trans_plane_data.dz[0] - propagate_distance);
            if (val < min_val) {
                min_val = val;
                index_z = i;
            }
        }

        result.Monitor1_y = Monitor1_y;
        result.Monitor_z = Monitor1_z;
        for (auto& row : Monitor) {
            result.Intensity.push_back(abs(row));
            result.Power.push_back(pow(abs(row), 2));
        }
        double amp_sum = std::accumulate(amp.begin(), amp.end(), 0.0);
        for (auto& row : result.Intensity) {
            for (double& val : row)
                val /= amp_sum;
        }
        for (auto& row : result.Power) {
            for (double& val : row)
                val /= pow(amp_sum, 2);
        }
        result.Strehl_ratio = *max_element(result.Power[index_z].begin(), result.Power[index_z].end()) / pow(amp_sum, 2);
    }

    // The second part of the function which corresponds to perspective_switch[1] == 1 is not translated here
    // This part of the function would require a similar process to the first part, but with different calculations

    return result;
}