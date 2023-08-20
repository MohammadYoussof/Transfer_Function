// line_spread_function is calculated given the parameters lambda, aperture, trans_plane_data, and focal_plane_position.
// It returns a structure containing the results.

#include <vector>
#include <complex>
#include <cmath>
#include <numeric>

using namespace std;

struct TransPlaneData {
    vector<double> y, OP;
};

struct LSFData {
    vector<double> Monitor_y, Intensity, Intensity_normalize, Power, Power_normalize;
    double Strehl_ratio;
};

LSFData line_spread_function(double lambda, double aperture, TransPlaneData trans_plane_data, double focal_plane_position) {
    double k0 = 2 * M_PI / lambda;
    int center_index = (trans_plane_data.y.size() + 1) / 2;

    for (int i = 0; i < trans_plane_data.OP.size(); i++) {
        if (std::isnan(trans_plane_data.OP[i]))
            trans_plane_data.OP[i] = 0;
        if (std::isnan(trans_plane_data.y[i]))
            trans_plane_data.y[i] = 0;
    }

    vector<complex<double>> phase_Mask(trans_plane_data.OP.size());
    for (int i = 0; i < trans_plane_data.OP.size(); i++) {
        double amp = trans_plane_data.OP[i];
        if (amp != 0)
            amp = 1;
        double phase = k0 * trans_plane_data.OP[i];
        phase_Mask[i] = amp * exp(complex<double>(0, 1) * phase);
    }

    double aperture_radius = aperture / 2;
    vector<double> Monitor_y(trans_plane_data.y.size());
    double linsp = -aperture_radius;
    double step = aperture / (Monitor_y.size() - 1);
    for (int i = 0; i < Monitor_y.size(); i++)
        Monitor_y[i] = linsp + i * step;

    double Distance = focal_plane_position;
    vector<complex<double>> Monitor(Monitor_y.size());
    for (int i = 0; i < Monitor_y.size(); i++) {
        for (int j = 0; j < trans_plane_data.y.size(); j++) {
            double R = sqrt(pow(Distance, 2) + pow(Monitor_y[i] - trans_plane_data.y[j], 2));
            Monitor[i] += phase_Mask[j] * exp(complex<double>(0, 1) * k0 * R);
        }
    }

    LSFData result;
    result.Monitor_y = Monitor_y;
    for (complex<double>& val : Monitor) {
        result.Intensity.push_back(abs(val));
        result.Power.push_back(pow(abs(val), 2));
    }
    double amp_sum = accumulate(trans_plane_data.OP.begin(), trans_plane_data.OP.end(), 0.0);
    for (double& val : result.Intensity)
        result.Intensity_normalize.push_back(val / amp_sum);
    for (double& val : result.Power)
        result.Power_normalize.push_back(val / pow(amp_sum, 2));
    result.Strehl_ratio = *max_element(result.Power.begin(), result.Power.end()) / pow(amp_sum, 2);

    return result;
}
