#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <fftw3.h>

using namespace std;

struct DiffractionLimit {
    vector<double> I, MTF, f, y;
    double cutoff_freq;
};

DiffractionLimit diffraction_limit(double lambda, double aperture, double BFL, double EFL, double focal_plane_position) {
    double z = focal_plane_position + (EFL - BFL);
    int size = round(2 * pow(aperture, 2) / (lambda * z)) * 20;
    if (size % 2 == 0)
        size += 1;

    vector<double> y(size);
    double linsp = -aperture / 2;
    double step = aperture / size;
    for (int i = 0; i < size; i++)
        y[i] = linsp + i * step;

    double f_number = z / aperture;
    double cutoff_freq = 1 / (lambda * f_number);

    double theta = atan(aperture / (2 * z));
    double linsp_theta = -theta;
    double step_theta = 2 * theta / size;
    vector<double> theta_all(size), I(size);
    for (int i = 0; i < size; i++) {
        theta_all[i] = linsp_theta + i * step_theta;
        I[i] = pow(sin((aperture / 2) * 2 * M_PI * sin(theta_all[i]) / lambda) / ((aperture / 2) * 2 * M_PI * sin(theta_all[i]) / lambda), 2);
        if (isnan(I[i]))
            I[i] = 1;
    }

    // FFT calculations require FFTW.
    
    vector<double> MTF(size);
    fftw_complex* data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
    fftw_plan plan = fftw_plan_dft_1d(size, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
    for (int i = 0; i < size; i++) {
        data[i][0] = I[i];
        data[i][1] = 0;
    }
    fftw_execute(plan);
    for (int i = 0; i < size; i++)
        MTF[i] = sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]);
    double max_val = *max_element(MTF.begin(), MTF.end());
    for (double& val : MTF)
        val /= max_val;
    fftw_destroy_plan(plan);
    fftw_free(data);

    vector<double> f(size);
    double linsp_f = -size / 2;
    double step_f = size / aperture;
    for (int i = 0; i < size; i++)
        f[i] = linsp_f + i * step_f;

    return {I, MTF, f, y, cutoff_freq};
}

//
// the FFTW library is equivalent to the fft function in MATLAB.
// This function assumes that the trigonometric functions operate similarly in MATLAB and C++.  
