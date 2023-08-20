// This function, light_source_setting, is to set up a light source given the parameters aperture, distance, cross_diameter_num, ang_x, and ang_y.
// It returns six arrays/matrices: s_x, s_y, s_z, L, M, and N.

#include <vector>
#include <cmath>

using namespace std;

vector<vector<double>> meshgrid(const vector<double>& a, const vector<double>& b) {
    vector<vector<double>> result(a.size(), vector<double>(b.size()));
    for (int i = 0; i < a.size(); i++)
        for (int j = 0; j < b.size(); j++)
            result[i][j] = a[i];
    return result;
}

tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>>
light_source_setting(double aperture, vector<double> distance, int cross_diameter_num, double ang_x, double ang_y) {
    if (cross_diameter_num % 2 == 0)
        cross_diameter_num += 1;

    vector<double> x(cross_diameter_num), y(cross_diameter_num);
    double linsp = -aperture / 2;
    double step = aperture / (cross_diameter_num - 1);
    for (int i = 0; i < cross_diameter_num; i++) {
        x[i] = linsp + i * step;
        y[i] = linsp + i * step;
    }

    vector<vector<double>> s_x = meshgrid(x, y);
    vector<vector<double>> s_y = meshgrid(y, x);
    for (int i = 0; i < cross_diameter_num; i++)
        for (int j = 0; j < cross_diameter_num; j++) {
            double r = sqrt(pow(s_x[i][j], 2) + pow(s_y[i][j], 2));
            if (r > aperture / 2) {
                s_x[i][j] = numeric_limits<double>::quiet_NaN();
                s_y[i][j] = numeric_limits<double>::quiet_NaN();
            }
        }

    vector<vector<double>> s_z(cross_diameter_num, vector<double>(cross_diameter_num, 0));
    vector<vector<double>> L(cross_diameter_num, vector<double>(cross_diameter_num, sin(ang_x * M_PI / 180)));
    vector<vector<double>> M(cross_diameter_num, vector<double>(cross_diameter_num, sin(ang_y * M_PI / 180)));
    vector<vector<double>> N(cross_diameter_num, vector<double>(cross_diameter_num));
    for (int i = 0; i < cross_diameter_num; i++)
        for (int j = 0; j < cross_diameter_num; j++)
            N[i][j] = sqrt(1 - (pow(L[i][j], 2) + pow(M[i][j], 2)));

    if (ang_x != 0 || ang_y != 0)
        for (int i = 0; i < cross_diameter_num; i++)
            for (int j = 0; j < cross_diameter_num; j++) {
                s_x[i][j] = s_x[i][j] + (L[i][j] / N[i][j]) * (-distance[0]);
                s_y[i][j] = s_y[i][j] + (M[i][j] / N[i][j]) * (-distance[0]);
            }

    return make_tuple(s_x, s_y, s_z, L, M, N);
}