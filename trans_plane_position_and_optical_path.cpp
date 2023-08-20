// This to calculate the position of the transmission plane and the optical path given the parameters surface_num, distance, material, Data, L, M, and N.
// It returns a structure containing the results.

#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

struct Data {
    vector<vector<vector<double>>> X_2, Y_2, Z_2;
};

struct TransPlaneData {
    vector<vector<double>> x, y, z, OP;
    vector<double> dz;
};

TransPlaneData trans_plane_position_and_optical_path(int surface_num, vector<double> distance, vector<double> material, Data data, vector<vector<double>> L, vector<vector<double>> M, vector<vector<double>> N) {
    vector<vector<double>> OP(data.X_2[0].size(), vector<double>(data.X_2[0][0].size()));
    for (int i = 0; i < surface_num; i++)
        for (int j = 0; j < data.X_2[0].size(); j++)
            for (int k = 0; k < data.X_2[0][0].size(); k++)
                OP[j][k] += sqrt(pow(data.X_2[1][i][j] - data.X_2[0][i][j], 2) + pow(data.Y_2[1][i][j] - data.Y_2[0][i][j], 2) + pow(data.Z_2[1][i][j] - data.Z_2[0][i][j], 2)) * material[i];

    vector<vector<double>> x0 = data.X_2[1][surface_num - 1];
    vector<vector<double>> y0 = data.Y_2[1][surface_num - 1];
    vector<vector<double>> z0 = data.Z_2[1][surface_num - 1];
    double max_z0 = -numeric_limits<double>::max();
    for (const auto& row : z0)
        for (double val : row)
            max_z0 = max(max_z0, val);
    vector<vector<double>> dz(x0.size(), vector<double>(x0[0].size()));
    for (auto& row : dz)
        for (double& val : row)
            val = max_z0 - val;

    vector<vector<double>> x1(x0.size(), vector<double>(x0[0].size()));
    vector<vector<double>> y1(y0.size(), vector<double>(y0[0].size()));
    vector<vector<double>> z1(z0.size(), vector<double>(z0[0].size(), max_z0));
    for (int i = 0; i < x0.size(); i++)
        for (int j = 0; j < x0[0].size(); j++) {
            x1[i][j] = x0[i][j] + (L[i][j] / N[i][j]) * dz[i][j];
            y1[i][j] = y0[i][j] + (M[i][j] / N[i][j]) * dz[i][j];
            OP[i][j] += sqrt(pow(x1[i][j] - x0[i][j], 2) + pow(y1[i][j] - y0[i][j], 2) + pow(z1[i][j] - z0[i][j], 2)) * material[surface_num];
        }

    vector<double> dz_sum(dz[0].size());
    for (double val : distance)
        for (double& val2 : dz_sum)
            val2 -= val;

    return {x1, y1, z1, dz_sum, OP};
}