#include <vector>
#include <unordered_map>

using namespace std;

struct Data {
    vector<vector<double>> X_1, Y_1, Z_1;
    vector<vector<vector<double>>> X_2, Y_2, Z_2;
};

Data data_reshape(vector<vector<double>> s_x_all, vector<vector<double>> s_y_all,
                  vector<vector<double>> s_z_all, int cross_diameter_num) {
    Data data;

    for (int i = 0; i < s_x_all.size(); i++) {
        vector<double> X_1_temp, Y_1_temp, Z_1_temp;
        for (int j = 0; j < cross_diameter_num; j++) {
            X_1_temp.push_back(s_x_all[i][j]);
            Y_1_temp.push_back(s_y_all[i][j]);
            Z_1_temp.push_back(s_z_all[i][j]);
        }

        vector<double> temp_x, temp_y, temp_z;
        for (int j = cross_diameter_num; j < s_x_all[i].size(); j++) {
            temp_x.push_back(s_x_all[i][j]);
            temp_y.push_back(s_y_all[i][j]);
            temp_z.push_back(s_z_all[i][j]);
        }

        data.X_1.push_back(X_1_temp);
        data.X_1.push_back(temp_x);
        data.Y_1.push_back(Y_1_temp);
        data.Y_1.push_back(temp_y);
        data.Z_1.push_back(Z_1_temp);
        data.Z_1.push_back(temp_z);

        data.X_2.push_back({X_1_temp, temp_x});
        data.Y_2.push_back({Y_1_temp, temp_y});
        data.Z_2.push_back({Z_1_temp, temp_z});
    }

    return data;
}

// The reshaping is done manually, which may not give the same results as MATLAB's reshape function if the sizes of the input arrays are not evenly divisible by cross_diameter_num.
//
// In C++, we often use the std::vector to represent dynamic arrays, and std::unordered_map for structures like MATLAB structs. 
// However, in this case, I've created a struct to make it easier to access the reshaped data.
//
// This function also assumes that s_x_all, s_y_all, and s_z_all are 2D arrays represented as vectors of vectors.
//
// No estimate for errors or edge effect cases. 
