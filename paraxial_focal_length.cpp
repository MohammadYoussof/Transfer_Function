 // paraxial_focal_length calculates the back focal length (BFL) and effective focal length (EFL) given the parameters surface_num, distance, material, and sur_radius.
 // The equation is taken from "Introdution to Optics Third Edition" by F.L. Pedrotti et al. (2006)
 // It returns two values: BFL and EFL. 

 #include <vector>

using namespace std;

pair<double, double> paraxial_focal_length(int surface_num, vector<double> distance, vector<double> material, vector<double> sur_radius) {
    vector<vector<double>> M = {{1, 0}, {0, 1}};

    for (int i = 0; i < surface_num; i++) {
        vector<vector<double>> M_traslation = {{1, distance[i]}, {0, 1}};
        vector<vector<double>> M_refraction = {{1, 0}, {(material[i] - material[i + 1]) / (material[i + 1] * sur_radius[i]), material[i] / material[i + 1]}};

        vector<vector<double>> temp(2, vector<double>(2));
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int l = 0; l < 2; l++)
                    temp[j][k] += M_refraction[j][l] * M_traslation[l][k];

        M = temp;
    }

    double A = M[0][0];
    double B = M[0][1];
    double C = M[1][0];
    double D = M[1][1];
    double q = -A / C;
    double s = (1 - A) / C;
    double f_s = q - s;

    return {q, f_s};
}
