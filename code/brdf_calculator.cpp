/*
Prerequisites:

- Clone the glm repository from https://github.com/g-truc/glm
- Set the path to the glm repository in the include directive below

General compilation and running commands:

g++ -std=c++17 -O3 -o brdf_calculator brdf_calculator.cpp -I/path/to/glm
./brdf_calculator

Compilation and running commands for Peter's setup:

cd ~/repos/EON_jcgt_paper/code
g++ -std=c++17 -O3 -o brdf_calculator brdf_calculator.cpp -I/Users/pkutz/repos/glm/
./brdf_calculator

Compilation and running script for Peter's setup (modify for your own setup):

./plot_brdfs.sh
*/

#include "brdf_common.h"

#include <fstream>
#include <string>
#include <vector>

class BRDFCalculator {
public:
    void computeAndOutputBRDFData() {
        vec3 rho = vec3(1.f);  // Albedo.

        std::vector<float> roughnesses = {0.1f, 0.666f}; // Roughness values.

        for (float roughness : roughnesses)
        {
            std::vector<float> thetas_i = {45.0f, 80.0f}; // Incident zenith angles in degrees.
            std::vector<float> phis = {0.0f, 45.0f};     // Relative azimuth angles in degrees.

            // Loop through each combination of theta_i and phi to generate separate CSV files.
            for (float theta_i : thetas_i) {
                for (float phi_o : phis) {
                    std::string filename = "brdf_theta_i_" + std::to_string((int)theta_i) + "_phi_" + std::to_string((int)phi_o) + "_roughness_" + std::to_string((roughness)) + ".csv";
                    std::ofstream file(filename);

                    // Write headers to the CSV file.
                    file << "theta_o,Lambert,ON,QON,footnoteQON,fujiiQON,FON,EON,GGX,MRM\n";

                    float angle_tol = 1.0e-4f;
                    // Compute BRDF values over the range of theta_o.
                    for (float theta_o = -90.0f+angle_tol; theta_o <= 90.0f-angle_tol; theta_o += 1.0f) {
                        vec3 wi_local = calculate_wi_local(theta_i);
                        vec3 wo_local = calculate_wo_local(theta_o, phi_o);

                        // Call each BRDF function with the generated directions.
                        vec3 brdfValue1 = f_lambert(rho);
                        vec3 brdfValue2 = f_ON(rho, roughness, wi_local, wo_local);
                        vec3 brdfValue3 = f_QON(rho, roughness, wi_local, wo_local, false);
                        vec3 brdfValue4 = f_QON(rho, roughness, wi_local, wo_local, true);
                        vec3 brdfValue5 = f_fujiiQON(rho, roughness, wi_local, wo_local);
                        vec3 brdfValue6 = f_FON(rho, roughness, wi_local, wo_local);
                        vec3 brdfValue7 = f_EON(rho, roughness, wi_local, wo_local, false);
                        vec3 brdfValue8 = f_GGX(rho, roughness, wi_local, wo_local);
                        vec3 brdfValue9 = f_MRM(rho, roughness, wi_local, wo_local);

                        // Extract the first element from each color, assuming all components are identical.
                        float scalarBrdfValue1 = brdfValue1[0];
                        float scalarBrdfValue2 = brdfValue2[0];
                        float scalarBrdfValue3 = brdfValue3[0];
                        float scalarBrdfValue4 = brdfValue4[0];
                        float scalarBrdfValue5 = brdfValue5[0];
                        float scalarBrdfValue6 = brdfValue6[0];
                        float scalarBrdfValue7 = brdfValue7[0];
                        float scalarBrdfValue8 = brdfValue8[0];
                        float scalarBrdfValue9 = brdfValue9[0];

                        // Write the BRDF values to the CSV file.
                        file << theta_o << ","
                            << scalarBrdfValue1 << ","
                            << scalarBrdfValue2 << ","
                            << scalarBrdfValue3 << ","
                            << scalarBrdfValue4 << ","
                            << scalarBrdfValue5 << ","
                            << scalarBrdfValue6 << ","
                            << scalarBrdfValue7 << ","
                            << scalarBrdfValue8 << ","
                            << scalarBrdfValue9 << "\n";
                    }
                    file.close();
                }
            }

        }

    }
};

int main() {
    BRDFCalculator calculator;
    calculator.computeAndOutputBRDFData();
    return 0;
}
