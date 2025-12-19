/*
Prerequisites:

- Clone the glm repository from https://github.com/g-truc/glm
- Set the path to the glm repository in the include directive below

General compilation and running commands:

g++ -std=c++17 -O3 -o brdf_profiler brdf_profiler.cpp -I/path/to/glm
./brdf_profiler

Compilation and running commands for Peter's setup:

cd ~/repos/EON_jcgt_paper/code
g++ -std=c++17 -O3 -o brdf_profiler brdf_profiler.cpp -I/Users/pkutz/repos/glm/
./brdf_profiler

Compilation and running script for Peter's setup (modify for your own setup):

./profile_brdfs.sh
*/

#include "brdf_common.h"

#include <algorithm>
#include <chrono>
#include <random>

class BRDFProfiler {
public:
    // This function profiles all of the BRDF functions.
    // Details:
    //     - It runs on a single thread.
    //     - It profiles only the BRDF evaluation (no importance sampling).
    //     - It interleaves different BRDFs to avoid bias due to system load varying over time.
    //     - It varies all BRDF inputs to make sure different performance characteristics are covered.
    //     - It calculates random directions outside the timed loop to avoid direction sampling overhead. 
    //     - It compares to a no-op baseline to compensate for timing and looping overhead.
    //     - It should be run with realistic compiler optimizations enabled.
    //     - It takes various measures to avoid code disappearing due to compiler optimizations (important).
    void profileBRDFs() {
        // Create timer instances for each BRDF.
        Timer noOpTimer("No-op"); // Baseline reference.
        Timer lambertTimer("Lambert");
        Timer fullOnTimer("Full ON");
        Timer qonTimer("QON");
        Timer footnoteQonTimer("Footnote QON");
        Timer fujiiQonTimer("Fujii QON");
        Timer fonTimer("FON");
        Timer eonExactTimer("EON (exact)");
        Timer eonApproxTimer("EON (approx)");

        // Accumulate and print the total BRDF value to avoid optimization.
        vec3 totalBrdfValue = vec3(0.0f);

        int numDirectionSamples = 40000; // Number of samples for profiling.

        for (int i = 0; i < numDirectionSamples; ++i) {
            // These direction are sampled outside the timer segments because they are expensive.
            vec3 wi = sample_unit_hemisphere_cosine();
            vec3 wo = sample_unit_hemisphere_cosine();

            // Define inner loops for the BRDF profiling.
            // The purpose is to do more within each timer segment to avoid timer overhead.
            // After each iteration we swap the directions for the next iteration to avoid optimization of lines involving the directions.
            vec3 rho;
            float roughness;
            int numRhoAndRSamples = 200;
#define BEGIN_INNER_LOOPS \
            for (int i = 0; i < numRhoAndRSamples; ++i) { \
                for (int j = 0; j < numRhoAndRSamples; ++j) { \
                    rho = vec3(float(i) * float(1.0f / numRhoAndRSamples)); \
                    roughness = float(j) * float(1.0f / numRhoAndRSamples);
#define END_INNER_LOOPS \
                    std::swap(wi, wo); /* Swap the direction to avoid optimization. */ \
                } \
            }

            // Accumulate the BRDF values to avoid optimization.
            vec3 brdfValue0 = vec3(0.0f);
            vec3 brdfValue1 = vec3(0.0f);
            vec3 brdfValue2 = vec3(0.0f);
            vec3 brdfValue3 = vec3(0.0f);
            vec3 brdfValue4 = vec3(0.0f);
            vec3 brdfValue5 = vec3(0.0f);
            vec3 brdfValue6 = vec3(0.0f);
            vec3 brdfValue7 = vec3(0.0f);
            vec3 brdfValue8 = vec3(0.0f);

            noOpTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue0 += vec3(0.1f); // Accumulate a constant value to avoid optimization.
            END_INNER_LOOPS
            noOpTimer.stop();

            lambertTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue1 += f_lambert(rho);
            END_INNER_LOOPS
            lambertTimer.stop();

            fullOnTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue2 += f_ON(rho, roughness, wi, wo);
            END_INNER_LOOPS
            fullOnTimer.stop();

            qonTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue3 += f_QON(rho, roughness, wi, wo, false);
            END_INNER_LOOPS
            qonTimer.stop();

            footnoteQonTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue4 += f_QON(rho, roughness, wi, wo, true);
            END_INNER_LOOPS
            footnoteQonTimer.stop();

            fujiiQonTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue5 += f_fujiiQON(rho, roughness, wi, wo);
            END_INNER_LOOPS
            fujiiQonTimer.stop();

            fonTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue6 += f_FON(rho, roughness, wi, wo);
            END_INNER_LOOPS
            fonTimer.stop();

            eonExactTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue7 += f_EON(rho, roughness, wi, wo, true);
            END_INNER_LOOPS
            eonExactTimer.stop();

            eonApproxTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue8 += f_EON(rho, roughness, wi, wo, false);
            END_INNER_LOOPS
            eonApproxTimer.stop();

            // Accumulate the BRDF values to avoid optimization.
            totalBrdfValue += brdfValue0 + brdfValue1 + brdfValue2 + brdfValue3 + brdfValue4 + brdfValue5 + brdfValue6 + brdfValue7 + brdfValue8;
        }

        std::cout << "BRDF profiling results:" << std::endl;
        std::cout << std::endl;

        // Print all the accumulated times.
        noOpTimer.printAccumulatedTime();
        lambertTimer.printAccumulatedTime();
        fullOnTimer.printAccumulatedTime();
        qonTimer.printAccumulatedTime();
        footnoteQonTimer.printAccumulatedTime();
        fujiiQonTimer.printAccumulatedTime();
        fonTimer.printAccumulatedTime();
        eonExactTimer.printAccumulatedTime();
        eonApproxTimer.printAccumulatedTime();

        std::cout << std::endl;
        std::cout << "BRDF profiling results adjusted for no-op baseline:" << std::endl;
        std::cout << std::endl;

        // Print all the accumulated times adjusted for the baseline accumulated time from the no-op timer.
        noOpTimer.printAccumulatedTimeAdjustedForBaseline(noOpTimer.getAccumulatedTime());
        lambertTimer.printAccumulatedTimeAdjustedForBaseline(noOpTimer.getAccumulatedTime());
        fullOnTimer.printAccumulatedTimeAdjustedForBaseline(noOpTimer.getAccumulatedTime());
        qonTimer.printAccumulatedTimeAdjustedForBaseline(noOpTimer.getAccumulatedTime());
        footnoteQonTimer.printAccumulatedTimeAdjustedForBaseline(noOpTimer.getAccumulatedTime());
        fujiiQonTimer.printAccumulatedTimeAdjustedForBaseline(noOpTimer.getAccumulatedTime());
        fonTimer.printAccumulatedTimeAdjustedForBaseline(noOpTimer.getAccumulatedTime());
        eonExactTimer.printAccumulatedTimeAdjustedForBaseline(noOpTimer.getAccumulatedTime());
        eonApproxTimer.printAccumulatedTimeAdjustedForBaseline(noOpTimer.getAccumulatedTime());

        // Print the total BRDF value to avoid optimization.
        std::cout << std::endl;
        std::cout << "Total BRDF value: " << totalBrdfValue[0] << ", " << totalBrdfValue[1] << ", " << totalBrdfValue[2] << std::endl;
    }
};

int main() {
    BRDFProfiler profiler;
    profiler.profileBRDFs();
    return 0;
}
