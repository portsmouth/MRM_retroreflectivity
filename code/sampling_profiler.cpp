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

        Timer lambertTimer("Lambert eval");
        Timer lambertCosSampleTimer("Lambert sample (cos)");

        Timer fullOnTimer("Full ON eval");
        Timer fullOnCosampleTimer("Full ON sample (cos)");

        Timer qonTimer("QON eval");
        Timer qonCosSampleTimer("QON sample (cos)");

        Timer footnoteQonTimer("Footnote QON eval");
        Timer footnoteQonCosSampleTimer("Footnote QON sample (cos)");

        Timer fujiiQonTimer("Fujii QON eval");
        Timer fujiiQonCosSampleTimer("Fujii QON sample (cos)");

        Timer fonTimer("FON eval");
        Timer fonCosSampleTimer("FON sample (cos)");

        Timer eonApproxTimer("EON (approx) eval");
        Timer eonApproxCosSampleTimer("EON (approx) sample (cos)");
        Timer eonApproxCLTCSampleTimer("EON (approx) sample (CLTC)");

        Timer eonExactTimer("EON (exact) eval");
        Timer eonExactCosSampleTimer("EON (exact) sample (cos)");
        Timer eonExactCLTCSampleTimer("EON (exact) sample (CLTC)");

        // Accumulate and print the total BRDF value to avoid optimization.
        vec3 totalBrdfValue = vec3(0.0f);
        vec4 totalSampleValue = vec4(0.0f);
        vec3 totalSampleBrdfValue = vec4(0.0f);

        int numDirectionSamples = 4000; // Number of samples for profiling.
        int numRhoAndRSamples = 200;
        int total_sample_count = numDirectionSamples * numRhoAndRSamples * numRhoAndRSamples;

        for (int i = 0; i < numDirectionSamples; ++i) {

            if (i % (numDirectionSamples/100)  == 0)
                std::cout << i << "/" << numDirectionSamples << std::endl;
            // These direction are sampled outside the timer segments because they are expensive.
            vec3 wi = sample_unit_hemisphere_cosine();
            vec3 wo = sample_unit_hemisphere_cosine();

            volatile float u1 = random_float();
            volatile float u2 = random_float();

            // Define inner loops for the BRDF profiling.
            // The purpose is to do more within each timer segment to avoid timer overhead.
            // After each iteration we swap the directions for the next iteration to avoid optimization of lines involving the directions.
            vec3 rho;
            float roughness;
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
            vec3 sampleBrdfValue0 = vec3(0.0f);
            vec4 sampleValue0 = vec4(0.0f);

            vec3 brdfValue1 = vec3(0.0f);
            vec3 sampleBrdfValue1 = vec3(0.0f);
            vec4 sampleValue1 = vec4(0.0f);

            vec3 brdfValue2 = vec3(0.0f);
            vec3 sampleBrdfValue2 = vec3(0.0f);
            vec4 sampleValue2 = vec4(0.0f);

            vec3 brdfValue3 = vec3(0.0f);
            vec3 sampleBrdfValue3 = vec3(0.0f);
            vec4 sampleValue3 = vec4(0.0f);

            vec3 brdfValue4 = vec3(0.0f);
            vec3 sampleBrdfValue4 = vec3(0.0f);
            vec4 sampleValue4 = vec4(0.0f);

            vec3 brdfValue5 = vec3(0.0f);
            vec3 sampleBrdfValue5 = vec3(0.0f);
            vec4 sampleValue5 = vec4(0.0f);

            vec3 brdfValue6 = vec3(0.0f);
            vec3 sampleBrdfValue6 = vec3(0.0f);
            vec4 sampleValue6 = vec4(0.0f);

            vec3 brdfValue7 = vec3(0.0f);
            vec3 sampleBrdfValue7 = vec3(0.0f);
            vec4 sampleValue7 = vec4(0.0f);

            vec3 brdfValue8 = vec3(0.0f);
            vec3 sampleBrdfValue8 = vec3(0.0f);
            vec4 sampleValue8 = vec4(0.0f);

            vec3 sampleBrdfValue9 = vec3(0.0f);
            vec4 sampleValue9 = vec4(0.0f);

            vec3 sampleBrdfValue10 = vec3(0.0f);
            vec4 sampleValue10 = vec4(0.0f);

            vec3 sampleBrdfValue11 = vec3(0.0f);
            vec4 sampleValue11 = vec4(0.0f);

            noOpTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue0 += vec3(0.1f); // Accumulate a constant value to avoid optimization.
            END_INNER_LOOPS
            noOpTimer.stop();

            // Lambert /////////////////////////////////////////
            lambertTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue1 += f_lambert(rho);
            END_INNER_LOOPS
            lambertTimer.stop();

            lambertCosSampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue1 += sample_lambert(rho, wo, u1, u2, sampleBrdfValue1);
            END_INNER_LOOPS
            lambertCosSampleTimer.stop();

            // fullON /////////////////////////////////////////
            fullOnTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue2 += f_ON(rho, roughness, wi, wo);
            END_INNER_LOOPS
            fullOnTimer.stop();

            fullOnCosampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue2 += sample_ON(rho, roughness, wo, u1, u2, sampleBrdfValue2);
            END_INNER_LOOPS
            fullOnCosampleTimer.stop();

            // QON  /////////////////////////////////////////
            qonTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue3 += f_QON(rho, roughness, wi, wo, false);
            END_INNER_LOOPS
            qonTimer.stop();

            qonCosSampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue3 += sample_QON(rho, roughness, wo, false, u1, u2, sampleBrdfValue3);
            END_INNER_LOOPS
            qonCosSampleTimer.stop();

            // footnote QON  /////////////////////////////////////////
            footnoteQonTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue4 += f_QON(rho, roughness, wi, wo, true);
            END_INNER_LOOPS
            footnoteQonTimer.stop();

            footnoteQonCosSampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue4 += sample_QON(rho, roughness, wo, true, u1, u2, sampleBrdfValue4);
            END_INNER_LOOPS
            footnoteQonCosSampleTimer.stop();

            // fujii QON /////////////////////////////////////////
            fujiiQonTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue5 += f_fujiiQON(rho, roughness, wi, wo);
            END_INNER_LOOPS
            fujiiQonTimer.stop();

            fujiiQonCosSampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue5 += sample_fujiiQON(rho, roughness, wo, u1, u2, sampleBrdfValue5);
            END_INNER_LOOPS
            fujiiQonCosSampleTimer.stop();

            // FON /////////////////////////////////////////
            fonTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue6 += f_FON(rho, roughness, wi, wo);
            END_INNER_LOOPS
            fonTimer.stop();

            fonCosSampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue6 += sample_FON(rho, roughness, wo, u1, u2, sampleBrdfValue6);
            END_INNER_LOOPS
            fonCosSampleTimer.stop();

            // EON approx /////////////////////////////////////////
            eonApproxTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue7 += f_EON(rho, roughness, wi, wo, false);
            END_INNER_LOOPS
            eonApproxTimer.stop();

            eonApproxCosSampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue7 += sample_EON_COS(rho, roughness, wo, false, u1, u2, sampleBrdfValue7);
            END_INNER_LOOPS
            eonApproxCosSampleTimer.stop();

            eonApproxCLTCSampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue8 += sample_EON_CLTC(rho, roughness, wo, false, u1, u2, sampleBrdfValue8);
            END_INNER_LOOPS
            eonApproxCLTCSampleTimer.stop();

            // EON exact /////////////////////////////////////////
            eonExactTimer.start();
            BEGIN_INNER_LOOPS
            brdfValue8 += f_EON(rho, roughness, wi, wo, true);
            END_INNER_LOOPS
            eonExactTimer.stop();

            eonExactCosSampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue9 += sample_EON_COS(rho, roughness, wo, true, u1, u2, sampleBrdfValue9);
            END_INNER_LOOPS
            eonExactCosSampleTimer.stop();

            eonExactCLTCSampleTimer.start();
            BEGIN_INNER_LOOPS
            sampleValue10 += sample_EON_CLTC(rho, roughness, wo, true, u1, u2, sampleBrdfValue10);
            END_INNER_LOOPS
            eonExactCLTCSampleTimer.stop();

            // Accumulate the BRDF values to avoid optimization.
            totalBrdfValue += brdfValue0 + brdfValue1 + brdfValue2 + brdfValue3 + brdfValue4 + brdfValue5 + brdfValue6 + brdfValue7 + brdfValue8;
            totalSampleValue += sampleValue0 + sampleValue1 + sampleValue2 + sampleValue3 + sampleValue4 + sampleValue5 + sampleValue6 + sampleValue7 + sampleValue8 + sampleValue9 + sampleValue10;
            totalSampleBrdfValue += sampleBrdfValue0 + sampleBrdfValue1 + sampleBrdfValue2 + sampleBrdfValue3 + sampleBrdfValue4 + sampleBrdfValue5 + sampleBrdfValue6 + sampleBrdfValue7 + sampleBrdfValue8 + sampleBrdfValue9 + sampleBrdfValue10;
        }

        std::cout << "BRDF sampling profiling results:" << std::endl;
        std::cout << "Total invocation count: " << total_sample_count << std::endl;
        std::cout << std::endl;

        noOpTimer.printTimePerInvocation(total_sample_count);

        lambertTimer.printTimePerInvocation(total_sample_count);
        lambertCosSampleTimer.printTimePerInvocation(total_sample_count);

        fullOnTimer.printTimePerInvocation(total_sample_count);
        fullOnCosampleTimer.printTimePerInvocation(total_sample_count);

        qonTimer.printTimePerInvocation(total_sample_count);
        qonCosSampleTimer.printTimePerInvocation(total_sample_count);

        footnoteQonTimer.printTimePerInvocation(total_sample_count);
        footnoteQonCosSampleTimer.printTimePerInvocation(total_sample_count);

        fujiiQonTimer.printTimePerInvocation(total_sample_count);
        fujiiQonCosSampleTimer.printTimePerInvocation(total_sample_count);

        fonTimer.printTimePerInvocation(total_sample_count);
        fonCosSampleTimer.printTimePerInvocation(total_sample_count);

        eonApproxTimer.printTimePerInvocation(total_sample_count);
        eonApproxCosSampleTimer.printTimePerInvocation(total_sample_count);
        eonApproxCLTCSampleTimer.printTimePerInvocation(total_sample_count);

        eonExactTimer.printTimePerInvocation(total_sample_count);
        eonExactCosSampleTimer.printTimePerInvocation(total_sample_count);
        eonExactCLTCSampleTimer.printTimePerInvocation(total_sample_count);

        // Print the total BRDF value to avoid optimization.
        std::cout << std::endl;
        std::cout << "Total BRDF value: " << totalBrdfValue[0] << ", " << totalBrdfValue[1] << ", " << totalBrdfValue[2] << std::endl;
        std::cout << "Total sample value: " << totalSampleValue[0] << ", " << totalSampleValue[1] << ", " << totalSampleValue[2] << std::endl;
    }
};

int main() {
    BRDFProfiler profiler;
    profiler.profileBRDFs();
    return 0;
}
