#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#define GLM_FORCE_CTOR_INIT
#include <glm/glm.hpp>

using namespace glm; // This allows the GLSL code to be used directly

#include "brdf_implementations.glsl"

class Timer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::duration<double, std::milli> accumulated_time = std::chrono::milliseconds(0);
    std::string name;

public:
    // Constructor to set the name of the timer for identification
    Timer(const std::string& name) : name(name) {}

    // Start the timer
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    // Stop the timer and accumulate the time elapsed since the last start
    void stop() {
        auto end_time = std::chrono::high_resolution_clock::now();
        accumulated_time += end_time - start_time;
    }

    // Reset the timer to zero
    void reset() {
        accumulated_time = std::chrono::milliseconds(0);
    }

    // Get the total accumulated time
    double getAccumulatedTime() const {
        return accumulated_time.count();
    }

    // Display the accumulated time
    void printAccumulatedTime() const {
        std::cout << name << " total time: " << getAccumulatedTime() << " milliseconds" << std::endl;
    }

    // Display the accumulated time adjusted for a baseline accumulated time
    void printAccumulatedTimeAdjustedForBaseline(double baseline) const {
        std::cout << name << " total time: " << (getAccumulatedTime() - baseline) << " milliseconds" << std::endl;
    }

    void printTimePerInvocation(unsigned int sample_count) const {
        std::cout << name << " time per invocation: " << getAccumulatedTime()*1.0e6f/float(sample_count) << " nanoseconds" << std::endl;
    }

};

vec3 calculate_wi_local(float theta_i) {
    return vec3(std::sin(radians(theta_i)), 0.0f, std::cos(radians(theta_i)));
}

vec3 calculate_wo_local(float theta_o, float phi) {
    return vec3(std::sin(radians(theta_o)) * std::cos(radians(phi)),
                std::sin(radians(theta_o)) * std::sin(radians(phi)),
                std::cos(radians(theta_o)));
}

// Generate a random float in [0.0f, 1.0f)
float random_float() {
    // Initialize a random engine and uniform distribution
    constexpr unsigned int seed = 12345; // Seed with a constant for reproducibility
    static std::mt19937 generator(seed);
    static std::uniform_real_distribution<float> distribution(0.0f, std::nextafter(1.0f, 0.0f)); // Max value is next to 1.0f but not including 1.0f

    return distribution(generator);
}

// Function to sample the upper z hemisphere with a cosine distribution
vec3 sample_unit_hemisphere_cosine() {
    constexpr float TwoPi = 6.28318531f;

    float phi = TwoPi * random_float();
    float z = sqrt(random_float()); // cos(theta)
    float r = sqrt(1.0f - z * z);   // sin(theta)
    return vec3(cos(phi) * r, sin(phi) * r, z);
}