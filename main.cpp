#include <iostream
#include <vector>
#include <string>
#include <math.h>

#define PI 3.14159265

int main() {
    std::vector<std::vector<std::pair<float, float, float>>> color;
    std::vector<std::vector<std::pair<float, float>>> velocity;

    const int width = 640;
    const int height = 480;

    for (int row = 0; row < height; row++) {
        std::vector<float> row;
        for (int col = 0; col < width; col++) {
            row.push_back(std::pair(sin(2 * PI * row), sin(2 * PI * col)));
        }
        velocity.push_back(row);
    }

    // meta params (non-physical simulation)
    int iters = 100;
    float timeStep = 0.01;

    // physical params
    float density = 1.0;

    for (int iter = 0; iter < iters; iter++) {
        std::vector<std::vector<std::pair<float, float>>> advectVelocity;
        for (int row = 0; row < height; row++) {
            std::vector<float> row;
            for (int col = 0; col < width; col++) {
                row.push_back(std::pair(sin(2 * PI * row), sin(2 * PI * col)));
            }
        }
    }

    return 0;
}
