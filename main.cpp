#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <math.h>

#define PI 3.14159265

using ScalarField = std::vector<std::vector<float> >;
using VectorField = std::vector<std::vector<std::vector<float> > >;

// meta params (non-physical simulation)
const int iters = 100;
const int jacobiIters = 10;
const float timeStep = 0.01;
const float epsilon = 10e-3;

// physical params
const float density = 1.0;

std::vector<float> bilinear(
        const VectorField& grid,
        const std::pair<float, float>& pos) {
    int height = grid.size();
    int width = grid[0].size();

    float x = std::clamp(pos.first, float(0.), float(width) - 1);
    float y = std::clamp(pos.second, float(0.), float(height) - 1);

    int x0 = std::clamp(int(x), 0, width - 1);
    int y0 = std::clamp(int(y), 0, height - 1);
    int x1 = std::clamp(x0 + 1, 0, width - 1);
    int y1 = std::clamp(y0 + 1, 0, height - 1);

    float xPercentile = (x - x0) / float(x1 - x0);
    float yPercentile = (y - y0) / float(y1 - y0);

    auto v00 = grid[y0][x0];
    auto v01 = grid[y0][x1];
    auto v10 = grid[y1][x0];
    auto v11 = grid[y1][x1];

    std::vector<float> result;
    for (int i = 0; i < v00.size(); i++) {
        float vb = v00[i] * (1 - xPercentile) + v01[i] * xPercentile;
        float vt = v10[i] * (1 - xPercentile) + v11[i] * xPercentile;
        float vf = vb * (1 - yPercentile) + vt * yPercentile;
        result.push_back(vf);
    }
    return result;
}

ScalarField calcDiv(const VectorField&  advectVelocity) {
    ScalarField divergence;
    int height = advectVelocity.size();
    int width = advectVelocity[0].size();

    for (int row = 0; row < height; row++) {
        std::vector<float> divRow;
        for (int col = 0; col < width; col++) {
            float div = -2.0 * epsilon * density / timeStep * (
                advectVelocity[row][(col + 1) % width][0] -
                advectVelocity[row][((col - 1) + width) % width][0] +
                advectVelocity[(row + 1) % height][col][1] -
                advectVelocity[((row - 1) + height) % height][col][1]
            );
            divRow.push_back(div);
        }
        divergence.push_back(divRow);
    }
    return divergence;
}

ScalarField calcPressure(const ScalarField& divergence) {
    ScalarField pressure;
    int height = divergence.size();
    int width = divergence[0].size();

    std::cout << height << " " << width << std::endl;

    for (int row = 0; row < height; row++) {
        std::vector<float> pressureRow(width, 0.0);
        pressure.push_back(pressureRow);
    }

    for (int jacobiIter = 0; jacobiIter < jacobiIters; jacobiIter++) {
        ScalarField updatedPressure;
        for (int row = 0; row < height; row++) {
            std::vector<float> pressureRow;
            for (int col = 0; col < width; col++) {
                float p = divergence[row][col] + 1.0 / 4.0 * (
                    pressure[(row + 2) % height][col] +
                    pressure[((row - 2) + height) % height][col] +
                    pressure[row][(col + 2) % width] +
                    pressure[row][((col - 2) + width) % width]
                );
                pressureRow.push_back(p);
            }
            updatedPressure.push_back(pressureRow);
        }
        pressure = updatedPressure;
    }
    return pressure;
}

int main() {
    const int width = 640;
    const int height = 480;

    VectorField color;
    VectorField velocity;
    for (int row = 0; row < height; row++) {
        std::vector<std::vector<float> > velRow;
        std::vector<std::vector<float> > colRow;

        for (int col = 0; col < width; col++) {
            velRow.push_back({
                static_cast<float>(sin(2 * PI * row)),
                static_cast<float>(sin(2 * PI * col))
            });
            colRow.push_back({
                static_cast<float>(abs(sin(2 * PI * row))),
                static_cast<float>(abs(cos(2 * PI * col))),
                static_cast<float>(abs(cos(2 * PI * col)))
            });
        }

        velocity.push_back(velRow);
        color.push_back(colRow);
    }

    for (int iter = 0; iter < iters; iter++) {
        VectorField advectVelocity;
        for (int row = 0; row < height; row++) {
            std::vector<std::vector<float> > advRow;
            for (int col = 0; col < width; col++) {
                std::pair<float, float> pos(col - velocity[row][col][0] * timeStep, row - velocity[row][col][1] * timeStep);
                advRow.push_back(bilinear(velocity, pos));
            }
            advectVelocity.push_back(advRow);
        }

        auto divergence = calcDiv(advectVelocity);
        auto pressure = calcPressure(divergence);

        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                float coeff = timeStep / (2.0 * density * epsilon);
                velocity[row][col] = {
                    advectVelocity[row][col][0] - coeff * (pressure[row][(col + 1) % width] - pressure[row][((col - 1) + width) % width]),
                    advectVelocity[row][col][1] - coeff * (pressure[(row + 1) % height][col] - pressure[((row - 1) + height) % height][col]),
                };

                std::pair<float, float> pos = {col - velocity[row][col][0] * timeStep, row - velocity[row][col][1] * timeStep};
                color[row][col] = bilinear(color, pos);
            }
        }
    }

    return 0;
}
