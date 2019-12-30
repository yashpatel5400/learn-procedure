#include <iostream
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
    int x0 = int(pos[0]);
    int y0 = int(pos[1]);
    int x1 = x0 + 1;
    int y1 = y0 + 1;

    float xPercentile = (pos[0] - x0) / float(x1 - x0);
    float yPercentile = (pos[1] - y0) / float(y1 - y0);

    float v00 = grid[y0][x0];
    float v01 = grid[y0][x1];
    float v10 = grid[y1][x0];
    float v11 = grid[y1][x1];

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
    for (int row = 0; row < height; row++) {
        std::vector<float> row;
        for (int col = 0; col < width; col++) {
            float div = -2.0 * epsilon * density / timeStep * (
                advectVelocity[col + 1][row] -
                advectVelocity[col - 1][row] +
                advectVelocity[col][row + 1] -
                advectVelocity[col][row - 1]
            );
            row.push_back(div);
        }
        divergence.push_back(row);
    }
    return divergence;
}

ScalarField calcPressure(const ScalarField& divergence) {
    ScalarField pressure;
    for (int row = 0; row < height; row++) {
        std::vector<float> row(width, 0.0);
        pressure.push_back(row);
    }

    for (int jacobiIter = 0; jacobiIter < jacobiIters; jacobiIter++) {
        std::vector<float> updatedPressure;
        for (int row = 0; row < height; row++) {
            std::vector<float> row;
            for (int col = 0; col < width; col++) {
                float p = divergence[row][col] + 1.0 / 4.0 * (
                    pressure[row + 2][col] +
                    pressure[row - 2][col] +
                    pressure[row][col + 2] +
                    pressure[row][col - 2]
                );
                row.push_back(p);
            }
            updatedPressure.push_back(row);
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
        std::vector<float> velRow;
        std::vector<float> colRow;

        for (int col = 0; col < width; col++) {
            velRow.emplace_back({sin(2 * PI * row), sin(2 * PI * col)});
            colRow.emplace_back({abs(sin(2 * PI * row)), abs(cos(2 * PI * col)), abs(cos(2 * PI * col))});
        }

        velocity.push_back(velRow);
        color.push_back(colRow);
    }

    for (int iter = 0; iter < iters; iter++) {
        VectorField advectVelocity;
        for (int row = 0; row < height; row++) {
            std::vector<float> row;
            for (int col = 0; col < width; col++) {
                std::pair<float, float> pos(col - velocity[row][col][0] * timeStep, row - velocity[row][col][1] * timeStep);
                row.push_back(bilinear(velocity, pos));
            }
            advectVelocity.push_back(row);
        }

        auto divergence = calcDiv(advectVelocity);
        auto pressure = calcPressure(divergence);

        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                float update = timeStep / (2.0 * density * epsilon) * (pressure[row][col + 1] - pressure[row][col - 1]);
                velocity[row][col] = {
                    advectVelocity[row][col][0] - update,
                    advectVelocity[row][col][1] - update,
                };

                std::pair<float, float> pos = (col - velocity[row][col][0] * timeStep, row - velocity[row][col][1] * timeStep);
                color[row][col] = bilinear(color, pos);
            }
        }
    }

    return 0;
}
