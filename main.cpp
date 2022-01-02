#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

static constexpr  double two_div_rand_max = 2.f/RAND_MAX;

struct fill_struct
{
    float value[8];
};

fill_struct frand() {
    return {
            static_cast<float>(std::rand() * two_div_rand_max - 1),
            static_cast<float>(std::rand() * two_div_rand_max - 1),
            static_cast<float>(std::rand() * two_div_rand_max - 1),
            static_cast<float>(std::rand() * two_div_rand_max - 1),
            static_cast<float>(std::rand() * two_div_rand_max - 1),
            static_cast<float>(std::rand() * two_div_rand_max - 1),
            static_cast<float>(std::rand() * two_div_rand_max - 1),
            static_cast<float>(std::rand() * two_div_rand_max - 1)
    };
}

struct Star {
    float px, py, pz;
    float vx, vy, vz;
    float mass;
    float padding;
};

std::array<Star, 48> stars;

void init() {
    for (int i = 0; i < 48; i++) {
        auto rand = frand();
        stars.at(i) = {
                rand.value[0],
                rand.value[1],
                rand.value[2],
                rand.value[3],
                rand.value[4],
                rand.value[5],
                rand.value[6],
                rand.value[7]
        };
    }
#pragma GCC unroll 8
    for(int i=0; i<48; ++i)
    {
        stars.at(i).mass +=1;
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;
auto eps2 = eps*eps;
void step() {
    for (auto &star: stars) {
        for (auto &other: stars) {
            float dx = other.px - star.px;
            float dy = other.py - star.py;
            float dz = other.pz - star.pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            auto dt_div_d2 = dt / d2;
            star.vx += dx * other.mass * G * dt_div_d2;
            star.vy += dy * other.mass * G * dt_div_d2;
            star.vz += dz * other.mass * G * dt_div_d2;
        }
    }
#pragma omp simd
    for (auto &star: stars) {
        star.px += star.vx * dt;
        star.py += star.vy * dt;
        star.pz += star.vz * dt;
    }
}

float calc() {
    float energy = 0;
    for (auto &star: stars) {
        float v2 = star.vx * star.vx + star.vy * star.vy + star.vz * star.vz;
        energy += star.mass * v2 * 0.5;
        for (auto &other: stars) {
            float dx = other.px - star.px;
            float dy = other.py - star.py;
            float dz = other.pz - star.pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            energy -= other.mass * star.mass * G / std::sqrt(d2) * 0.5;
        }
    }
    return energy;
}

template <class Func>
long benchmark(Func const &func) {
    auto t0 = std::chrono::steady_clock::now();
    func();
    auto t1 = std::chrono::steady_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    return dt.count();
}

int main() {
    init();
    printf("Initial energy: %f\n", calc());
    auto dt = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
