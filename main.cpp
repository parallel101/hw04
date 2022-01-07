#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

constexpr float sDoubleInvRandMax = 2.0f / RAND_MAX;

float frand() {
    return (float)rand() * sDoubleInvRandMax - 1;
}

template <size_t N>
struct Star {
    float px[N];
    float py[N];
    float pz[N];
    float vx[N];
    float vy[N];
    float vz[N];
    float mass[N];
};

constexpr size_t sStarCount = 48;
Star<sStarCount> stars;

void init() {
#pragma GCC unroll 16
    for (size_t i = 0; i < sStarCount; ++i)
    {
        stars.px[i] = frand(); 
        stars.py[i] = frand(); 
        stars.pz[i] = frand(); 
        stars.vx[i] = frand(); 
        stars.vy[i] = frand(); 
        stars.vz[i] = frand(); 
        stars.mass[i] = frand() + 1;
    }
}

constexpr float G = 0.001;
constexpr float eps = 0.001;
constexpr float eps2 = eps * eps;
constexpr float dt = 0.01;

void step() {
#pragma omp simd
    for(size_t i = 0; i < sStarCount; ++i)
    {
        float sumx = 0.0f;
        float sumy = 0.0f;
        float sumz = 0.0f;

        for(size_t j = 0; j < sStarCount; ++j)
        {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float G_dt_inv_d2_mass = (G * dt) * std::pow(dx * dx + dy * dy + dz * dz + eps2, -1.5f) * stars.mass[j];

            sumx += G_dt_inv_d2_mass * dx;
            sumy += G_dt_inv_d2_mass * dy;
            sumz += G_dt_inv_d2_mass * dz;
        }

        stars.vx[i] += sumx;
        stars.vy[i] += sumy;
        stars.vz[i] += sumz;
    }

    for(size_t i = 0; i < sStarCount; ++i)
    {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for(size_t i = 0; i < sStarCount; ++i)
    {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 * 0.5f;

        float G_mass_i_half = (G * 0.5f) * stars.mass[i];

        for(size_t j = 0; j < sStarCount; ++j)
        {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            energy -= stars.mass[j] * G_mass_i_half / std::sqrt(d2);
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
