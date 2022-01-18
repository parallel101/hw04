#include <cstdio>
#include <cstdlib>
#include <array>
#include <chrono>
#include <cmath>

constexpr size_t g_len = 48;

float G = 0.001;
float eps = 0.001;
float dt = 0.01;

struct StarVec {
    std::array<float, g_len> px, py, pz;
    std::array<float, g_len> vx, vy, vz;
    std::array<float, g_len> mass;
};

StarVec stars;

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

void init() {
    printf("SOA mode\n");
    for (int i = 0; i < g_len; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1;
    }
}

void step() {
    float eps2 = eps * eps;
    for (size_t i=0;i<g_len;i++) {
        float px = stars.px[i], py = stars.py[i], pz = stars.pz[i];
        float vx = .0f, vy = .0f, vz = .0f;

        for (size_t j=0;j<g_len;j++) {
            float dx = stars.px[j] - px;
            float dy = stars.py[j] - py;
            float dz = stars.pz[j] - pz;

            float d2 = 1.0 / (dx * dx + dy * dy + dz * dz + eps2);
            d2 *= std::sqrt(d2);
            float dt_G_mass = dt * G * stars.mass[j];
            
            vx += dx * d2 * dt_G_mass;
            vy += dy * d2 * dt_G_mass;
            vz += dz * d2 * dt_G_mass;
        }

        stars.vx[i] += vx;
        stars.vy[i] += vy;
        stars.vz[i] += vz;
    }

    for (size_t i=0;i<g_len;i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for (size_t i=0;i<g_len;i++) {
        float px = stars.px[i], py = stars.py[i], pz = stars.pz[i];
        float vx = stars.vx[i], vy = stars.vy[i], vz = stars.vz[i];
        float v2 = vx * vx + vy * vy + vz * vz;
        energy += stars.mass[i] * v2 / 2;
        for (size_t j=0;j<g_len;j++) {
            float dx = stars.px[j] - px;
            float dy = stars.py[j] - py;
            float dz = stars.pz[j] - pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars.mass[j] * stars.mass[i] * G / std::sqrt(d2) / 2;
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