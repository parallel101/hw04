#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

static float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct alignas(64) Star {
    float px[48], py[48], pz[48];
    float vx[48], vy[48], vz[48];
    float mass[48];
};

Star stars;

void init() {
    #pragma GCC unroll 4
    for (int i = 0; i < 48; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1;
    }
}

float G = 0.001;
float eps = 0.001;
float eps_square = eps * eps;
float dt = 0.01;

void step() {
    #pragma GCC unroll 4
    for (int i = 0; i < 48; i++) {
        float px = stars.px[i], py = stars.py[i], pz = stars.pz[i];
        float vx_plus = 0, vy_plus = 0, vz_plus = 0;
        #pragma GCC unroll 4
        for (int j = 0; j < 48; j++) {
            float dx = stars.px[j] - px;
            float dy = stars.py[j] - py;
            float dz = stars.pz[j] - pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps_square;
            d2 *= std::sqrt(d2);
            float value = stars.mass[j] * (G * dt) / d2;
            vx_plus += dx * value;
            vy_plus += dy * value;
            vz_plus += dz * value;
        }
        stars.vx[i] += vx_plus;
        stars.vy[i] += vy_plus;
        stars.vz[i] += vz_plus;
    }
    #pragma GCC unroll 4
    for (int i = 0; i < 48; i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for (int i = 0; i < 48; i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        for (int j = 0; j < 48; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars.mass[j] * stars.mass[i] * G / sqrt(d2) / 2;
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
        #pragma GCC unroll 64
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
