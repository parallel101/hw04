#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

static float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct Star {
    float px[48], py[48], pz[48];
    float vx[48], vy[48], vz[48];
    float mass[48];
};

Star stars;

constexpr size_t pnum = 48;

void init() {
    for (size_t i = 0; i < pnum; i++) {
        stars.px[i]=frand();
        stars.py[i]=frand();
        stars.pz[i]=frand();
        stars.vx[i]=frand();
        stars.vy[i]=frand();
        stars.vz[i]=frand();
        stars.mass[i]=frand()+1;
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;
float eps2 = eps*eps;
float Gpdt = G * dt;

void step() {
#pragma GCC ivdep
// #pragma omp simd
    for (size_t i=0; i<pnum; i++) {
        for (size_t j=0; j<pnum; j++) {
            float dx = stars.px[j]-stars.px[i];
            float dy = stars.py[j]-stars.py[i];
            float dz = stars.pz[j]-stars.pz[i];
            float d2 = dx*dx + dy*dy + dz*dz + eps2;
            d2 *= std::sqrt(d2);
            stars.vx[i] += dx * stars.mass[j] * Gpdt / d2;
            stars.vy[i] += dy * stars.mass[j] * Gpdt / d2;
            stars.vz[i] += dz * stars.mass[j] * Gpdt / d2;
        }
    }
    for (size_t i = 0; i<pnum; i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
    
}

float calc() {
    float energy = 0;
#pragma GCC ivdep
// #pragma omp simd
    for (size_t i = 0; i<pnum; i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        for (size_t j = 0; j<pnum; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx*dx + dy*dy + dz*dz + eps2;
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