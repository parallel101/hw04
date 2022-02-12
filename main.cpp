#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

constexpr int size = 48;

struct Stars {
    float px[size], py[size], pz[size];
    float vx[size], vy[size], vz[size];
    float mass[size];
};
Stars stars;

void init() {
    for (int i = 0; i < size; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand()+1;
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;
float Gdt = G*dt;
float epssqure = eps*eps;

void step() {
#pragma omp simd
    for (size_t i = 0; i != size; i++) {
        for (size_t j = 0; j != size; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + epssqure;
            d2 *= std::sqrt(d2);
            float x = stars.mass[j] * Gdt / d2;
            stars.vx[i] += dx * x;
            stars.vy[i] += dy * x;
            stars.vz[i] += dz * x;
        }
    }
    for (size_t i = 0; i != size; i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for (size_t i = 0; i != size; i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        for (size_t j = 0; j != size; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
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
