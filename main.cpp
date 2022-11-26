#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

const int MAXLEN = 48;

struct Star {
    float px[MAXLEN], py[MAXLEN], pz[MAXLEN];
    float vx[MAXLEN], vy[MAXLEN], vz[MAXLEN];
    float mass[MAXLEN];
} stars;

void init() {
    for (int i = 0; i < MAXLEN; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1;
    }
}

const float G = 0.001;
const float eps = 0.001;
const float dt = 0.01;

const float Gdt = G * dt;

void step() {
    for (size_t i = 0; i < MAXLEN; ++i)
    {
        for (size_t j = 0; j < MAXLEN; ++j)
        {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            d2 = 1 / std::sqrt(d2);
            d2 = d2 * d2 * d2;
            stars.vx[i] += dx * stars.mass[j] * Gdt * d2;
            stars.vy[i] += dy * stars.mass[j] * Gdt * d2;
            stars.vz[i] += dz * stars.mass[j] * Gdt * d2;
        }
    }
    
    // for (auto &star: stars) {
    //     for (auto &other: stars) {
    //         float dx = other.px - star.px;
    //         float dy = other.py - star.py;
    //         float dz = other.pz - star.pz;
    //         float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
    //         d2 = 1 / std::sqrt(d2);
    //         d2 = d2 * d2 * d2;
            
    //     }
    // }

    for(size_t i = 0; i < MAXLEN; ++i) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for (size_t i = 0; i < MAXLEN; ++i)
    {
        float starvx = stars.vx[i], starvy = stars.vy[i], starvz = stars.vz[i];
        float starpx = stars.px[i], starpy = stars.py[i], starpz = stars.pz[i];
        float starmass = stars.mass[i];
        float v2 = starvx * starvx + starvy * starvy + starvz * starvz;
        energy += starmass * v2 / 2;
        for (size_t j = 0; j < MAXLEN; ++j) {
            float dx = stars.px[j] - starpx;
            float dy = stars.py[j] - starpy;
            float dz = stars.pz[j] - starpz;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars.mass[j] * starmass * G / sqrt(d2) / 2;
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
