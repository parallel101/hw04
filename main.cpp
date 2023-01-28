#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

const size_t N = 48;

float frand() {
    return (float)std::rand() / RAND_MAX * 2 - 1;
}

struct Star {
    float px[N], py[N], pz[N];
    float vx[N], vy[N], vz[N];
    float mass[N];
};

Star stars;

void init() {
    #pragma GCC unroll 16
    for (size_t i = 0; i < N; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1.f;
    }
}

const float G = 0.001;
const float eps = 0.001;
const float dt = 0.01;
const float eps2 = eps * eps;
const float dt2 = dt * dt;
const float G_plus_dt = G * dt;

void step() {
    for (size_t i = 0; i < N; i++) {
        #pragma GCC unroll 16
        for (size_t j = 0; j < N; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            float inv_d2 = 1.f / d2;
            stars.vx[i] += dx * stars.mass[j] * G_plus_dt * inv_d2;
            stars.vy[i] += dy * stars.mass[j] * G_plus_dt * inv_d2;
            stars.vz[i] += dz * stars.mass[j] * G_plus_dt * inv_d2;
        }
    }
    #pragma GCC unroll 16
    for (size_t i = 0; i < N; i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    float half_G = G * 0.5f;
    for (size_t i = 0; i < N; i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 * 0.5f;
        #pragma GCC unroll 16
        for (size_t j = 0; j < N; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            float inv_sqrt_d2 = 1.f / std::sqrt(d2);
            energy -= stars.mass[j] * stars.mass[i] * half_G * inv_sqrt_d2;
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
