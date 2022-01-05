#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

float frand() {
    return (float)std::rand() / RAND_MAX * 2 - 1;
}

struct Star {
    alignas(16) std::array<float, 48> px, py, pz;
    alignas(16) std::array<float, 48> vx, vy, vz;
    alignas(16) std::array<float, 48> mass;
};

Star stars;

void init() {
    for (size_t i = 0; i < 48; i++) {
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
float eps2 = eps * eps;
float dt = 0.01;
float G_dt = G * dt;

void step() {
    for (size_t i = 0; i < 48; i++) {
        float px = stars.px[i];
        float py = stars.py[i];
        float pz = stars.pz[i];
        float vx_temp = 0.0f;
        float vy_temp = 0.0f;
        float vz_temp = 0.0f;
        for (size_t j = 0; j < 48; j++) {
            float dx = stars.px[j] - px;
            float dy = stars.py[j] - py;
            float dz = stars.pz[j] - pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            float mass_g_dt_inv_d2 = stars.mass[j] * G_dt / d2;
            vx_temp += dx * mass_g_dt_inv_d2;
            vy_temp += dy * mass_g_dt_inv_d2;
            vz_temp += dz * mass_g_dt_inv_d2;
        }
        stars.vx[i] += vx_temp;
        stars.vy[i] += vy_temp;
        stars.vz[i] += vz_temp;
    }

    for (size_t i = 0; i < 48; i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for (size_t i = 0; i < 48; i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 * 0.5f;
        float px = stars.px[i];
        float py = stars.py[i];
        float pz = stars.pz[i];
        float mass_G = stars.mass[i] * G;
        for (size_t j = 0; j < 48; j++) {
            float dx = stars.px[j] - px;
            float dy = stars.py[j] - py;
            float dz = stars.pz[j] - pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars.mass[j] * mass_G / std::sqrt(d2) * 0.5f;
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
