#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

constexpr size_t star_size = 48;
constexpr float G = 0.001;
constexpr float eps = 0.001;
constexpr float dt = 0.01;
constexpr float eps2 = eps * eps;
constexpr float Gdt = G * dt;

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

template<size_t N>
struct Star {
    float px[N], py[N], pz[N];
    float vx[N], vy[N], vz[N];
    float mass[N];
};

Star<star_size> star;

void init() {
    for (int i = 0; i < star_size; i++) {
        star.px[i] = frand();
        star.py[i] = frand();
        star.pz[i] = frand();
        star.vx[i] = frand();
        star.vy[i] = frand();
        star.vz[i] = frand();
        star.mass[i] = frand() + 1;
    }
}

void step() {
    for (size_t i = 0; i < star_size; ++i) {
        float sumx = 0.0f, sumy = 0.0f, sumz = 0.0f;
        for (size_t j = 0; j < star_size; ++j) {
            float dx = star.px[j] - star.px[i];
            float dy = star.py[j] - star.py[i];
            float dz = star.pz[j] - star.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            sumx += dx * star.mass[j] * Gdt / d2;
            sumy += dy * star.mass[j] * Gdt / d2;
            sumz += dz * star.mass[j] * Gdt / d2;
        }
        star.vx[i] += sumx;
        star.vy[i] += sumy;
        star.vz[i] += sumz;
    }
    for (size_t i = 0; i < star_size; ++i) {
        star.px[i] += star.vx[i] * dt;
        star.py[i] += star.vy[i] * dt;
        star.pz[i] += star.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for (size_t i = 0; i < star_size; ++i) {
        float v2 = star.vx[i] * star.vx[i] + star.vy[i] * star.vy[i] + star.vz[i] * star.vz[i];
        energy += star.mass[i] * v2 / 2;
        for (size_t j = 0; j < star_size; ++j) {
            float dx = star.px[j] - star.px[i];
            float dy = star.py[j] - star.py[i];
            float dz = star.pz[j] - star.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            energy -= star.mass[j] * star.mass[i] * G / std::sqrt(d2) / 2;
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
