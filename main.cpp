#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

const int N = 48;

struct StarField {
    std::array<float, N> px, py, pz;
    std::array<float, N> vx, vy, vz;
    std::array<float, N> mass;
};

StarField stars;

void init() {
    for (int i = 0; i < N; i++) {
        stars.px.fill(frand());
        stars.py.fill(frand());
        stars.pz.fill(frand());
        stars.vx.fill(frand());
        stars.vy.fill(frand());
        stars.vz.fill(frand());
        stars.mass.fill(frand() + 1);
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;

void step() {
    float Gdt = G * dt;
    float eps_square = eps * eps;

    // 修改为外层是变化的
    for (size_t i = 0; i < N; i++) {
        float px = stars.px[i], py = stars.py[i], pz = stars.pz[i];
        float vx = stars.vx[i], vy = stars.vy[i], vz = stars.vz[i];

        #pragma GCC unroll 2
        for (size_t j = 0; j < N; j++) {
            float dx = stars.px[j] - px;
            float dy = stars.py[j] - py;
            float dz = stars.pz[j] - pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps_square;
            d2 *= std::sqrt(d2);
            vx += dx * stars.mass[j] * Gdt / d2;
            vy += dy * stars.mass[j] * Gdt / d2;
            vz += dz * stars.mass[j] * Gdt / d2;
        }

        stars.vx[i] = vx;
        stars.vy[i] = vy;
        stars.vz[i] = vz;
    }

    #pragma GCC unroll 3
    for (size_t i = 0; i < N; i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;

    for (int i = 0; i < N; i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 * 0.5;
        float px = stars.px[i];
        float py = stars.py[i];
        float pz = stars.pz[i];
        float mass_G_half = stars.mass[i] * G * 0.5;
        for (int j = 0; j < N; j++) {
            float dx = stars.px[j] - px;
            float dy = stars.py[j] - py;
            float dz = stars.pz[j] - pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars.mass[j] * mass_G_half / std::sqrt(d2);
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
