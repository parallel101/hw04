#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

static float frand() {
    return (float)std::rand() / RAND_MAX * 2 - 1;
}

const size_t N = 48;
const float G = 0.001f;
const float eps2 = 0.001f * 0.001f;
const float dt = 0.01f;

struct Star {
    std::array<float, N> px, py, pz;
    std::array<float, N> vx, vy, vz;
    std::array<float, N> mass;
};

Star stars;

void init() {
    for (size_t i = 0; i < N; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1.0f;
    }
}

void step() {
    float Gdt = G * dt;
    for (size_t i = 0; i < N; i++)
    {
        float px = stars.px[i];
        float py = stars.py[i];
        float pz = stars.pz[i];
        float vx = 0.0f, vy = 0.0f, vz = 0.0f;
        for (size_t j = 0; j < N; j++)
        {
            float dx = stars.px[j] - px;
            float dy = stars.py[j] - py;
            float dz = stars.pz[j] - pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            float acc = stars.mass[j] * Gdt / d2;
            vx += dx * acc;
            vy += dy * acc;
            vz += dz * acc;
        }
        stars.vx[i] += vx;
        stars.vy[i] += vy;
        stars.vz[i] += vz;
    }
    for (size_t i = 0; i < N; i++)
    {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0.0f;
    for (size_t i = 0; i < N; i++)
    {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2.0f;
        float minus = 0.0f;
        for (size_t j = 0; j < N; j++)
        {
            if(i == j) continue;
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            minus += stars.mass[j] * stars.mass[i] * G / std::sqrt(d2) / 2.0f;
        }
        energy -= minus;
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
