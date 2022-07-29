#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

#define N 48

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct alignas(16) Star {
    std::array<float,N> px, py, pz;
    std::array<float,N> vx, vy, vz;
    std::array<float,N> mass;
};

Star stars;

// 原计划的初始化，但是完全没有矢量化
// void init() {
//     for (size_t i = 0; i < N; ++i) {
//         stars.px[i] = frand();
//         stars.py[i] = frand();
//         stars.pz[i] = frand();
//         stars.vx[i] = frand();
//         stars.vy[i] = frand();
//         stars.vz[i] = frand();
//         stars.mass[i] = frand() + 1;
//     }
// }

template <int M = 0>
void sub_init(std::array<float,N> &arr) {
    #pragma GCC unroll 4
    for (size_t i = 0; i < N; ++i) {
        arr[i] = frand() + M;
    }
}

void init() {
    sub_init<0>(stars.px);
    sub_init<0>(stars.py);
    sub_init<0>(stars.pz);
    sub_init<0>(stars.vx);
    sub_init<0>(stars.vy);
    sub_init<0>(stars.vz);
    sub_init<1>(stars.mass);
}


float G = 0.001;
float eps = 0.001;
float dt = 0.01;

float eps2 = eps * eps;
float Gdt = G * dt;

void step() {
    for (size_t i = 0; i < N; ++i) {
        float dvx = 0, dvy = 0, dvz = 0;
        #pragma GCC unroll 4
        for (size_t j = 0; j < N; ++j) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);

            float tmp = stars.mass[j] * Gdt / d2;
            dvx += dx * tmp;
            dvy += dy * tmp;
            dvz += dz * tmp;
        }
        stars.vx[i] += dvx;
        stars.vy[i] += dvy;
        stars.vz[i] += dvz;
    }
    #pragma GCC unroll 4
    for (size_t i = 0; i < N; ++i) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    #pragma GCC unroll 4
    for (size_t i = 0; i < N; ++i) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
    }
    for (size_t i = 0; i < N; ++i) {
        float tmp = stars.mass[i] * G / 2;
        #pragma GCC unroll 4
        for (size_t j = 0; j < N; ++j) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            float reverse_sqrt_d2 = 1 / std::sqrt(d2);
            energy -= stars.mass[j] * tmp *reverse_sqrt_d2;
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
