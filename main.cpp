#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

float frand() {
    return (float)std::rand() / RAND_MAX * 2 - 1;
}

const int N = 48;

struct Star {
    std::array<float, N> px, py, pz;
    std::array<float, N> vx, vy, vz;
    std::array<float, N> mass;
};

Star stars;

void init() {
    for (std::size_t i = 0; i < N; i++) {
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
float eps2 = eps * eps;
float dt = 0.01;

void step() {
    for (std::size_t i = 0; i < stars.px.size(); i++) {

        float vx_tmp = 0;
        float vy_tmp = 0;
        float vz_tmp = 0;
        for (std::size_t j = 0; j < stars.px.size(); j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            float d = 1 / d2;
            vx_tmp += dx * stars.mass[j] * G * dt * d;
            vy_tmp += dy * stars.mass[j] * G * dt * d;
            vz_tmp += dz * stars.mass[j] * G * dt * d;
        }
        stars.vx[i] += vx_tmp;
        stars.vy[i] += vy_tmp;
        stars.vz[i] += vz_tmp;
    }


    for (std::size_t i = 0; i < stars.px.size(); i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }

}

float calc() {
    float energy = 0;
    for (std::size_t i = 0; i < stars.px.size(); i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 * 0.5;
        for (std::size_t j = 0; j < stars.px.size(); j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            energy -= stars.mass[i] * stars.mass[j] * G / std::sqrt(d2) * 0.5;
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
