#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>
float frand() {
    return (float)std::rand() /RAND_MAX * 2- 1;
}

struct Star {
    std::array<float, 48> px, py, pz;
    std::array<float, 48>  vx, vy, vz;
    std::array<float, 48>  mass;
};

Star stars;

void init() {
    for (std::size_t i = 0; i < 48; i++) {
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
float dt = 0.01;

void step() {
    for (std::size_t i = 0; i < stars.px.size(); i++) {
        float temp_vx = 0;
        float temp_vy = 0;
        float temp_vz = 0;

        for (std::size_t j = 0; j < stars.px.size(); j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            d2 *= std::sqrt(d2);
            temp_vx += dx * stars.mass[j] * G * dt / d2;
            temp_vy += dy * stars.mass[j] * G * dt / d2;
            temp_vz += dz * stars.mass[j] * G * dt / d2;
        }
        stars.vx[i] += temp_vx;
        stars.vy[i] += temp_vy;
        stars.vz[i] += temp_vz;
    }

    float temp_px = 0;
    float temp_py = 0;
    float temp_pz = 0;
    for (std::size_t i = 0; i < stars.px.size(); i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }


}

float calc() {
    float energy = 0;
    for (std::size_t i = 0; i < stars.px.size(); i++)  {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        for  (std::size_t j = 0; j < stars.px.size(); j++) {
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
