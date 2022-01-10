#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>
float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}
const size_t num = 48;
struct Star {
    float px[num], py[num], pz[num];
    float vx[num], vy[num], vz[num];
    float mass[num];
};

Star stars;

void init() {
    for (int i = 0; i < num; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand()+1;
        // stars[i] = {
        //     frand(), frand(), frand(),
        //     frand(), frand(), frand(),
        //     frand() + 1,
        // };
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;
float eps2 = eps * eps;
void step() {
    float G_dt = G*dt;
    // for (auto &star: stars) {
    //     for (auto &other: stars) {
    //         float dx = other.px - star.px;
    //         float dy = other.py - star.py;
    //         float dz = other.pz - star.pz;
    //         float d2 = dx * dx + dy * dy + dz * dz + eps2;
    //         d2 *= std::sqrt(d2);
    //         float inv_d2 = 1.0f/d2;
    //         star.vx += dx * other.mass * G_dt * inv_d2;
    //         star.vy += dy * other.mass * G_dt * inv_d2;
    //         star.vz += dz * other.mass * G_dt * inv_d2;
    //     }
    // }
    for (size_t i = 0;i<num;i++) {
        for (size_t j =0;j<num;j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            float inv_d2 = 1.0f/d2;
            stars.vx[i] += dx * stars.mass[j] * G_dt * inv_d2;
            stars.vy[i] += dy * stars.mass[j] * G_dt * inv_d2;
            stars.vz[i] += dz * stars.mass[j] * G_dt * inv_d2;
        }
    }
    // for (auto &star: stars) {
    //     star.px += star.vx * dt;
    //     star.py += star.vy * dt;
    //     star.pz += star.vz * dt;
    // }
    for (size_t i = 0;i<num;i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    // for (auto &star: stars) {
    //     float v2 = star.vx * star.vx + star.vy * star.vy + star.vz * star.vz;
    //     energy += star.mass * v2 * 0.5;
    //     for (auto &other: stars) {
    //         float dx = other.px - star.px;
    //         float dy = other.py - star.py;
    //         float dz = other.pz - star.pz;
    //         float d2 = dx * dx + dy * dy + dz * dz + eps;
    //         energy -= other.mass * star.mass * G / std::sqrt(d2) * 0.5;
    //     }
    // }
    for (size_t i = 0;i<num;i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 * 0.5;
        for (size_t j =0;j<num;j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            energy -= stars.mass[j] * stars.mass[i] * G / std::sqrt(d2) * 0.5;
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
