#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

inline float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}
// AoS:Array of struct
//struct Star {
//    float px, py, pz;
//    float vx, vy, vz;
//    float mass;
//};
// SoA:struct of array
const int N = 48;
template<int N=48>
struct Star {
    std::array<float,N> px,py,pz;
    std::array<float,N> vx,vy,vz,mass;
};

//std::vector<Star> stars;
Star<N> stars;

void init() {
//    for (size_t i = 0; i < 48; i++) {
//        stars.push_back({
//            frand(), frand(), frand(),
//            frand(), frand(), frand(),
//            frand() + 1,
//        });
//    }
#pragma GCC unroll 4
    for(size_t i = 0;i < 48;i++){
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

void step() {
//    for (auto &star: stars) {
//        for (auto &other: stars) {
//            float dx = other.px - star.px;
//            float dy = other.py - star.py;
//            float dz = other.pz - star.pz;
//            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
//            d2 *= sqrt(d2);
//            float inve_d2 = 1.0 / d2;
////            star.vx += dx * other.mass * G * dt / d2;
////            star.vy += dy * other.mass * G * dt / d2;
////            star.vz += dz * other.mass * G * dt / d2;
//            star.vx += dx * other.mass * G * dt * inve_d2;
//            star.vy += dy * other.mass * G * dt * inve_d2;
//            star.vz += dz * other.mass * G * dt * inve_d2;
//        }
//    }
float gt =  G * dt;
#pragma GCC unroll 4
    for(size_t i = 0;i < N;i++){
        for(size_t j = 0;j < N;j++){
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            d2 *= std::sqrt(d2);
            float inve_d2 = gt / d2;
//            star.vx += dx * other.mass * G * dt / d2;
//            star.vy += dy * other.mass * G * dt / d2;
//            star.vz += dz * other.mass * G * dt / d2;
            stars.vx[i] += dx * stars.mass[j] * inve_d2;
            stars.vy[i] += dy * stars.mass[j] * inve_d2;
            stars.vz[i] += dz * stars.mass[j] * inve_d2;
        }
    }


//    for (auto &star: stars) {
//        star.px += star.vx * dt;
//        star.py += star.vy * dt;
//        star.pz += star.vz * dt;
//    }
#pragma GCC unroll 4
    for(size_t k = 0;k < N;k++){
        stars.px[k] += stars.vx[k] * dt;
        stars.py[k] += stars.vy[k] * dt;
        stars.pz[k] += stars.vz[k] * dt;
    }


}

float calc() {
    float energy = 0;
//    for (auto &star: stars) {
//        float v2 = star.vx * star.vx + star.vy * star.vy + star.vz * star.vz;
////        energy += star.mass * v2 / 2;
//        energy += star.mass * v2 * 0.5;
//        for (auto &other: stars) {
//            float dx = other.px - star.px;
//            float dy = other.py - star.py;
//            float dz = other.pz - star.pz;
//            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
////            energy -= other.mass * star.mass * G / sqrt(d2) / 2;
//            energy -= other.mass * star.mass * G / sqrt(d2) * 0.5;
//        }
//    }
//    return energy;
const float eps2 = eps * eps;
#pragma GCC unroll 4
    for(size_t i = 0;i < N;i++){
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 * 0.5;
        for(size_t j = 0;j < N;j++) {
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
        for (size_t i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
