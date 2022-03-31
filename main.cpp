#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>
float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct Star {
    float px, py, pz;
    float vx, vy, vz;
    float mass;
    
};

std::array<Star,48> stars;

void init() {
  //  #pragma omp simd
    #pragma GCC ivdep
    for (int i = 0; i < 48; i++) {
        stars[i] = (Star{
            frand(), frand(), frand(),
            frand(), frand(), frand(),
            frand() + 1,
        });
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;

void step() {
    float G = 0.001;
float eps = 0.001;
float dt = 0.01;
float epss = eps*eps;
float gpdt = G*dt;
//#pragma omp simd
#pragma GCC ivdep
    for (int i=0;i<48;i++) {
        for (int j=0;j<48;j++) {
            float dx = stars[j].px - stars[i].px;
            float dy = stars[j].py - stars[i].py;
            float dz = stars[j].pz - stars[i].pz;
            float d2 = dx * dx + dy * dy + dz * dz + epss;
            d2 *= sqrt(d2);
            stars[i].vx += dx * stars[j].mass * gpdt / d2;
            stars[i].vy += dy * stars[j].mass * gpdt / d2;
            stars[i].vz += dz * stars[j].mass * gpdt / d2;
        }
    }
    for (int i=0;i<48;i++) {
        stars[i].px += stars[i].vx * dt;
        stars[i].py += stars[i].vy * dt;
        stars[i].pz += stars[i].vz * dt;
    }
}

float calc() {
    float G = 0.001;
float eps = 0.001;
float dt = 0.01;
float epss = eps*eps;
    float energy = 0;
    //#pragma omp simd
    #pragma GCC ivdep
    for (int i=0;i<48;i++) {
        float v2 = stars[i].vx * stars[i].vx + stars[i].vy * stars[i].vy + stars[i].vz * stars[i].vz;
        energy += stars[i].mass * v2 / 2;
        for (int j=0;j<48;j++) {
            float dx = stars[j].px - stars[i].px;
            float dy = stars[j].py - stars[i].py;
            float dz = stars[j].pz - stars[i].pz;
            float d2 = dx * dx + dy * dy + dz * dz + epss;
            energy -= stars[j].mass * stars[i].mass * G / sqrt(d2) / 2;
        }
    }
    return energy;
}

template <class Func>
constexpr long benchmark(Func const &func) {
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
	#pragma omp simd
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
