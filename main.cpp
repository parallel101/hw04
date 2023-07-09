#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

int randMax = 1 / RAND_MAX;

const size_t N = 48;

float frand() {
    return (float)rand() *randMax * 2 - 1;
}

struct Star {
    float px[N], py[N], pz[N];
    float vx[N], vy[N], vz[N];
    float mass[N];
};

Star stars;


void init() {
    float tmp = frand();
    #pragma GCC unroll 16
    for (int i = 0; i < 48; i++) {
        stars.px[i] = tmp;
        stars.py[i] = tmp;
        stars.pz[i] = tmp;
        stars.vx[i] = tmp;
        stars.vy[i] = tmp;
        stars.vz[i] = tmp;
        stars.mass[i] = tmp+1;
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;
float Gdt = G * dt;
float eps2 = eps * eps;

void step() {
    for(size_t i = 0; i < 48; i++) {
        #pragma GCC unroll 16
        for(size_t j = 0; j < 48; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            float dao = 1 / d2;

            stars.vx[i] += dx * stars.mass[j] * (Gdt * dao);
            stars.vy[i] += dy * stars.mass[j] * (Gdt * dao);
            stars.vz[i] += dz * stars.mass[j] * (Gdt * dao);
        }
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    Star other = stars;
    Star star = stars;  
    for (std::size_t i = 0;i <48;i++) {
        float v2 = star.vx[i] * star.vx[i] + star.vy[i] * star.vy[i] + star.vz[i] * star.vz[i];
        energy += star.mass[i] * v2 *0.5f;
        #pragma GCC unroll 16
        for (std::size_t j = 0;j < 48;j++) {
            float dx = other.px[j] - star.px[i];
            float dy = other.py[j] - star.py[i];
            float dz = other.pz[j] - star.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            float tmp = 1 / sqrt(d2);
            energy -= other.mass[j] * star.mass[i] * G * tmp *0.5f;
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
