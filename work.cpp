#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>
#include <immintrin.h>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct Star {
    float px, py, pz;
    float vx, vy, vz;
    float mass;
};

constexpr int NSTAR = 48;
struct StarSOA
{
    float data[7][NSTAR];
    inline float* operator[] (size_t pos) { return data[pos];}
    inline const float* operator[] (const size_t pos) const { return &(*data[pos]);}
};
StarSOA stars2;
std::vector<Star> stars;

float G = 0.001;
float eps = 0.001;
float dt = 0.01;



void init() {
    for (int i = 0; i < 48; i++) {
        stars.push_back({
            frand(), frand(), frand(),
            frand(), frand(), frand(),
            frand() + 1,
        });
    }
}

void step() {
    for (auto &star: stars) {
        for (auto &other: stars) {
            float dx = other.px - star.px;
            float dy = other.py - star.py;
            float dz = other.pz - star.pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            d2 *= sqrt(d2);
            star.vx += dx * other.mass * G * dt / d2;
            star.vy += dy * other.mass * G * dt / d2;
            star.vz += dz * other.mass * G * dt / d2;
        }
    }
    for (auto &star: stars) {
        star.px += star.vx * dt;
        star.py += star.vy * dt;
        star.pz += star.vz * dt;
    }
}


float calc() {
    float energy = 0;
    for (auto &star: stars) {
        float v2 = star.vx * star.vx + star.vy * star.vy + star.vz * star.vz;
        energy += star.mass * v2 / 2;
        for (auto &other: stars) {
            float dx = other.px - star.px;
            float dy = other.py - star.py;
            float dz = other.pz - star.pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= other.mass * star.mass * G / sqrt(d2) / 2;
        }
    }
    return energy;
}
void init2() {
    // for (int i=0; i<NSTAR; i++) {
    //     for (int j=5; j>=0; j--) {
    //         stars2[j][i] = frand();
    //     }
    //     stars2[6][i] = frand() + 1;
    // }
    for (int i=0; i<NSTAR; i++) {
        stars2[0][i] = stars[i].px;
        stars2[1][i] = stars[i].py;
        stars2[2][i] = stars[i].pz;
        stars2[3][i] = stars[i].vx;
        stars2[4][i] = stars[i].vy;
        stars2[5][i] = stars[i].vz;
        stars2[6][i] = stars[i].mass;
    }
}

// 借鉴雷神之锤3的速算法
inline float Q_rsqrt(float number)
{
    constexpr float threehalfs = 1.5F;
    float x2 = number * 0.5F;
    float y = number;
    int i = *reinterpret_cast<int*>(&y);           // evil floating point bit level hacking
    i = 0x5f3759df - (i >> 1); 
    y = *reinterpret_cast<float*>(&i);
    y = y * (threehalfs - (x2 * y * y)); // 1st iteration
    y = y * (threehalfs - (x2 * y * y)); // 2nd iteration, this can be removed

    return y;
}

void step2() {
    auto other = stars2;
    for (int i=0; i<NSTAR; i++) {
        float vx{0}, vy{0}, vz{0};
        for (int j=0; j<NSTAR; j++) {
            float dx = other[0][j] - stars2[0][i];
            float dy = other[1][j] - stars2[1][i];
            float dz = other[2][j] - stars2[2][i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            // d2 = 1./(std::sqrt(d2)*d2);     // 292 ms
            d2 = Q_rsqrt(d2)/d2;            // 273ms
            vx += dx * other[6][j] * (G * dt) * d2;
            vy += dy * other[6][j] * (G * dt) * d2;
            vz += dz * other[6][j] * (G * dt) * d2;
        }
        stars2[3][i] += vx;
        stars2[4][i] += vy;
        stars2[5][i] += vz;
        stars2[0][i] += stars2[3][i] * dt;
        stars2[1][i] += stars2[4][i] * dt;
        stars2[2][i] += stars2[5][i] * dt;
    }
}


float calc2() {
    float energy = 0;
    for (int i=0; i<NSTAR; i++) {
        float v2 = stars2[3][i] * stars2[3][i] + stars2[4][i] * stars2[4][i] + stars2[5][i] * stars2[5][i];
        energy += stars2[6][i] * v2 /2;
        for (int j=0; j<NSTAR; j++) {
            float dx = stars2[0][j] - stars2[0][i];
            float dy = stars2[1][j] - stars2[1][i];
            float dz = stars2[2][j] - stars2[2][i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars2[6][j] * stars2[6][j] * G / sqrt(d2) / 2;            
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
    init2();
    printf("original:\n");
    printf("Initial energy: %f\n", calc());
    auto dt = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);

    printf("my:\n");
    printf("Initial energy: %f\n", calc2());
    auto dt2 = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step2();
    });
    printf("Final energy: %f\n", calc2());
    printf("Time elapsed: %ld ms\n", dt2);    
    return 0;
}