#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

// 减少初始化时的乘法
constexpr float RAND_MAX2 = 1.f / RAND_MAX * 2;
float frand() {
    return (float)rand() * RAND_MAX2 - 1;
}

constexpr size_t length = 48;
struct Star {
    float px[length], py[length], pz[length];
    float vx[length], vy[length], vz[length];
    float mass[length];
};

Star stars;

void init() {
    #pragma GCC unroll 8
    for (size_t i = 0; i < length; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1;
    }
}

constexpr float G = 0.001;
constexpr float eps = 0.001;
constexpr float dt = 0.01;
// 提前计算, 减少不必要的乘法
constexpr float Gdt = G * dt;
constexpr float eps2 = eps * eps;
constexpr float G2 = G / 2;

void step() {
    for (size_t i = 0; i < length; i++) {
        // 减少不必要的内存访问
        float spxi = stars.px[i];
        float spyi = stars.py[i];
        float spzi = stars.pz[i];

        // 先累加到初始化为0的局部变量
        float tmp_vxi = 0;
        float tmp_vyi = 0;
        float tmp_vzi = 0;
        
        #pragma GCC unroll 8
        for (size_t j = 0; j < length; j++) {
            float dx = stars.px[j] - spxi;
            float dy = stars.py[j] - spyi;
            float dz = stars.pz[j] - spzi;
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            // Gdt = G * dt 放到for循环外部
            d2 = 1.f / d2;
            // 乘法变加法
            tmp_vxi += dx * stars.mass[j] * d2;
            tmp_vyi += dy * stars.mass[j] * d2;
            tmp_vzi += dz * stars.mass[j] * d2;
        }
        // 累加结束后再写入到全局变量中
        stars.vx[i] += tmp_vxi * Gdt;
        stars.vy[i] += tmp_vyi * Gdt;
        stars.vz[i] += tmp_vzi * Gdt;
    }

    #pragma GCC unroll 8
    for (size_t i = 0; i < length; i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for (size_t i = 0; i < length; i++) {
        // 减少不必要的内存访问
        float pxi = stars.px[i];
        float pyi = stars.py[i];
        float pzi = stars.pz[i];
        float massi = stars.mass[i];
        
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += massi * v2 / 2;

        // 先累加到初始化为0的局部变量
        float tmp = 0;
        #pragma GCC unroll 8
        for (size_t j = 0; j < length; j++) {
            float dx = stars.px[j] - pxi;
            float dy = stars.py[j] - pyi;
            float dz = stars.pz[j] - pzi;
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            // 将massi = stars.mass[i]和G2 = G / 2放到for循环外部
            // 减少乘法次数
            tmp += stars.mass[j] / std::sqrt(d2);
        }
        // 累加结束后写入
        energy -= tmp * massi * G2;
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
