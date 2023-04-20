#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

constexpr float daoshu = 1 / RAND_MAX;
float frand() {
    return (float)rand() * daoshu * 2 - 1;
}


constexpr int N = 48;
struct Star {
    float px[N];
    float py[N];
    float pz[N];
    float vx[N];
    float vy[N];
    float vz[N];
    float mass[N];
};

Star stars;

void init() {
    for (size_t i = 0; i < 48; i++) {
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
    
        for(size_t i = 0; i < 48; i++) {
            for(size_t j = 0; j < 48; j++) {
                float dx = stars.px[j] - stars.px[i];
                float dy = stars.py[j] - stars.py[i];
                float dz = stars.pz[j] - stars.pz[i];
                float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
                d2 *= std::sqrt(d2);
                float dao = 1 / d2;
                
                stars.vx[i] += dx * stars.mass[j] * (G * dt * dao);
                stars.vy[i] += dy * stars.mass[j] * (G * dt * dao);
                stars.vz[i] += dz * stars.mass[j] * (G * dt * dao);
            }
            stars.px[i] += stars.vx[i] * dt;
            stars.py[i] += stars.vy[i] * dt;
            stars.pz[i] += stars.vz[i] * dt;
        }
    }

float calc() {
    float energy = 0;
    for(size_t i = 0; i < 48; ++i) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 * 0.5f;
        for(size_t j = 0; j < 48; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + (eps * eps);
            energy -= stars.mass[j] * stars.mass[i]  / std::sqrt(d2) * (0.5f * G);
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
