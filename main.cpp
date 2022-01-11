#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

constexpr float speedup = 1.0 / RAND_MAX;

float frand() {
    return (float)rand() / speedup * 2 - 1;
}
constexpr int length = 48;

struct Star {
    float px[length]; 
    float py[length]; 
    float pz[length];
    float vx[length]; 
    float vy[length]; 
    float vz[length];
    float mass[length];
};

Star stars;

void init() {
    for (int i = 0; i < 48; i++) {
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
    size_t len = length;
    float eps2 = eps * eps;
    float  gdt = G * dt;
    #pragma GCC unroll 16
    for (size_t i = 0 ; i < len; i++) {
            float dxs[length];
            float dys[length];
            float dzs[length];
            float d2s[length];
            float ivf_d2s[length];
            #pragma opm simd
            for(size_t j=0; j < len; j++)
            {
                dxs[j] = stars.px[j] - stars.px[i];
            }
            #pragma opm simd
            for(size_t j=0; j < len; j++)
            {
                dys[j] = stars.py[j] - stars.py[i];
            }
            #pragma opm simd
            for(size_t j=0; j < len; j++)
            {
                dzs[j] = stars.pz[j] - stars.pz[i];
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++)
            {
                d2s[j] = dxs[j] * dxs[j] + dys[j] * dys[j] + dzs[j] * dzs[j] + eps2;
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++){
                ivf_d2s[j] = 1.0 / (d2s[j] * std::sqrt(d2s[j]));
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++){
                stars.vx[i] += dxs[j] * stars.mass[j] * (gdt * ivf_d2s[j]);
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++){
                stars.vy[i] += dys[j] * stars.mass[j] * (gdt * ivf_d2s[j]);
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++){
                stars.vz[i] += dzs[j] * stars.mass[j] * (gdt * ivf_d2s[j]);
            }
        }
        #pragma opm simd
        for(size_t i=0; i<len; i++)
        {
            stars.px[i] += stars.vx[i] * dt ;
        }
        #pragma opm simd
        for(size_t i=0; i < len; i++)
        {
            stars.py[i] += stars.vy[i] * dt ;
        }
        #pragma opm simd
        for(size_t i = 0; i < len; i++)
        {
            stars.pz[i] += stars.vz[i] * dt;
        }
}

float calc() {
    float energy = 0;
    size_t len = length;
    for (size_t i = 0; i < len; i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i]* stars.vy[i]+ stars.vz[i]* stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        #pragma GCC unroll 32
        for (size_t j=0; j < len; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            float ivf_d2 = 1.0 / (std::sqrt(d2) * 2);
            energy -= stars.mass[j] * stars.mass[j] * (G * ivf_d2);
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