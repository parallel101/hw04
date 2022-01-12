/*

Initial energy: -13.787754
Final energy: -13.404387
Time elapsed: 215 ms

*/

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

constexpr float g_ivf = 1.0 / RAND_MAX;
float frand() {
    return (float)rand() * g_ivf * 2 - 1;
}
constexpr int g_len = 48;

struct StarVec{
    float pxs[g_len];
    float pys[g_len];
    float pzs[g_len];
    float vxs[g_len];
    float vys[g_len];
    float vzs[g_len];
    float masses[g_len];
};

StarVec stars;

void init() {
    for (int i = 0; i < 48; i++) {
        stars.pxs[i] = frand();
        stars.pys[i] = frand();
        stars.pzs[i] = frand();
        stars.vxs[i] = frand();
        stars.vys[i] = frand();
        stars.vzs[i] = frand();
        stars.masses[i] = frand()+1;
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;


void step() {
    size_t len = g_len;
    float eps2 = eps * eps;
    float gdt = G * dt;
    #pragma GCC unroll 16
    for(size_t i=0; i<len; i++){
        float dxs[g_len];
        float dys[g_len];
        float dzs[g_len];
        float d2s[g_len];
        float ivf_d2s[g_len];
        #pragma opm simd        
        for(size_t j=0; j<len; j++){
            dxs[j] = stars.pxs[j] - stars.pxs[i];
        }
        #pragma opm simd
        for(size_t j=0; j<len; j++){
            dys[j] = stars.pys[j] - stars.pys[i];
        }
        #pragma opm simd
        for(size_t j=0; j<len; j++){
            dzs[j] = stars.pzs[j] - stars.pzs[i];
        }
        #pragma opm simd
        for(size_t j=0; j<len; j++){
            d2s[j] = dxs[j] * dxs[j] + dys[j] * dys[j] + dzs[j] * dzs[j] + eps2;
        }
        #pragma opm simd
        for(size_t j=0; j<len; j++){
            ivf_d2s[j] = 1.0 / (d2s[j] * std::sqrt(d2s[j]));
        }
        #pragma opm simd
        for(size_t j=0; j<len; j++){
            stars.vxs[i] += dxs[j] * stars.masses[j] * (gdt * ivf_d2s[j]);
        }
        #pragma opm simd
        for(size_t j=0; j<len; j++){
            stars.vys[i] += dys[j] * stars.masses[j] * (gdt * ivf_d2s[j]);
        }
        #pragma opm simd
        for(size_t j=0; j<len; j++){
            stars.vzs[i] += dzs[j] * stars.masses[j] * (gdt * ivf_d2s[j]);
        }
    }

    #pragma opm simd
    for(size_t i=0; i<len; i++){
        stars.pxs[i] += stars.vxs[i] * dt;
    }
    #pragma opm simd
    for(size_t i=0; i<len; i++){
        stars.pys[i] += stars.vys[i] * dt;
    }
    #pragma opm simd
    for(size_t i=0; i<len; i++){
        stars.pzs[i] += stars.vzs[i] * dt;
    }
}

float calc() {
    float energy = 0;
    size_t len = g_len;
    for(size_t i=0; i<len; i++){
        float v2 = stars.vxs[i] * stars.vxs[i] + stars.vys[i] * stars.vys[i] + stars.vzs[i] * stars.vzs[i];
        energy += stars.masses[i] * v2 / 2;
        #pragma GCC unroll 32
        for(size_t j=0; j<len; j++){
            float dx = stars.pxs[j] - stars.pxs[i];
            float dy = stars.pys[j] - stars.pys[i];
            float dz = stars.pzs[j] - stars.pzs[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            float ivf_d2 = 1.0 / (std::sqrt(d2) * 2);
            energy -= stars.masses[j] * stars.masses[j] * (G * ivf_d2);
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
    printf("-----------opt-------------\n");
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
