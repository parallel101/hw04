#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

#include <mmintrin.h> //mmx
#include <xmmintrin.h> //sse
#include <emmintrin.h> //sse2
#include <pmmintrin.h> //sse3

#pragma GCC target("avx")
#pragma GCC optimize(3)
#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-ffast-math")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("inline-functions")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("no-stack-protector")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("inline-small-functions")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("inline-functions-called-once")
#pragma GCC optimize("-fdelete-null-pointer-checks")
#pragma G++ ivdep
#pragma G++ unroll 4

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct Star {
    float px, py, pz;
    float vx, vy, vz;
    float mass, mass1;
    Star(float _px,float _py,float _pz,float _vx,float _vy,float _vz,float _mass):
        px(_px),py(_py),pz(_pz),vx(_vx),vy(_vy),vz(_vz),mass(_mass){}
    Star() {}
}stars[48];

void init() {
    for (uint32_t i = 0; i < 48; i++) {
        stars[i]= Star( frand(), frand(), frand(), frand(), frand(), frand(), frand() + 1);
    }
}

constexpr float G = 0.001;
constexpr float eps = 0.001;
constexpr float dt = 0.01;

void step() {
    for (size_t i=0;i<(uint32_t)48;++i) {
        for (size_t j=0;j<(uint32_t)48;++j) {
            #pragma opm simd
            float dx = stars[j].px - stars[i].px;
            float dy = stars[j].py - stars[i].py;
            float dz = stars[j].pz - stars[i].pz;
            float d2 = dx * dx + dy * dy + dz * dz + (eps * eps);
            d2 *= sqrt(d2);
            d2 = stars[j].mass * G * dt / d2;
            stars[i].vx += dx * d2;
            stars[i].vy += dy * d2;
            stars[i].vz += dz * d2;
        }
    }
    for(size_t i=0;i<48; ++i){
        #pragma opm simd
        stars[i].px += stars[i].vx * dt;
        stars[i].py += stars[i].vy * dt;
        stars[i].pz += stars[i].vz * dt;
    }
}

float calc() {
    float energy = 0;
    for (size_t i=0;i<48;++i) {
        #pragma opm simd
        float v2 = stars[i].vx * stars[i].vx + stars[i].vy * stars[i].vy + stars[i].vz * stars[i].vz;
        energy += stars[i].mass * v2 * 0.5;
        for (size_t j=0;j<48;++j) {
            #pragma opm simd
            float dx = stars[j].px - stars[i].px;
            float dy = stars[j].py - stars[i].py;
            float dz = stars[j].pz - stars[i].pz;
            float d2 =  (dx * dx + dy * dy + dz * dz + (eps * eps));
            energy -= stars[j].mass * stars[i].mass * 0.0005 / sqrt(d2);
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
