#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <x86intrin.h>
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


float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

__declspec(align(16)) float px[48],py[48],pz[48];
__declspec(align(16)) float vx[48],vy[48],vz[48];
__declspec(align(16)) float mass[48];

void init() {
    for (uint32_t i = 0; i < 48; i++) {
        px[i] = frand();py[i]=frand();pz[i] = frand();
        vx[i] = frand();vy[i]=frand();vz[i]=frand();
        mass[i] = frand()+1;
    }
}

constexpr float G = 0.001;
constexpr float eps = 0.001;
constexpr float dt = 0.01;

void step() {
    float dx,dy,dz,d2;
    for (size_t i=0;i<(uint32_t)48;++i) {
        for (size_t j=0;j<(uint32_t)48;++j) {
            #pragma opm simd
            dx = px[j] - px[i];
            dy = py[j] - py[i];
            dz = pz[j] - pz[i];
            d2 = dx * dx + dy * dy + dz * dz + (eps * eps);
            d2 *= sqrt(d2);
            d2 = mass[j] * G * dt / d2;
            vx[i] += dx * d2;
            vy[i] += dy * d2;
            vz[i] += dz * d2;
        }
    }
    for(size_t i=0;i<48; ++i){
        #pragma opm simd
        px[i] += vx[i] * dt;
        py[i] += vy[i] * dt;
        pz[i] += vz[i] * dt;
    }
}

float calc() {
    float dx,dy,dz,d2;
    float energy = 0;
    for (size_t i=0;i<48;++i) {
        #pragma opm simd
        float v2 = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
        energy += mass[i] * v2 * 0.5;
        for (size_t j=0;j<48;++j) {
            #pragma opm simd
            dx = px[j] - px[i];
            dy = py[j] - py[i];
            dz = pz[j] - pz[i];
            d2 =  (dx * dx + dy * dy + dz * dz + (eps * eps));
            energy -= mass[j] * mass[i] * 0.0005 / sqrt(d2);
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
