/*

Initial energy: -10.980484
Final energy: -10.608335
Time elapsed: 1594 ms

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


struct StarVec{
    std::vector<float> pxs;
    std::vector<float> pys;
    std::vector<float> pzs;
    std::vector<float> vxs;
    std::vector<float> vys;
    std::vector<float> vzs;
    std::vector<float> masses;
};

StarVec stars;

void init() {
    stars.pxs.resize(48);
    stars.pys.resize(48);
    stars.pzs.resize(48);
    stars.vxs.resize(48);
    stars.vys.resize(48);
    stars.vzs.resize(48);
    stars.masses.resize(48);
    // for (int i = 0; i < 48; i++) {
    //     stars.pxs[i] = frand();
    //     stars.pys[i] = frand();
    //     stars.pzs[i] = frand();
    //     stars.vxs[i] = frand();
    //     stars.vys[i] = frand();
    //     stars.vzs[i] = frand();
    //     stars.masses[i] = frand()+1;
    // }

    for (int i = 0; i < 48; i++) {
        stars.pxs[i] = frand();
    }
    for (int i = 0; i < 48; i++) {
        stars.pys[i] = frand();
    }
    for (int i = 0; i < 48; i++) {
        stars.pzs[i] = frand();
    }
    for (int i = 0; i < 48; i++) {
        stars.vxs[i] = frand();
    }
    for (int i = 0; i < 48; i++) {
        stars.vys[i] = frand();
    }
    for (int i = 0; i < 48; i++) {
        stars.vzs[i] = frand();
    }
    for (int i = 0; i < 48; i++) {
        stars.masses[i] = frand();
    }
    #pragma GCC unroll 4
    for (int i = 0; i < 48; i++) {
        stars.masses[i] += 1;
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;

void step() {
    size_t len = stars.pxs.size() / 4 * 4;
    for(size_t i=0; i<len; i++){
        for(size_t j=0; j<len; j++){
            float dx = stars.pxs[j] - stars.pxs[i];
            float dy = stars.pys[j] - stars.pys[i];
            float dz = stars.pzs[j] - stars.pzs[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            float ivf_d2 = 1.0 / (d2 * std::sqrt(d2));
            stars.vxs[i] += dx * stars.masses[j] * (G * dt * ivf_d2);
            stars.vys[i] += dy * stars.masses[j] * (G * dt * ivf_d2);
            stars.vzs[i] += dz * stars.masses[j] * (G * dt * ivf_d2);
        }
    }
    // #pragma opm simd
    // for(size_t i=0; i<len; i++){
    //     stars.pxs[i] += stars.vxs[i] * dt;
    //     stars.pys[i] += stars.vys[i] * dt;
    //     stars.pzs[i] += stars.vzs[i] * dt;
    // }
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
    size_t len = stars.pxs.size() / 4 * 4;
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
