#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

float frand() {
    return (float)std::rand() / RAND_MAX * 2 - 1;
}

const std::size_t  N = 48;
struct Star {
    std::array<float,N> px, py, pz;
    std::array<float,N> vx, vy, vz;
    std::array<float,N> mass;
};



Star stars;


void init() {
    for (std::size_t i = 0; i < N; i++) {
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
const float eps_2 = eps*eps;

void step() {
    for (std::size_t i = 0; i < N; i++) {

        float vx_t = 0;
        float vy_t = 0;
        float vz_t = 0;
        for (std::size_t j = 0; j < N; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps_2;
            d2 *= std::sqrt(d2);
            d2 = 1 / d2;
            vx_t += dx * stars.mass[j] * G * dt * d2;
            vy_t += dy * stars.mass[j] * G * dt * d2;
            vz_t += dz * stars.mass[j] * G * dt * d2;
        }
        stars.vx[i] += vx_t;
        stars.vy[i] += vx_t;
        stars.vz[i] += vx_t;

    }
    for(std::size_t i =0; i<N ; i++){
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }

}

float calc() {
    float energy = 0;
    for(std::size_t i = 0 ; i < N; i++){
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2*0.5;
        for(std::size_t j = 0; j<48 ; j++){
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            d2 =  1/std::sqrt(d2);
            energy -= 0.5 * d2 * stars.mass[j] * stars.mass[i] * G ;

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
