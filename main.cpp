#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>

constexpr float rand_tmp = (float)1/RAND_MAX*2;
float frand() {
    return (float)rand()*rand_tmp - 1;
}

struct alignas(32) Star {
    float px[48], py[48], pz[48];
    float vx[48], vy[48], vz[48];
    float mass[48];
};

Star stars;

void init() {
    float rand_array[48*7];
    for(size_t i=0;i<48*7;i++)
        rand_array[i]=frand();
    for(size_t i=0;i<48;i++)
    {
        stars.px[i]=rand_array[i];
    }
    for(size_t i=0;i<48;i++)
    {
        stars.py[i]=rand_array[48+i];
    }
    for(size_t i=0;i<48;i++)
    {
        stars.pz[i]=rand_array[48*2+i];
    }
    for(size_t i=0;i<48;i++)
    {
        stars.vx[i]=rand_array[48*3+i];
    }
    for(size_t i=0;i<48;i++)
    {
        stars.vy[i]=rand_array[48*4+i];
    }
    for(size_t i=0;i<48;i++)
    {
        stars.vz[i]=rand_array[48*5+i];
    }
    for(size_t i=0;i<48;i++)
    {
        stars.mass[i]=rand_array[48*6+i]+1;
    }
}

constexpr const float G = 0.001;
constexpr const float eps = 0.001;
constexpr const float dt = 0.01;

void step() {
    float Gdt=G * dt;
    #pragma omp simd
    for (size_t i=0;i<48;i++) {
        float dx[48],dy[48],dz[48],d2[48],tmp[48];
        #pragma omp simd
        for(size_t j=0;j<48;j++)
        {
            dx[j] = stars.px[j] - stars.px[i];
            dy[j] = stars.py[j] - stars.py[i];
            dz[j] = stars.pz[j] - stars.pz[i];
        }
        #pragma omp simd
        for(size_t j=0;j<48;j++)
        {
            d2[j]=eps * eps;
            d2[j] += dx[j] * dx[j];
            d2[j] += dy[j] * dy[j];
            d2[j] += dz[j] * dz[j];
        }
        #pragma omp simd
        for(size_t j=0;j<48;j++)
        {
            d2[j] *= std::sqrt(d2[j]);
            tmp[j] = 1/d2[j];
        }

        #pragma omp simd
        for (size_t j=0;j<48;j++) {
            stars.vx[i] += dx[j] * stars.mass[j] * Gdt * tmp[j];
            stars.vy[i] += dy[j] * stars.mass[j] * Gdt * tmp[j];
            stars.vz[i] += dz[j] * stars.mass[j] * Gdt * tmp[j];
        }
    }
    #pragma omp simd
    for (size_t i=0;i<48;i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    float v2[48],dx[48],dy[48],dz[48],d2[48],tmp[48];
    #pragma omp simd
    for(size_t i=0;i<48;i++)
    {
        v2[i] = stars.vx[i] * stars.vx[i];
        v2[i] += stars.vy[i] * stars.vy[i];
        v2[i] += stars.vz[i] * stars.vz[i];
    }
    float tmp_energy=0;
    for(size_t i=0;i<48;i++)
    {
        tmp_energy += stars.mass[i] * v2[i];
    }
    energy+=tmp_energy*0.5;
    #pragma omp simd
    for (size_t i=0;i<48;i++) {
        #pragma omp simd
        for (size_t j=0;j<48;j++) {
            dx[j] = stars.px[j] - stars.px[i];
            dy[j] = stars.py[j] - stars.py[i];
            dz[j] = stars.pz[j] - stars.pz[i];
        }
        #pragma omp simd
        for(size_t j=0;j<48;j++)
        {
            d2[j] = eps * eps; 
            d2[j] += dx[j] * dx[j];
            d2[j] += dy[j] * dy[j];
            d2[j] += dz[j] * dz[j];
        }
        #pragma omp simd
        for(size_t j=0;j<48;j++)
        {
            tmp[j] = 1/std::sqrt(d2[j])*0.5;
            energy -= stars.mass[j] * stars.mass[i] * G * tmp[j];
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
        for (size_t i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
/*
baseline
Initial energy: -8.571528
Final energy: -8.511633
Time elapsed: 1546 ms

simd optimized
Initial energy: -9.936085
Final energy: -9.926659
Time elapsed: 198 ms
*/
