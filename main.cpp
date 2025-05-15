#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

static float frand()
{
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct alignas(32) Star
{
    float px, py, pz;
    float vx, vy, vz;
    float mass;
};

std::vector<Star> stars;

void init()
{
    for (size_t i = 0; i < 48; i++)
    {
        stars.push_back({
            frand(),
            frand(),
            frand(),
            frand(),
            frand(),
            frand(),
            frand() + 1,
        });
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;
float eps2 = eps * eps;
float mm = G * dt;

void step()
{
#pragma GCC ivdep
    // #pragma omp simd
    for (auto &star : stars)
    {
        for (auto &other : stars)
        {
            float dx = other.px - star.px;
            float dy = other.py - star.py;
            float dz = other.pz - star.pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            star.vx += dx * other.mass * mm / d2;
            star.vy += dy * other.mass * mm / d2;
            star.vz += dz * other.mass * mm / d2;
        }
    }
    for (auto &star : stars)

    {
        star.px += star.vx * dt;
        star.py += star.vy * dt;
        star.pz += star.vz * dt;
    }
}

float calc()
{
    float energy = 0;
    for (auto &star : stars)
    {
        float v2 = star.vx * star.vx + star.vy * star.vy + star.vz * star.vz;
        energy += star.mass * v2 * 0.5f;
        for (auto &other : stars)
        {
            float dx = other.px - star.px;
            float dy = other.py - star.py;
            float dz = other.pz - star.pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            energy -= other.mass * star.mass * G / std::sqrt(d2) * 0.5f;
        }
    }
    return energy;
}

template <class Func>
long benchmark(Func const &func)
{
    auto t0 = std::chrono::steady_clock::now();
    func();
    auto t1 = std::chrono::steady_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    return dt.count();
}

int main()
{
    init();
    printf("Initial energy: %f\n", calc());
    auto dt = benchmark([&]
                        {
        for (size_t i = 0; i < 100000; i++)
            step(); });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
