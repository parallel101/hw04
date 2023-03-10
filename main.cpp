#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

float frand()
{
    return (float)rand() / RAND_MAX * 2 - 1;
}

static constexpr size_t k_count = 48;

template <size_t N>
struct alignas(16) Stars
{
    std::array<float, N> px;
    std::array<float, N> py;
    std::array<float, N> pz;
    std::array<float, N> vx;
    std::array<float, N> vy;
    std::array<float, N> vz;
    std::array<float, N> mass;
};
Stars<k_count> stars;

void init()
{
    for (int i = 0; i < k_count; i++)
    {
        stars.px[i]   = frand();
        stars.py[i]   = frand();
        stars.pz[i]   = frand();
        stars.vx[i]   = frand();
        stars.vy[i]   = frand();
        stars.vz[i]   = frand();
        stars.mass[i] = frand() + 1;
    }
}

float G   = 0.001;
float eps = 0.001;
float dt  = 0.01;

void step()
{
    float eps2 = eps * eps;
    float gt   = G * dt;

    for (size_t i = 0; i < k_count; ++i)
    {
        float dvx = 0.0f;
        float dvy = 0.0f;
        float dvz = 0.0f;

#pragma clang loop unroll_count(4)
        for (size_t j = 0; j < k_count; ++j)
        {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];

            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);

            float mgtd2 = stars.mass[j] * gt / d2;

            dvx += mgtd2 * dx;
            dvy += mgtd2 * dy;
            dvz += mgtd2 * dz;
        }

        stars.vx[i] += dvx;
        stars.vy[i] += dvy;
        stars.vz[i] += dvz;
    }

#pragma clang loop unroll_count(4)
    for (size_t i = 0; i < k_count; ++i)
    {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc()
{
    float half_g = G * 0.5f;
    float eps2   = eps * eps;

    float energy = 0;

#pragma clang loop unroll_count(4)
    for (size_t i = 0; i < k_count; ++i)
    {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] +
                   stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 * 0.5f;
    }

    for (size_t i = 0; i < k_count; ++i)
    {
        float mg = stars.mass[i] * half_g;
        
#pragma clang loop unroll_count(4)
        for (size_t j = 0; j < k_count; ++j)
        {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];

            float d2 = dx * dx + dy * dy + dz * dz + eps2;

            energy -= stars.mass[j] * mg / std::sqrt(d2);
        }
    }
    return energy;
}

template <class Func>
long benchmark(const Func& func)
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
    auto dt = benchmark(
        [&]
        {
            for (int i = 0; i < 100000; i++) step();
        });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}