#include <cstdio>
#include <cstdlib>
#include <array>
#include <chrono>
#include <cmath>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

template <std::size_t N>
struct Stars
{
    std::array<float, N> px, py, pz;
    std::array<float, N> vx, vy, vz;
    std::array<float, N> mass;

    constexpr std::size_t size()
    {
        return N;
    }
};

Stars<48> stars;

void init() {
    const auto size{ stars.size() };
    for (std::size_t i{ }; i < size; ++i)
    {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1;
    }
}

static constexpr float G{ 0.001f };
static constexpr float eps{ 0.001f };
static constexpr float eps2{ eps * eps };
static constexpr float dt{ 0.01f };
static constexpr float G_dt{ G * dt };

void step() {
    const auto size{ stars.size() };
    for (std::size_t i{ }; i < size; ++i)
    {
        auto px{ stars.px[i] };
        auto py{ stars.py[i] };
        auto pz{ stars.pz[i] };
        float vx{ };
        float vy{ };
        float vz{ };
        for (std::size_t j{ }; j < size; ++j)
        {
            auto dx{ stars.px[j] - px };
            auto dy{ stars.py[j] - py };
            auto dz{ stars.pz[j] - pz };
            auto d2{ dx * dx + dy * dy + dz * dz + eps2 };
            d2 *= std::sqrt(d2);
            auto tmp0{ stars.mass[j] * G_dt / d2 };
            vx += dx * tmp0;
            vy += dy * tmp0;
            vz += dz * tmp0;
        }
        stars.vx[i] += vx;
        stars.vy[i] += vy;
        stars.vz[i] += vz;
    }
    for (std::size_t i{ }; i < size; ++i)
    {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    const auto size{ stars.size() };
    for (std::size_t i{ }; i < size; ++i)
    {
        auto px{ stars.px[i] };
        auto py{ stars.py[i] };
        auto pz{ stars.pz[i] };
        auto vx{ stars.vx[i] };
        auto vy{ stars.vy[i] };
        auto vz{ stars.vz[i] };
        auto mass{ stars.mass[i] };
        auto v2 = vx * vx + vy * vy + vz * vz;
        energy += mass * v2 * 0.5f;
        for (std::size_t j{ }; j < size; ++j)
        {
            auto dx{ stars.px[j] - px };
            auto dy{ stars.py[j] - py };
            auto dz{ stars.pz[j] - pz };
            auto d2{ dx * dx + dy * dy + dz * dz + eps2 };
            energy -= stars.mass[j] * mass * G / std::sqrt(d2) * 0.5f;
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
