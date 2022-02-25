#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <tuple>
#include <vector>

#define VECTOR_SIZE 48

/* static float frand() {
    float INV_RAND_MAX = 1.0 / RAND_MAX;
    return (float)rand() * (INV_RAND_MAX * 2) - 1;
} */
static std::tuple<float, float, float, float, float, float, float> frand7() {
    float INV_RAND_MAX = 1.0 / RAND_MAX;
    float a = (float)rand() * (INV_RAND_MAX * 2) - 1;
    float b = (float)rand() * (INV_RAND_MAX * 2) - 1;
    float c = (float)rand() * (INV_RAND_MAX * 2) - 1;
    float d = (float)rand() * (INV_RAND_MAX * 2) - 1;
    float e = (float)rand() * (INV_RAND_MAX * 2) - 1;
    float f = (float)rand() * (INV_RAND_MAX * 2) - 1;
    float g = (float)rand() * (INV_RAND_MAX * 2) - 1;
    return std::make_tuple(a, b, c, d, e, f, g);
}

struct alignas(32) Star {
    // float px, py, pz;
    // float vx, vy, vz;
    // float mass;
    float px[VECTOR_SIZE];
    float py[VECTOR_SIZE];
    float pz[VECTOR_SIZE];
    float vx[VECTOR_SIZE];
    float vy[VECTOR_SIZE];
    float vz[VECTOR_SIZE];
    float mass[VECTOR_SIZE];
};

// std::vector<Star> stars;
Star stars;

static void init() {
#pragma GCC unroll 8
    for (size_t i = 0; i < VECTOR_SIZE; i++) {
        // stars.push_back({
        //     frand(), frand(), frand(),
        //     frand(), frand(), frand(),
        //     frand() + 1,
        // });
        // stars.px.emplace_back(frand());
        // stars.py.emplace_back(frand());
        // stars.pz.emplace_back(frand());
        // stars.vx.emplace_back(frand());
        // stars.vy.emplace_back(frand());
        // stars.vz.emplace_back(frand());
        // stars.mass.emplace_back(frand() + 1);
        auto [a, b, c, d, e, f, g] = frand7();
        stars.px[i] = a;
        stars.py[i] = b;
        stars.pz[i] = c;
        stars.vx[i] = d;
        stars.vy[i] = e;
        stars.vz[i] = f;
        stars.mass[i] = g + 1;
    }
}

constexpr float G = 0.001;
constexpr float eps = 0.001;
constexpr float dt = 0.01;
constexpr float eps_times_eps = eps * eps;  // eps * eps

static void step() {
#pragma GCC ivdep
    for (int i = 0; i < VECTOR_SIZE; ++i) {
        for (int j = 0; j < VECTOR_SIZE; ++j) {
            if (i == j)
                continue;
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps_times_eps;
            d2 *= std::sqrt(d2);
            float inv_d2 = 1.0 / d2;
            float tmp = stars.mass[j] * (G * dt) * inv_d2;
            stars.vx[i] += dx * tmp;
            stars.vy[i] += dy * tmp;
            stars.vz[i] += dz * tmp;
        }
    }
#pragma GCC unroll 16
    for (int i = 0; i < VECTOR_SIZE; ++i) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

static float calc() {
    float energy = 0;
#pragma GCC ivdep
    for (int i = 0; i < VECTOR_SIZE; ++i) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        for (int j = 0; j < VECTOR_SIZE; ++j) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps_times_eps;
            float inv_sqrt_d2 = 1 / std::sqrt(d2);
            energy -= stars.mass[j] * stars.mass[i] * (G * 0.5f) * inv_sqrt_d2;
        }
    }
    return energy;
}

template <class Func>
static long benchmark(Func const& func) {
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
