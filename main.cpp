#include <array>
#include <cstdio>
#include <chrono>
#include <cmath>

#include <immintrin.h>

float frand() {
    return (float)rand() / float(RAND_MAX) * 2 - 1;
}

constexpr const auto stars_size = 48;
constexpr const auto alignment = 64; // cache lines width

struct alignas(alignment) Stars {
    using Element = std::array<float, stars_size>;

    alignas(alignment) Element vx, vy, vz;
    alignas(alignment) Element px, py, pz;
    alignas(alignment) Element mass;
};

Stars stars;

void init() {
    for (int i = 0; i < stars_size; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1;
    }
}

constexpr float G = 0.001;
constexpr float eps = 0.001;
constexpr float dt = 0.01;

auto epsx = _mm256_set1_ps(eps * eps);
auto Gdt = _mm256_set1_ps(G * dt);
auto dtx = _mm256_set1_ps(dt);
alignas(alignment) __m256 mgdt[stars_size / 8];
alignas(alignment) __m256 gvx[stars_size];
alignas(alignment) __m256 gvy[stars_size];
alignas(alignment) __m256 gvz[stars_size];

float reduce_sum(const __m256& x) noexcept {
    __m256 s0 = _mm256_hadd_ps(x, x);
    s0 = _mm256_hadd_ps(s0, s0);

    __m128 s1 = _mm256_extractf128_ps(s0, 1);
    s1 = _mm_add_ps(_mm256_extractf128_ps(s0, 0), s1);

    return _mm_cvtss_f32(s1);
}

void step() noexcept {
    for (unsigned i = 0; i < stars_size; i += 8) {
        mgdt[i >> 3] = _mm256_mul_ps(Gdt, _mm256_load_ps(&stars.mass[i]));
    }
    for (unsigned i = 0; i < stars_size; ++i) {
        auto vx = _mm256_setzero_ps();
        auto vy = _mm256_setzero_ps();
        auto vz = _mm256_setzero_ps();
        auto px = _mm256_set1_ps(stars.px[i]);
        auto py = _mm256_set1_ps(stars.py[i]);
        auto pz = _mm256_set1_ps(stars.pz[i]);
        for (unsigned j = 0; j < stars_size; j += 8) {
            auto dx = _mm256_sub_ps(_mm256_load_ps(&stars.px[j]), px);
            auto dy = _mm256_sub_ps(_mm256_load_ps(&stars.py[j]), py);
            auto dz = _mm256_sub_ps(_mm256_load_ps(&stars.pz[j]), pz);
            auto d2 = _mm256_fmadd_ps(dx, dx, _mm256_fmadd_ps(dy, dy, _mm256_fmadd_ps(dz, dz, epsx)));
            d2 = _mm256_mul_ps(mgdt[j >> 3], _mm256_rcp_ps(_mm256_mul_ps(d2, _mm256_sqrt_ps(d2))));
            vx = _mm256_fmadd_ps(dx, d2, vx);
            vy = _mm256_fmadd_ps(dy, d2, vy);
            vz = _mm256_fmadd_ps(dz, d2, vz);
        }
        gvx[i] = vx;
        gvy[i] = vy;
        gvz[i] = vz;
    }
    for (unsigned i = 0; i < stars_size; i += 1) {
        stars.vx[i] += reduce_sum(gvx[i]);
        stars.vy[i] += reduce_sum(gvy[i]);
        stars.vz[i] += reduce_sum(gvz[i]);
    }
    for (unsigned i = 0; i < stars_size; i += 8) {
        auto x = _mm256_add_ps(_mm256_load_ps(&stars.px[i]), _mm256_mul_ps(_mm256_load_ps(&stars.vx[i]), dtx));
        auto y = _mm256_add_ps(_mm256_load_ps(&stars.py[i]), _mm256_mul_ps(_mm256_load_ps(&stars.vy[i]), dtx));
        auto z = _mm256_add_ps(_mm256_load_ps(&stars.pz[i]), _mm256_mul_ps(_mm256_load_ps(&stars.vz[i]), dtx));
        _mm256_store_ps(&stars.px[i], x);
        _mm256_store_ps(&stars.py[i], y);
        _mm256_store_ps(&stars.pz[i], z);
    }
}

float calc() {
    float energy = 0;
    for (auto i = 0; i < stars_size; ++i) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;

        for (auto j = 0; j < stars_size; ++j) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars.mass[j] * stars.mass[i] * G / sqrt(d2) / 2;
        }
    }
    return energy;
}

template <class Func>
long benchmark(Func const& func) {
    auto t0 = std::chrono::steady_clock::now();
    func();
    auto t1 = std::chrono::steady_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    return static_cast<long>(dt.count());
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
