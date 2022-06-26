#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <immintrin.h>

static float frand() { return (float)std::rand() / RAND_MAX * 2 - 1; }

struct alignas(64) Star {
    float px[48], py[48], pz[48];
    float vx[48], vy[48], vz[48];
    float mass[48];
};

static Star stars[48];

static void init() {
    for (int i = 0; i < 48; i++) {
        stars->px[i] = frand();
        stars->py[i] = frand();
        stars->pz[i] = frand();
        stars->vx[i] = frand();
        stars->vy[i] = frand();
        stars->vz[i] = frand();
        stars->mass[i] = frand() + 1;
    }
}

static constexpr float eps = 0.001;
static constexpr float eps2 = eps * eps;
static constexpr float G = 0.001;
static constexpr float dt = 0.01;
static constexpr float gdt = G * dt;

void step() {
    for (size_t i=0; i<48; ++i) {
        for (size_t j=0; j<48; ++j) {
            float dx = stars->px[j] - stars->px[i];
            float dy = stars->py[j] - stars->py[i];
            float dz = stars->pz[j] - stars->pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);

            stars->vx[i] += dx * stars->mass[j] * gdt / d2;
            stars->vy[i] += dy * stars->mass[j] * gdt / d2;
            stars->vz[i] += dz * stars->mass[j] * gdt / d2;
        }
    }

    for (size_t i=0; i<48; ++i) {
        stars->px[i] += stars->vx[i] * dt;
    }
    for (size_t i=0; i<48; ++i) {
        stars->py[i] += stars->vy[i] * dt;
    }
    for (size_t i=0; i<48; ++i) {
        stars->pz[i] += stars->vz[i] * dt;
    }
}

void step_avx512() {
    const auto eps2_vec = _mm512_set1_ps(eps2);
    const auto gdt_vec = _mm512_set1_ps(gdt);
    for (size_t i=0; i<48; ++i) {
        auto star_px_vec = _mm512_set1_ps (stars->px[i]);
        auto star_py_vec = _mm512_set1_ps (stars->py[i]);
        auto star_pz_vec = _mm512_set1_ps (stars->pz[i]);

        auto vx_vec = _mm512_setzero_ps();
        auto vy_vec = _mm512_setzero_ps();
        auto vz_vec = _mm512_setzero_ps();
        for (size_t j=0; j<48; j+=16) {
            auto other_px_vec = _mm512_loadu_ps(&stars->px[j]);
            auto other_py_vec = _mm512_loadu_ps(&stars->py[j]);
            auto other_pz_vec = _mm512_loadu_ps(&stars->pz[j]);
            auto other_mass_vec = _mm512_loadu_ps(&stars->mass[j]);

            auto dx_vec = _mm512_sub_ps(other_px_vec, star_px_vec);
            auto dy_vec = _mm512_sub_ps(other_py_vec, star_py_vec);
            auto dz_vec = _mm512_sub_ps(other_pz_vec, star_pz_vec);
            auto d2_vec = _mm512_fmadd_ps(dx_vec, dx_vec, _mm512_fmadd_ps(dy_vec, dy_vec, _mm512_fmadd_ps(dz_vec, dz_vec, eps2_vec)));
            auto mul_vec = _mm512_div_ps(_mm512_mul_ps(other_mass_vec, gdt_vec), _mm512_sqrt_ps(d2_vec));

            vx_vec = _mm512_fmadd_ps(dx_vec,mul_vec, vx_vec);
            vy_vec = _mm512_fmadd_ps(dx_vec,mul_vec, vx_vec);
            vz_vec = _mm512_fmadd_ps(dx_vec,mul_vec, vx_vec);
        }

        auto delta_vx = _mm512_reduce_add_ps(vx_vec);
        auto delta_vy = _mm512_reduce_add_ps(vy_vec);
        auto delta_vz = _mm512_reduce_add_ps(vz_vec);

        stars->vx[i] += delta_vx;
        stars->vy[i] += delta_vy;
        stars->vz[i] += delta_vz;
    }

    const auto dt_vec = _mm512_set1_ps(dt);
    for (size_t i=0; i<48; i+=16) {
        auto px_vev = _mm512_loadu_ps(&stars->px[i]);
        auto py_vev = _mm512_loadu_ps(&stars->py[i]);
        auto pz_vev = _mm512_loadu_ps(&stars->pz[i]);

        auto vx_vev = _mm512_loadu_ps(&stars->vx[i]);
        auto vy_vev = _mm512_loadu_ps(&stars->vy[i]);
        auto vz_vev = _mm512_loadu_ps(&stars->vz[i]);

        px_vev = _mm512_fmadd_ps(vx_vev, dt_vec, px_vev);
        py_vev = _mm512_fmadd_ps(vy_vev, dt_vec, py_vev);
        pz_vev = _mm512_fmadd_ps(vz_vev, dt_vec, pz_vev);

        _mm512_storeu_ps(&stars->px[i], px_vev);
        _mm512_storeu_ps(&stars->py[i], py_vev);
        _mm512_storeu_ps(&stars->pz[i], pz_vev);
    }
}

float calc() {
    float energy = 0;
    for (size_t i = 0; i < 48; ++i) {
        float v2 = stars->vx[i] * stars->vx[i] + stars->vy[i] * stars->vy[i] + stars->vz[i] * stars->vz[i];
        energy += stars->mass[i] * v2 / 2;
        for (size_t j = 0; j < 48; ++j) {
            float dx = stars->px[j]- stars->px[i];
            float dy = stars->py[j] - stars->py[i];
            float dz = stars->pz[j] - stars->pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars->mass[j] * stars->mass[i] * G / std::sqrt(d2) / 2;
        }
    }
    return energy;
}

template <class Func> long benchmark(Func const &func) {
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
            step_avx512();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
