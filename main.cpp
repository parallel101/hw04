#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <iostream>
#include <immintrin.h>
#include <cassert>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

constexpr std::size_t N = 48 ;
constexpr float G = 0.001f;
constexpr float eps = 0.001f;
constexpr float dt = 0.01f;
constexpr float eps_sqr = eps * eps ;
constexpr float G_dt = G * dt ;

template<std::size_t N>
struct Star {
    alignas(64) float px[N] , py[N] , pz[N] ;
    alignas(64) float vx[N] , vy[N] , vz[N] ;
    alignas(64) float mass[N];
};

Star<N> stars;

void init(){
    for(std::size_t i {} ; i < N ; ++i) {
        stars.px[i] = frand() , stars.py[i] = frand() , stars.pz[i] = frand() ;
        stars.vx[i] = frand() , stars.vy[i] = frand() , stars.vz[i] = frand() ;
        stars.mass[i] = frand() + 1;
    }
}

using fsimd_t = __m256;
inline auto set1    (float f)               { return _mm256_set1_ps(f);}
inline auto load    (const float * f)       { return _mm256_load_ps(f);}
inline void store   (float * f , fsimd_t a) {        _mm256_store_ps(f , a);}
inline auto add     (fsimd_t a , fsimd_t b) { return _mm256_add_ps(a ,b) ;}
inline auto sub     (fsimd_t a , fsimd_t b) { return _mm256_sub_ps(a ,b) ;}
inline auto mul     (fsimd_t a , fsimd_t b) { return _mm256_mul_ps(a ,b) ;}
inline auto div     (fsimd_t a , fsimd_t b) { return _mm256_div_ps(a ,b) ;}
inline auto sqrt    (fsimd_t a)             { return _mm256_sqrt_ps(a); }
inline auto rsqrt   (fsimd_t a)             { return _mm256_rsqrt_ps(a);}
inline auto rcp     (fsimd_t a)             { return _mm256_rcp_ps(a);  }

void step_avx() {

    auto eps_sqr8 = set1(eps_sqr);

    for(std::size_t j{} ; j < N ; ++j ){
        auto mg_dt = set1(G * stars.mass[j] * dt);
        auto xj = set1(stars.px[j]);
        auto yj = set1(stars.py[j]);
        auto zj = set1(stars.pz[j]);

        auto unroll_body = [&](std::size_t i) {
            auto xi = load(&stars.px[i]);
            auto yi = load(&stars.py[i]);
            auto zi = load(&stars.pz[i]);

            auto dx = sub(xj , xi);
            auto dy = sub(yj , yi);
            auto dz = sub(zj , zi);

            auto x2 = mul(dx , dx);
            auto y2 = mul(dy , dy);
            auto z2 = mul(dz , dz);

            auto d2 = add(add(x2 , y2), add(z2 , eps_sqr8));
            auto inv_d2 = rcp(d2);
            auto inv_d = rsqrt(d2);
            auto mg_dt_invd3 = mul(mul(mg_dt , inv_d2) , inv_d);

            auto vx = load(&stars.vx[i]);
            auto vy = load(&stars.vy[i]);
            auto vz = load(&stars.vz[i]);

            auto new_vx = add(mul(mg_dt_invd3 , dx) , vx);
            auto new_vy = add(mul(mg_dt_invd3 , dy) , vy);
            auto new_vz = add(mul(mg_dt_invd3 , dz) , vz);

            store(&stars.vx[i] , new_vx);
            store(&stars.vy[i] , new_vy);
            store(&stars.vz[i] , new_vz);
        };

        unroll_body(0);
        unroll_body(8);
        unroll_body(16);
        unroll_body(24);
        unroll_body(32);
        unroll_body(40);

        // for(std::size_t i{} ; i < N ; i += 8){
        //     unroll_body(i);
        // }
    }

    for(std::size_t i {} ; i < N; i += 8) {
        auto dt8 = set1(dt);
        auto vx = load(&stars.vx[i]);
        auto vy = load(&stars.vy[i]);
        auto vz = load(&stars.vz[i]);

        auto new_px = add(load(&stars.px[i]) , mul(vx , dt8));
        auto new_py = add(load(&stars.py[i]) , mul(vy , dt8));
        auto new_pz = add(load(&stars.pz[i]) , mul(vz , dt8));

        store(&stars.px[i] , new_px);
        store(&stars.py[i] , new_py);
        store(&stars.pz[i] , new_pz);
    }
}

float calc() {
    float energy = 0;
    for(std::size_t i {} ; i < N;  ++i ){
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2; 
        
        for(std::size_t j {} ; j < N ; ++j) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps ;
            energy -= stars.mass[j] * stars.mass[i] *  G / std::sqrt(d2) / 2 ;
        }
    }
    return energy;
}

template <class Func>
auto benchmark(Func const &func) {
    auto t0 = std::chrono::steady_clock::now();
    func();
    auto t1 = std::chrono::steady_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    return dt.count();
}

int main() {
    init();
    std::cout << "Initial energy: "<< calc() << std::endl;
    auto dt = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step_avx();
    });
    std::cout << "Final energy: " << calc() << std::endl;
    std::cout << "Time elapsed: " << dt << " ms" << std::endl;
    return 0;
}
