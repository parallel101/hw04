#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <iostream>
#include <immintrin.h>
#include <cassert>
#include <Eigen/Eigen>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

constexpr std::size_t N = 48 ;
constexpr float G = 0.001f;
constexpr float eps = 0.001f;
constexpr float dt = 0.01f;
constexpr float eps_sqr = eps * eps ;
constexpr float G_dt = G * dt ;

// using P = Eigen::Vector4f; 
// using V = Eigen::Vector4f;

struct alignas(16) P{
    float x , y , z , w;
    auto data() noexcept -> float * {return &x;}
};

struct alignas(16) V{
    float x , y , z , w; 
    auto data() noexcept -> float * { return &x ;}
};

struct Star {
    alignas(16) P p[N] ;
    alignas(16) V v[N] ;
    float mass[N];

    std::size_t tot{};

    void push_back(
        float _px , float _py , float _pz ,
        float _vx , float _vy , float _vz ,
        float _m){

        auto i = tot++ ;
        p[i] = {_px , _py , _pz , eps};
        v[i] = {_vx , _vy , _vz , 0.f};
        mass[i] = _m;        
    }

    constexpr std::size_t size() const {return N ;}

}stars;

void init() {
    for (std::size_t i {}; i < N ; ++i) {
        stars.push_back(
            frand(), frand(), frand(),
            frand(), frand(), frand(),
            frand() + 1
        );
    }
}

void step_sse() {

    for(std::size_t j{} ; j < N ; ++j ){
        float m_G = stars.mass[j] * G_dt ;
        auto pj = _mm_load_ps((stars.p[j].data()));

        for(std::size_t i {} ; i < N ; ++i ) {
            auto pi = _mm_load_ps((stars.p[i].data()));
            auto v = _mm_load_ps((stars.v[i].data()));

            auto d = _mm_sub_ps(pj , pi);      // pj - pi 
            auto d_sqr = _mm_dp_ps(d, d , 0x7f);
            float d2 = d_sqr.m128_f32[0] + eps_sqr ;
            d2 *= std::sqrt(d2);
            float m_G_inv_d2 = m_G / d2 ;
            
            auto res = _mm_mul_ps(d , _mm_set_ps(0.f , m_G_inv_d2 , m_G_inv_d2 , m_G_inv_d2)) ;
            res = _mm_add_ps(v , res);
            _mm_store_ps((stars.v[i].data()) , res);
        }
    }

    for(std::size_t i {} ; i < N; ++i){
        auto v = _mm_load_ps(stars.v[i].data());
        auto p = _mm_load_ps(stars.p[i].data());
        __m128 t{dt , dt , dt , 0.f};
        _mm_store_ps(stars.p[i].data() , _mm_add_ps(p , _mm_mul_ps(v, t)));
    }
}

float calc() {
    float energy = 0;
    for(std::size_t i {} ; i < N;  ++i ){
        float v2 = stars.v[i].x * stars.v[i].x + stars.v[i].y * stars.v[i].y + stars.v[i].z * stars.v[i].z;
        energy += stars.mass[i] * v2 / 2; 
        
        for(std::size_t j {} ; j < N ; ++j) {
            float dx = stars.p[j].x - stars.p[i].x;
            float dy = stars.p[j].y - stars.p[i].y;
            float dz = stars.p[j].z - stars.p[i].z;
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
            step_sse();
    });
    std::cout << "Final energy: " << calc() << std::endl;
    std::cout << "Time elapsed: " << dt << " ms" << std::endl;
    return 0;
}
