#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>
#include <immintrin.h>
#include <xmmintrin.h>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct Star {
    float px, py, pz;
    float vx, vy, vz;
    float mass;
};

constexpr int NSTAR = 48;
struct StarSOA
{
    float data[7][NSTAR];
    inline float* operator[] (size_t pos) { return data[pos];}
    inline const float* operator[] (const size_t pos) const { return &(*data[pos]);}
};
StarSOA stars2;
std::vector<Star> stars;

float G = 0.001;
float eps = 0.001;
float dt = 0.01;



void init() {
    for (int i = 0; i < 48; i++) {
        stars.push_back({
            frand(), frand(), frand(),
            frand(), frand(), frand(),
            frand() + 1,
        });
    }
}

void step() {
    for (auto &star: stars) {
        for (auto &other: stars) {
            float dx = other.px - star.px;
            float dy = other.py - star.py;
            float dz = other.pz - star.pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            d2 *= sqrt(d2);
            star.vx += dx * other.mass * G * dt / d2;
            star.vy += dy * other.mass * G * dt / d2;
            star.vz += dz * other.mass * G * dt / d2;
        }
    }
    for (auto &star: stars) {
        star.px += star.vx * dt;
        star.py += star.vy * dt;
        star.pz += star.vz * dt;
    }
}


float calc() {
    float energy = 0;
    for (auto &star: stars) {
        float v2 = star.vx * star.vx + star.vy * star.vy + star.vz * star.vz;
        energy += star.mass * v2 / 2;
        for (auto &other: stars) {
            float dx = other.px - star.px;
            float dy = other.py - star.py;
            float dz = other.pz - star.pz;
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= other.mass * star.mass * G / sqrt(d2) / 2;
        }
    }
    return energy;
}
void init2() {
    // for (int i=0; i<NSTAR; i++) {
    //     for (int j=5; j>=0; j--) {
    //         stars2[j][i] = frand();
    //     }
    //     stars2[6][i] = frand() + 1;
    // }
    for (int i=0; i<NSTAR; i++) {
        stars2[0][i] = stars[i].px;
        stars2[1][i] = stars[i].py;
        stars2[2][i] = stars[i].pz;

        stars2[3][i] = stars[i].vx;
        stars2[4][i] = stars[i].vy;
        stars2[5][i] = stars[i].vz;

        stars2[6][i] = stars[i].mass;
    }
}

// 借鉴雷神之锤3的速算法
inline float Q_rsqrt(float number)
{
    constexpr float threehalfs = 1.5F;
    float x2 = number * 0.5F;
    float y = number;
    int i = *reinterpret_cast<int*>(&y);           // evil floating point bit level hacking
    i = 0x5f3759df - (i >> 1); 
    y = *reinterpret_cast<float*>(&i);
    y = y * (threehalfs - (x2 * y * y)); // 1st iteration
    y = y * (threehalfs - (x2 * y * y)); // 2nd iteration, this can be removed

    return y;
}


// #define HAVE_AVX512
#ifdef HAVE_AVX512
#define NSTEP 16
#define LOAD_PS(_P) _mm512_load_ps(_P)
#define STORE_PS(_P, _F) _mm512_store_ps(_P, _F)
#define SET1_PS(_F) _mm512_set1_ps(_F)
#define MUL_PS(_a, _b) _mm512_mul_ps(_a, _b)
#define SUB_PS(_a, _b) _mm512_sub_ps(_a, _b)
#define ADD_PS(_a, _b) _mm512_add_ps(_a, _b)
#define DIV_PS(_a, _b) _mm512_div_ps(_a, _b)
#define RSQRT_PS(_a) _mm512_rsqrt14_ps(_a)
#else
#define NSTEP 8
#define LOAD_PS(_P) _mm256_load_ps(_P)
#define STORE_PS(_P, _F) _mm256_store_ps(_P, _F)
#define SET1_PS(_F) _mm256_set1_ps(_F)
#define MUL_PS(_a, _b) _mm256_mul_ps(_a, _b)
#define SUB_PS(_a, _b) _mm256_sub_ps(_a, _b)
#define ADD_PS(_a, _b) _mm256_add_ps(_a, _b)
#define DIV_PS(_a, _b) _mm256_div_ps(_a, _b)
#define RSQRT_PS(_a) _mm256_rsqrt_ps(_a)
#endif

void step2() {
    auto other = stars2;
    for (int i=0; i<NSTAR; i+=NSTEP) {
        auto star_px = LOAD_PS(&stars2[0][i]);
        auto star_py = LOAD_PS(&stars2[1][i]);
        auto star_pz = LOAD_PS(&stars2[2][i]);

        auto vxx = LOAD_PS(&stars2[3][i]);
        auto vyy = LOAD_PS(&stars2[4][i]);
        auto vzz = LOAD_PS(&stars2[5][i]);

        for (int j=0; j<NSTAR; j++) {
            auto other_px = SET1_PS(other[0][j]);
            auto other_py = SET1_PS(other[1][j]);
            auto other_pz = SET1_PS(other[2][j]);

            auto dxx = SUB_PS(other_px, star_px);
            auto dyy = SUB_PS(other_py, star_py);
            auto dzz = SUB_PS(other_pz, star_pz);

            auto dd2 = MUL_PS(dxx, dxx) + MUL_PS(dyy, dyy) + 
                       MUL_PS(dzz, dzz) + SET1_PS(eps*eps);

            dd2 = DIV_PS(RSQRT_PS(dd2), dd2);
            auto inter = SET1_PS(other[6][j]*(G*dt));
            inter = MUL_PS(inter, dd2);

            vxx = ADD_PS(vxx, MUL_PS(dxx, inter));
            vyy = ADD_PS(vyy, MUL_PS(dyy, inter));
            vzz = ADD_PS(vzz, MUL_PS(dzz, inter));
        }
        STORE_PS(&stars2[3][i], vxx);
        STORE_PS(&stars2[4][i], vyy);
        STORE_PS(&stars2[5][i], vzz);

        static auto dtt = SET1_PS(dt);
        vxx = MUL_PS(vxx, dtt);
        vyy = MUL_PS(vyy, dtt);
        vzz = MUL_PS(vzz, dtt);

        star_px =  ADD_PS(star_px, vxx);
        star_py =  ADD_PS(star_py, vyy);
        star_pz =  ADD_PS(star_pz, vzz);

        STORE_PS(&stars2[0][i], star_px);
        STORE_PS(&stars2[1][i], star_py);
        STORE_PS(&stars2[2][i], star_pz);
    }
}


float calc2() {
    float energy = 0;
    for (int i=0; i<NSTAR; i++) {
        float v2 = stars2[3][i] * stars2[3][i] + stars2[4][i] * stars2[4][i] + stars2[5][i] * stars2[5][i];
        energy += stars2[6][i] * v2 /2;
        for (int j=0; j<NSTAR; j++) {
            float dx = stars2[0][j] - stars2[0][i];
            float dy = stars2[1][j] - stars2[1][i];
            float dz = stars2[2][j] - stars2[2][i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars2[6][j] * stars2[6][i] * G / sqrt(d2) / 2;            
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
    init2();
    printf("original:\n");
    printf("Initial energy: %f\n", calc());
    auto dt = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);

    printf("my:\n");
    printf("Initial energy: %f\n", calc2());
    auto dt2 = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step2();
    });
    printf("Final energy: %f\n", calc2());
    printf("Time elapsed: %ld ms\n", dt2);    
    return 0;
}