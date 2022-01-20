#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>


//#pragma omp simd //3th: SIMD：没有明显的提升。。。可能是没有打开gcc -fopenmp -O3. 开了也没用。。

constexpr int N = 48;

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

// 2nd: SOA 反而慢了？？？。。。N太小。。
struct Star {
    std::vector<float> px,py,pz;
    std::vector<float> vx,vy,vz;
    std::vector<float> mass;
};
Star stars;

void init() {
    for (int i = 0; i < N; i++) {
        stars.px.push_back(frand());
        stars.py.push_back(frand());
        stars.pz.push_back(frand());
        stars.vx.push_back(frand());
        stars.vy.push_back(frand());
        stars.vz.push_back(frand());
        stars.mass.push_back(frand());
    }
}
constexpr float G = 0.001;
constexpr float eps = 0.000001;
constexpr float dt = 0.01;
constexpr float G_dt = G * dt;

void step()
{
    /*//不要展开这个10W次的大循环，反而变慢
    for (int counter = 0; counter < 100000; counter++) //9th：把循环调用变成调用循环，没用。
    {*/
        //10th: 展开小循环，好像反而慢了一丢丢, 放弃治疗
        #if defined(__GNUC__) || defined(__clang__)
            //#pragma GCC unroll 48
        #elif defined(_MSC_VER)
            #pragma unroll 48
        #endif

        #pragma GCC ivdep
        for (int i = 0; i < N; i++)
        {
            //8th: 嵌套加法优化, 快了一丢丢
            float dvx = 0;
            float dvy = 0;
            float dvz = 0;
        
            for (int j = 0; j < N; j++)
            {
                float dx = stars.px[j] - stars.px[i];
                float dy = stars.py[j] - stars.py[i];
                float dz = stars.pz[j] - stars.pz[i];
                float d2 = dx * dx + dy * dy + dz * dz + eps; //11th: 提前算好eps**2， 节省1次乘法
                d2 *= std::sqrt(d2);                   //7th: 改称std::sqrt配合ffast-math,快了一点点（大受震撼.jpg。 截至目前，快了三倍多。。。
                dvx += dx * stars.mass[j] * G_dt / d2; //1st: 省下一次乘法, 然并卵。
                dvy += dy * stars.mass[j] * G_dt / d2;
                dvz += dz * stars.mass[j] * G_dt / d2;
            }
            stars.vx[i] += dvx;
            stars.vy[i] += dvy;
            stars.vz[i] += dvz;
        }
        #pragma GCC ivdep //12th: 好像快了一星半点？
            
        for (int i = 0; i < N; i++)
        {
            stars.px[i] += stars.vx[i] * dt;
            stars.py[i] += stars.vy[i] * dt;
            stars.pz[i] += stars.vz[i] * dt;
        }
    /*}*/
}

float calc()
{
    float energy = 0;
    for (int i = 0; i < N; i++)
    {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        for (int j = 0; j < N; j++)
        {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps;
            energy -= stars.mass[j] * stars.mass[i] * G / std::sqrt(d2) / 2;
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
        for(int i = 0;i<100000;i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
