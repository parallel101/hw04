#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

//-march=native 让编译器自动判断当前硬件支持的指令
//强迫编译器在编译期求值
//用 constexpr 函数迫使编译器进行常量折叠！
constexpr float speedup = 1.0 / RAND_MAX;

float frand() {
    return (float)rand() / speedup * 2 - 1;
}
constexpr int length = 48;

//存储在栈上无法动态扩充大小，这就是为什么 vector 
//这种数据结构要存在堆上，而固定长度的 array 可以存在栈上

struct Star {
    float px[length]; 
    float py[length]; 
    float pz[length];
    float vx[length]; 
    float vy[length]; 
    float vz[length];
    float mass[length];
};

//SOA：分离存储多个属性
//不符合面向对象编程 (OOP) 的习惯，但常常有利于性能。又称之为面向数据编程

//AOS：紧凑存储多个属性
//符合一般面向对象编程 (OOP) 的习惯，但常常不利于性能
Star stars;

void init() {
    for (int i = 0; i < 48; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1;
    }
}
//循环中的不变量：挪到外面来
float G = 0.001;
float eps = 0.001;
float dt = 0.01;

void step() {
    size_t len = length;
    float eps2 = eps * eps;
    float  gdt = G * dt;
    #pragma GCC unroll 16
    for (size_t i = 0 ; i < len; i++) {
            float dxs[length];
            float dys[length];
            float dzs[length];
            float d2s[length];
            float ivf_d2s[length];
            #pragma opm simd
            for(size_t j=0; j < len; j++)
            {
                dxs[j] = stars.px[j] - stars.px[i];
            }
            #pragma opm simd
            for(size_t j=0; j < len; j++)
            {
                dys[j] = stars.py[j] - stars.py[i];
            }
            #pragma opm simd
            for(size_t j=0; j < len; j++)
            {
                dzs[j] = stars.pz[j] - stars.pz[i];
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++)
            {
                d2s[j] = dxs[j] * dxs[j] + dys[j] * dys[j] + dzs[j] * dzs[j] + eps2;
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++){
                ivf_d2s[j] = 1.0 / (d2s[j] * std::sqrt(d2s[j]));
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++){
                stars.vx[i] += dxs[j] * stars.mass[j] * (gdt * ivf_d2s[j]);
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++){
                stars.vy[i] += dys[j] * stars.mass[j] * (gdt * ivf_d2s[j]);
            }
            #pragma opm simd
            for(size_t j=0; j<len; j++){
                stars.vz[i] += dzs[j] * stars.mass[j] * (gdt * ivf_d2s[j]);
            }
        }
//pragma omp simd 
//C/C++ 的缺点：指针的自由度过高，允许多个 immutable reference 指向同一个对象
//而 Rust 从语法层面禁止，从而让编译器放心大胆优化。
        #pragma opm simd
        for(size_t i=0; i<len; i++)
        {
            stars.px[i] += stars.vx[i] * dt ;
        }
        #pragma opm simd
        for(size_t i=0; i < len; i++)
        {
            stars.py[i] += stars.vy[i] * dt ;
        }
        #pragma opm simd
        for(size_t i = 0; i < len; i++)
        {
            stars.pz[i] += stars.vz[i] * dt;
        }
}
//结论：不管是编译器还是 CPU，都喜欢顺序的连续访问
float calc() {
    float energy = 0;
    size_t len = length;
    for (size_t i = 0; i < len; i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i]* stars.vy[i]+ stars.vz[i]* stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        #pragma GCC unroll 32
//小的循环体进行 unroll 可能是划算的，但最好不要 unroll 大的循环体，否则会造成指令缓存的压力反而变慢！
        for (size_t j=0; j < len; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            float ivf_d2 = 1.0 / (std::sqrt(d2) * 2);
            energy -= stars.mass[j] * stars.mass[j] * (G * ivf_d2);
            //数学优化：除法变乘法
        }
    }
    return energy;
}
//-ffast-math 选项让 GCC 更大胆地尝试浮点运算的优化，有时能带来 2 倍左右的提升。作为代价，他对 NaN 和无穷大的处理，可能会和 IEEE 标准（腐朽的）规定的不一致。
//如果你能保证，程序中永远不会出现 NaN 和无穷大，那么可以放心打开 -ffast-math。
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