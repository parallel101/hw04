#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <array>
#include<omp.h>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct Star {
    float px, py, pz;
    float vx, vy, vz;
    float mass;
};

std::vector<Star> stars;

void init() {
    for (int i = 0; i < 48; i++) {
        stars.push_back({
            frand(), frand(), frand(),
            frand(), frand(), frand(),
            frand() + 1,
        });
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;

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

static int SIZE=48;

// SOA
struct StarSOA {
    // std::vector<float> px,py,pz,vx,vy,vz,mass;
    std::array<float,48> px,py,pz,vx,vy,vz,mass;
    // float px[48],py[48],pz[48],vx[48],vy[48],vz[48],mass[48];
};

StarSOA stars_soa;

void initSOA(){
    // vector
    // for(auto i=0;i<48;i++){
    //     stars_soa.px.push_back(stars[i].px);
    //     stars_soa.py.push_back(stars[i].py);
    //     stars_soa.pz.push_back(stars[i].pz);
    //     stars_soa.vx.push_back(stars[i].vx);
    //     stars_soa.vy.push_back(stars[i].vy);
    //     stars_soa.vz.push_back(stars[i].vz);
    //     stars_soa.mass.push_back(stars[i].mass);
    // }

    // array
    for(auto i=0;i<48;i++){
        stars_soa.px[i]=stars[i].px;
        stars_soa.py[i]=stars[i].py;
        stars_soa.pz[i]=stars[i].pz;
        stars_soa.vx[i]=stars[i].vx;
        stars_soa.vy[i]=stars[i].vy;
        stars_soa.vz[i]=stars[i].vz;
        stars_soa.mass[i]=stars[i].mass;
    }
}

void stepSOA(){
    auto size=48;
    for(size_t i=0;i<size;i++){
        // 先读局部变量，累加完，在写入
        auto px=stars_soa.px[i],py=stars_soa.py[i],pz=stars_soa.pz[i];
        auto vx=stars_soa.vx[i],vy=stars_soa.vy[i],vz=stars_soa.vz[i];
        auto mass=stars_soa.mass[i];
        for(size_t j=0;j<size;j++){
            auto dx=stars_soa.px[j]-px;
            auto dy=stars_soa.py[j]-py;
            auto dz=stars_soa.pz[j]-pz;
            auto d2=dx * dx + dy * dy + dz * dz + eps * eps;
            d2=d2*std::sqrt(d2); //使用std里的函数
            vx+=dx * stars_soa.mass[j] * (G*dt) / d2; //给不变量添加括号
            vy+=dy * stars_soa.mass[j] * (G*dt) / d2;
            vz+=dz * stars_soa.mass[j] * (G*dt) / d2;
        }
        stars_soa.vx[i]=vx,stars_soa.vy[i]=vy,stars_soa.vz[i]=vz;
        stars_soa.px[i]+=vx*dt,stars_soa.py[i]+=vy*dt,stars_soa.pz[i]+=vz*dt;
    }
}


float calcSOA(){
    float energy=0;
    auto size=48;
    for(auto i=0;i<size;i++){
        float v2=stars_soa.vx[i]*stars_soa.vx[i]+stars_soa.vy[i]*stars_soa.vy[i]+stars_soa.vz[i]*stars_soa.vz[i];
        energy+=stars_soa.mass[i]*v2/2;
        float px=stars_soa.px[i], py=stars_soa.py[i], pz=stars_soa.pz[i];
        float mass=stars_soa.mass[i];
        for(auto j=0;j<size;j++){
            float dx=stars_soa.px[j]-px;
            float dy=stars_soa.py[j]-py;
            float dz=stars_soa.pz[j]-pz;
            float d2= dx * dx + dy * dy + dz * dz + eps * eps;
            energy-=stars_soa.mass[j]*mass*(G*0.5) / std::sqrt(d2);
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
    initSOA();

    printf("Initial energy: %f\n", calc());
    auto dt = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);

    printf("Initial energy: %f\n", calcSOA());
    auto dt_soa = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            stepSOA();
    });
    printf("Final energy: %f\n", calcSOA());
    printf("Time elapsed: %ld ms\n", dt_soa);
    return 0;
}
