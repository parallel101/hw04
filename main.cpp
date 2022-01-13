#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <chrono>
#include <cmath>

constexpr const float G = 0.001;
constexpr const float eps = 0.001;
constexpr const float dt = 0.01;
constexpr const float multi_eps = 0.001f * 0.001f;
constexpr const int N = 48;

static float frand() {
  constexpr const float tmp = 1.0f / RAND_MAX;
  return (float)std::rand() * tmp  * 2 - 1;
}

struct Star {
  std::array<float, N> px, py, pz;
  std::array<float, N> vx, vy, vz;
  std::array<float, N> mass;
    
};
Star stars;

void init() {
#pragma omp simd 
  for (size_t i = 0; i < N; i++) {
    stars.px[i] = frand();
    stars.py[i] = frand();
    stars.pz[i] = frand();
    stars.vx[i] = frand();
    stars.vy[i] = frand();
    stars.vz[i] = frand();
    stars.mass[i] = frand() + 1.0f;
  }
}


void step() {
#pragma GCC unroll 16
  for(size_t i = 0; i != N; i++){
   float dx{0}, dy{0}, dz{0}, d2{0};
   float vx{0}, vy{0}, vz{0};
    for(size_t j = 0; j != N; j++){
      float mass = stars.mass[j];
      dx = stars.px[j] - stars.px[i];
      dy = stars.py[j] - stars.py[i];
      dz = stars.pz[j] - stars.pz[i];
      d2 = dx * dx + dy * dy + dz * dz + multi_eps;
      d2 *= std::sqrt(d2);
      vx += dx * mass * G * dt / d2;
      vy += dy * mass * G * dt / d2;
      vz += dz * mass * G * dt / d2;
    }
    stars.vx[i] += vx;
    stars.vy[i] += vy;
    stars.vz[i] += vz;
    
  }
    
  for(size_t i = 0; i != N; i++){
    stars.px[i] += stars.vx[i] * dt;
    stars.py[i] += stars.vy[i] * dt;
    stars.pz[i] += stars.vz[i] * dt;
  }


    
    
}

float calc() {

    float energy = 0;
#pragma GCC unroll 8
    for(size_t i = 0; i!=N; i++){
      float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
      energy += stars.mass[i] * v2 / 2;
      float dx{0}, dy{0}, dz{0}, d2{0};
      for(size_t j = 0; j != N; j++){
	 if(i == j) continue;
	dx = stars.px[j] - stars.px[i];
	dy = stars.py[j] - stars.py[i];
	dz = stars.pz[j] - stars.pz[i];
	d2 = dx * dx + dy * dy + dz * dz + multi_eps;
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
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
