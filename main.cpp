#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <array>


const int N = 48;

static float frand() {
    return (float)std::rand() / RAND_MAX * 2 - 1;
}

/* 
struct Star {
    float px, py, pz;
    float vx, vy, vz;
    float mass;
};
*/

struct Star{
	std::array<float, N> px, py, pz;
	std::array<float, N> vx, vy, vz;
	std::array<float, N> mass;
};

Star stars;

/*
void init() {
    for (int i = 0; i < 48; i++) {
        stars.push_back({
            frand(), frand(), frand(),
            frand(), frand(), frand(),
            frand() + 1,
        });
    }
}
*/

void init() {
	for (int i = 0; i < N; i ++) {
		stars.px[i] = frand();
		stars.py[i] = frand();
		stars.pz[i] = frand();
		stars.vx[i] = frand();
		stars.vy[i] = frand();
		stars.vz[i] = frand();
		stars.mass[i] = frand();
	}
}


float G = 0.001;
float eps = 0.001;
float dt = 0.01;

/*
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
*/

void step() {
	const float t = G * dt;
	const float epss = eps * eps;
	for (int i = 0; i < N; i ++) {
		float px = stars.px[i], py = stars.py[i], pz = stars.pz[i];
		float vx = 0.0f, vy = 0.0f, vz = 0.0f;
		for (int j = 0; j < N; j ++) {
			float dx = stars.px[j] - px;
			float dy = stars.py[j] - py;
			float dz = stars.pz[j] - pz;
			float d2 = dx * dx + dy * dy + dz * dz + epss;
			d2 *= std::sqrt(d2);
			float xx = (1 / d2) * t * stars.mass[j];
			vx += dx * xx;
			vy += dy * xx;
			vz += dz * xx;
		}
		stars.vx[i] += vx;
		stars.vy[i] += vy;
		stars.vz[i] += vz; 
	}

	for (int i = 0; i < N; i ++) {
		stars.px[i] += stars.vx[i] * dt;
		stars.py[i] += stars.vy[i] * dt;
		stars.pz[i] += stars.vz[i] * dt;
	}
}

/*
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
*/

float calc() {
	float energy = 0;
	const float epss = eps * eps;
	for (int i = 0; i < N; i ++) {
		float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
		energy += stars.mass[i] * v2 * 0.5f;
		float px = stars.px[i], py = stars.py[i], pz = stars.pz[i];
		for (int j = 0; j < N; j ++) {
			float dx = stars.px[j] - px;
			float dy = stars.py[j] - py;
			float dz = stars.py[j] - pz;
			float d2 = dx * dx + dy * dy + dz * dz + epss;
			float s_d2 = 1 / std::sqrt(d2);
			energy -= stars.mass[j] * stars.mass[i] * G * 0.5; 
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
