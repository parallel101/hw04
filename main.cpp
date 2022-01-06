#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct Star {
	std::vector<float> px;
	std::vector<float> py;
	std::vector<float> pz;
	std::vector<float> vx;
	std::vector<float> vy;
	std::vector<float> vz;
	std::vector<float> mass;
    //float px, py, pz;
    //float vx, vy, vz;
    //float mass;
};

//std::vector<Star> stars;
Star stars;

void init() {
    for (size_t i = 0; i < 48; i++) {
        //stars.push_back({
        //    frand(), frand(), frand(),
        //    frand(), frand(), frand(),
        //    frand() + 1,
        //});
		stars.px.push_back(frand());
		stars.py.push_back(frand());
		stars.pz.push_back(frand());
		stars.vx.push_back(frand());
		stars.vy.push_back(frand());
		stars.vz.push_back(frand());
		stars.mass.push_back(frand() + 1);
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;

void step() {
	size_t star_size = stars.mass.size();
    const float epsMulEps = eps * eps;
    const float dtMulG = dt * G;
	for (size_t i = 0; i < star_size; ++i)
	{
        for (size_t j = 0; j < star_size; ++j)
        {
			float dx = stars.px[j] - stars.px[i];
			float dy = stars.py[j] - stars.py[i];
			float dz = stars.pz[j] - stars.pz[i];
			float d2 = dx * dx + dy * dy + dz * dz + epsMulEps;
			d2 *= std::sqrt(d2);

			stars.vx[i] += dx * stars.mass[j] * dtMulG / d2;
			stars.vy[i] += dy * stars.mass[j] * dtMulG / d2;
			stars.vz[i] += dz * stars.mass[j] * dtMulG / d2;
        }
	}

    for (size_t i = 0; i < star_size; ++i)
    {
		//star.px += star.vx * dt;
		//star.py += star.vy * dt;
		//star.pz += star.vz * dt;
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }

    //for (auto &star: stars) {
    //    for (auto &other: stars) {
    //        float dx = other.px - star.px;
    //        float dy = other.py - star.py;
    //        float dz = other.pz - star.pz;
    //        float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
    //        d2 *= std::sqrt(d2);
    //        star.vx += dx * other.mass * G * dt / d2;
    //        star.vy += dy * other.mass * G * dt / d2;
    //        star.vz += dz * other.mass * G * dt / d2;
    //    }
    //}
    //for (auto &star: stars) {
    //    star.px += star.vx * dt;
    //    star.py += star.vy * dt;
    //    star.pz += star.vz * dt;
    //}
}

float calc() {
    float energy = 0;
	size_t star_size = stars.mass.size();
	const float epsMulEps = eps * eps;
	const float dtMulG = dt * G;

    for (size_t i = 0; i < star_size; ++i)
    {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2.0f;
        for (size_t j = 0; j < star_size; ++j)
        {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + epsMulEps;
            energy -= stars.mass[j] * stars.mass[i] * G / std::sqrt(d2) / 2;
        }
    }

    //for (auto &star: stars) {
    //    float v2 = star.vx * star.vx + star.vy * star.vy + star.vz * star.vz;
    //    energy += star.mass * v2 / 2;
    //    for (auto &other: stars) {
    //        float dx = other.px - star.px;
    //        float dy = other.py - star.py;
    //        float dz = other.pz - star.pz;
    //        float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
    //        energy -= other.mass * star.mass * G / std::sqrt(d2) / 2;
    //    }
    //}
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
