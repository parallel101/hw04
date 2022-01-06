#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <immintrin.h>
#include <utility>
#include <array>

constexpr size_t STAR_NUM = 48;
constexpr size_t SIMD_WIDTH = 16;  //16 floats
constexpr size_t STAR_16_NUM = (STAR_NUM + SIMD_WIDTH - 1) / SIMD_WIDTH;

constexpr float G = 0.001;
constexpr float eps = 0.001;
constexpr float dt = 0.01;
constexpr float G_dt = G * dt;

const __m512 eps2 = _mm512_set1_ps(eps * eps);
const __m512 vec_dt = _mm512_set1_ps(dt);
const __m512 vec_G_dt = _mm512_set1_ps(G_dt);

//Layout: AoSoA
struct Star_16 {
	alignas(64) __m512 vx;
	alignas(64) __m512 vy;
	alignas(64) __m512 vz;

	alignas(64) __m512 px;
	alignas(64) __m512 py;
	alignas(64) __m512 pz;

	alignas(64) __m512 mass;
};

std::array<Star_16, STAR_16_NUM> stars;

static float frand() {
	return (float)rand() / RAND_MAX * 2 - 1;
}

void init() {
	for (int i = 0; i < stars.size(); i++) {
		for (int j = 0; j < SIMD_WIDTH; j++) {
			stars[i].px.m512_f32[j] = frand();
			stars[i].py.m512_f32[j] = frand();
			stars[i].pz.m512_f32[j] = frand();
			stars[i].vx.m512_f32[j] = frand();
			stars[i].vy.m512_f32[j] = frand();
			stars[i].vz.m512_f32[j] = frand();
			stars[i].mass.m512_f32[j] = frand() + 1;
		}
	}
}

template<class Fn, size_t... M>
__forceinline static void unroll_impl(Fn fn, size_t L, std::integer_sequence<size_t, M...> iter) {
	constexpr auto S = sizeof...(M);
	if (L == 1) {
		((fn(M)), ...);
	}
	else {
		for (size_t i = 0; i < L; i++) {
			((fn(M + i * S)), ...);
		}
	}
}

template<size_t N, size_t S = N, class Fn>   // N: total iterations, S = iterations per loop. S==N implies unrolling completely
__forceinline constexpr static void UNROLL(Fn fn) {
	static_assert(N % S == 0);
	unroll_impl(fn, N / S, std::make_index_sequence<S>());
}

void step() {
	std::array<__m512, STAR_16_NUM> d_vx{};
	std::array<__m512, STAR_16_NUM> d_vy{};
	std::array<__m512, STAR_16_NUM> d_vz{};
	for (size_t j = 0; j < stars.size(); ++j) {
		auto& star_j = stars[j];
		for (size_t k = 0; k < SIMD_WIDTH; ++k) {
			__m512 px_jk = _mm512_set1_ps(star_j.px.m512_f32[k]);
			__m512 py_jk = _mm512_set1_ps(star_j.py.m512_f32[k]);
			__m512 pz_jk = _mm512_set1_ps(star_j.pz.m512_f32[k]);
			__m512 mass_jk = _mm512_set1_ps(star_j.mass.m512_f32[k]);

			UNROLL<stars.size()>([&](size_t i) {  //unrolling size: 3/1
				__m512 dx = _mm512_sub_ps(px_jk, stars[i].px);
				__m512 dy = _mm512_sub_ps(py_jk, stars[i].py);
				__m512 dz = _mm512_sub_ps(pz_jk, stars[i].pz);
				//float dx = px[j] - px[i];
				//float dy = py[j] - py[i];
				//float dz = pz[j] - pz[i];

				__m512 prod = _mm512_mul_ps(dx, dx);
				prod = _mm512_fmadd_ps(dy, dy, prod);
				prod = _mm512_fmadd_ps(dz, dz, prod);
				prod = _mm512_add_ps(eps2, prod);
				//float prod = dx * dx + dy * dy + dz * dz + eps * eps;

				__m512 inverse_sqrt = _mm512_rsqrt14_ps(prod);  // an approximate version of _mm512_invsqrt_ps 
				__m512 inverse_sqrt2 = _mm512_mul_ps(inverse_sqrt, inverse_sqrt);
				__m512 inverse_sqrt3 = _mm512_mul_ps(inverse_sqrt2, inverse_sqrt);
				//float inverse_sqrt3 = 1 / sqrt(prod)^3;

				__m512 factor = _mm512_mul_ps(inverse_sqrt3, mass_jk);
				//float factor = mass[j] / sqrt(prod)^3;

				d_vx[i] = _mm512_fmadd_ps(dx, factor, d_vx[i]);
				d_vy[i] = _mm512_fmadd_ps(dy, factor, d_vy[i]);
				d_vz[i] = _mm512_fmadd_ps(dz, factor, d_vz[i]);
				//d_vx += dx * factor;
				//d_vy += dy * factor; 
				//d_vz += dz * factor; 
				});
			}
		}

	UNROLL<stars.size()>([&](size_t i) {   //unrolling size: 3/1
		stars[i].vx = _mm512_fmadd_ps(d_vx[i], vec_G_dt, stars[i].vx);
		stars[i].vy = _mm512_fmadd_ps(d_vy[i], vec_G_dt, stars[i].vy);
		stars[i].vz = _mm512_fmadd_ps(d_vz[i], vec_G_dt, stars[i].vz);
		//vx[i] += d_vx * G * dt;
		//vy[i] += d_vy * G * dt;
		//vz[i] += d_vz * G * dt;

		stars[i].px = _mm512_fmadd_ps(stars[i].vx, vec_dt, stars[i].px);
		stars[i].py = _mm512_fmadd_ps(stars[i].vy, vec_dt, stars[i].py);
		stars[i].pz = _mm512_fmadd_ps(stars[i].vz, vec_dt, stars[i].pz);
		//px[i] += vx[i] * dt;
		//py[i] += vy[i] * dt;
		//pz[i] += vz[i] * dt;
		});
}

float calc() {
	float energy = 0.f;
	const float half_G = G * 0.5;
	for (size_t i = 0; i < stars.size(); i++) {
		auto const& star_i = stars[i];

		for (size_t k = 0; k < SIMD_WIDTH; k++) {
			const float px = star_i.px.m512_f32[k];
			const float py = star_i.py.m512_f32[k];
			const float pz = star_i.pz.m512_f32[k];
			const float vx = star_i.vx.m512_f32[k];
			const float vy = star_i.vy.m512_f32[k];
			const float vz = star_i.vz.m512_f32[k];
			const float mass_ik = star_i.mass.m512_f32[k];

			const float v2 = vx * vx + vy * vy + vz * vz;
			energy += mass_ik * v2 / 2;
			float delta_e = 0.f;

			for (size_t j = 0; j < stars.size(); j++) {
				auto const& star_j = stars[j];

				for (size_t m = 0; m < SIMD_WIDTH; m++) {
					if (!(i == j && k == m)) {
						const float dx = star_j.px.m512_f32[m] - px;
						const float dy = star_j.py.m512_f32[m] - py;
						const float dz = star_j.pz.m512_f32[m] - pz;
						const float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
						delta_e += star_j.mass.m512_f32[m] / sqrt(d2);
					}
				}
			}
			energy -= delta_e * mass_ik * half_G;
		}
	}
	return energy;
}

template <class Func>
long long benchmark(Func const& func) {
	auto t0 = std::chrono::high_resolution_clock::now();
	func();
	auto t1 = std::chrono::high_resolution_clock::now();
	auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
	return dt.count();
}

int main() {
	init();
	printf("Initial energy: %f\n", calc());  // Initial energy: -8.571527
	auto const dt = benchmark([&] {
		for (size_t i = 0; i < 100000; i++)
			step();
		});
	printf("Final energy: %f\n", calc());  // Final energy: -8.562095
	printf("Time elapsed: %lld ms\n", dt);
	//system("pause");
	return 0;
}
