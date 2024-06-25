#pragma once

#include "include.hpp"

#define GSITER 2
// #define DEBUG_MODE
// #define SAVE_PARAM_CSV

#ifdef SAVE_PARAM_CSV
#include "mylib.hpp"
#endif // SAVE_PARAM_CSV

using u32 = uint_fast32_t;
using i32 = int_fast32_t;

// Identifier variaces hierachical gaussian markov random feield
namespace HGMRF {
	template<typename Type>
	class dvhgmrf_od
	{
		using vec = std::vector<Type>;
		using mat = std::vector<std::vector<Type>>;
		using mat3D = std::vector<std::vector<std::vector<Type>>>;

	public:
		// ============================================================
		dvhgmrf_od() {
			this->lambda = 1e-11;
			this->alpha = 1e-8;
			this->gammma2 = 1e-8;
			this->sigma2 = 5e-01;
			this->maxepoch = 1000;
			this->eps = 1e-10;
			this->lambda_rate = 1e-13;
			this->alpha_rate = 1e-9;
			this->gammma2_rate = 1e-9;
		}

		~dvhgmrf_od() {}

		// ============================================================
		// Processing denoising for IVGMRF
		vec denoising(const mat& noise)
		{
			// Image setting
			this->enumerate = static_cast<u32>(noise.size());
			this->n = static_cast<u32>(noise[0].size());
			this->eigen = calc_eigen(noise[0]);
			this->expect = calc_expect(noise);
			this->avg_img = averaging(noise);
			const auto _sigma2 = this->sigma2;
			for (u32 k = 0; k < this->enumerate; ++k) this->vec_sigma2.emplace_back(_sigma2);

			// Variace for image denoising
			const auto _eps = this->eps;
			const auto _mepoch = this->maxepoch;
			const auto _noise = std::move(centerling(this->avg_img));
			const auto _noise2D = std::move(centerling2D(noise));
			vec u = _noise;
			vec v = _noise;
			vec w = _noise;

			u32 _epoch;
			Type error = 0.0;
			gauss_seidel_method(_noise2D, u, v, w, error);
			pred_parameters(_noise2D, u, v, w);

			// denoise algorithm
			for (_epoch = 0; _epoch < _mepoch; ++_epoch) {
				gauss_seidel_method(_noise2D, u, v, w, error);
#ifdef DEBUG_MODE 
				std::cout << "error : " << error << std::endl;
				std::cout << "alpha : " << this->alpha << std::endl;
				std::cout << "lambda : " << this->lambda << std::endl;
				std::cout << "gammma2 : " << this->gammma2 << std::endl;
				for (u32 k = 0; k < this->enumerate; ++k) {
					std::cout << "sigma2[" << k << "] : " << this->vec_sigma2[k] << std::endl;
				}
				std::cout << "==============" << std::endl;
#endif
				// std::cout << "error : " << error << std::endl;
				if (error < _eps) break;
				pred_parameters(_noise2D, u, v, w);
			}

			// save variaces
			this->epoch = _epoch;
			this->v_final = std::move(v);
			this->w_final = std::move(w);

			return decenterling(u);
		}

		// ============================================================
		// accessor
		Type get_lambda() { return this->lambda; }
		Type get_alpha() { return this->alpha; }
		Type get_gammma2() { return this->gammma2; }
		vec get_vec_sigma2() { return this->vec_sigma2; }
		u32 get_epoch() { return this->epoch; }
		vec get_avg_img() { return this->avg_img; }
		vec get_v() { return this->v_final; }
		vec get_w() { return this->w_final; }

		// ============================================================
	private:
		// HGMRF parameters
		Type lambda;
		Type alpha;
		Type gammma2;
		Type sigma2;
		vec vec_sigma2;

		// Algorithm parameters
		u32 maxepoch;
		Type eps;
		u32 epoch;

		// parameters for estimating the optimal Gaussian distribution
		Type lambda_rate;
		Type alpha_rate;
		Type gammma2_rate;
		Type sigma2_rate;

		// Image parameters
		u32 enumerate;
		u32 n;

		// Variable for holding image data
		Type expect;
		vec eigen;
		vec avg_img;
		mat centered_signals;
		vec v_final;
		vec w_final;

		// ================================================================================
		// MAP estimation
		void gauss_seidel_method(const mat& noise, vec& u, vec& v, vec& w, Type& err) noexcept
		{
			const auto _enumerate = this->enumerate;
			const auto _n = static_cast<i32>(this->n);
			const auto _lambda = this->lambda;
			const auto _alpha = this->alpha;
			const auto _gammma2 = this->gammma2;
			const auto _vec_sigma2 = this->vec_sigma2;

			Type inv_sigma2 = 0;
			for (auto& vs : _vec_sigma2) inv_sigma2 += 1 / vs;

			// Update u and v
			for (u32 iter = 0; iter < GSITER; ++iter) {
				err = 0;
				for (i32 i = 0; i < _n; ++i) {
					// Numerator for u and v
					Type u_numerator = _gammma2 * v[i];
					for (u32 k = 0; k < _enumerate; ++k) u_numerator += noise[k][i] / _vec_sigma2[k];
					Type v_numerator = _lambda * u[i];
					// Denominator
					Type denominator = _lambda;

					// Square grid graph structure
					if (i + 1 < _n) {
						u_numerator += _alpha * u[i + 1];
						v_numerator += _alpha * (u[i] + v[i + 1] - u[i + 1]);
						denominator += _alpha;
					}
					if (i - 1 >= 0) {
						u_numerator += _alpha * u[i - 1];
						v_numerator += _alpha * (u[i] + v[i - 1] - u[i - 1]);
						denominator += _alpha;
					}

					const Type u_current = u_numerator / (denominator + inv_sigma2);
					const Type v_current = v_numerator / (denominator + _gammma2);
					err += std::abs(u[i] - u_current);
					u[i] = u_current;
					v[i] = v_current;
				}
				err /= _n;
			}

			// update w
			for (u32 iter = 0; iter < GSITER; ++iter) {
				for (i32 i = 0; i < _n; ++i) {
					// Numerator and denominator for w
					Type w_numerator = v[i];
					Type w_denominator = _lambda;

					// path graph structure
					if (i + 1 < _n) {
						w_numerator += _alpha * w[i + 1];
						w_denominator += _alpha;
					}
					if (i - 1 >= 0) {
						w_numerator += _alpha * w[i - 1];
						w_denominator += _alpha;
					}
					w[i] = w_numerator / w_denominator;
				}
			}
		}

		// parameters estimation
		void pred_parameters(const mat& noise, const vec& u, const vec& v, const vec& w) noexcept
		{
			const auto _enumerate = this->enumerate;
			const auto _n = this->n;
			const auto _lambda = this->lambda;
			const auto _alpha = this->alpha;
			const auto _gammma2 = this->gammma2;
			const auto _eigen = this->eigen;
			const auto _vec_sigma2 = this->vec_sigma2;

			Type inv_sigma2 = 0;
			for (auto& vs : _vec_sigma2) inv_sigma2 += 1 / vs;

			// variances for gradient
			Type lambda_grad = 0.0;
			Type alpha_grad = 0.5 * _gammma2 * _gammma2 * smooth_term(w, w) - 0.5 * smooth_term(u, u);
			Type gammma2_grad = 0.0;
			vec sigma2_strict(_enumerate, 0);

			// parameters gradient estimation
			for (u32 i = 0; i < _n; ++i) {
				const Type first = _lambda + _alpha * _eigen[i];
				const Type second = _gammma2 + first;

				const Type psi = first * first / second;
				const Type chi = inv_sigma2 + psi;

				lambda_grad += 0.5 * _gammma2 * _gammma2 * w[i] * w[i] - u[i] * u[i] + 0.5 * (2 / first - 1 / second) * inv_sigma2 * 1 / chi;
				alpha_grad += 0.5 * (2 / first - 1 / second) * _eigen[i] * inv_sigma2 / chi;
				gammma2_grad += 0.5 * v[i] * v[i] - 0.5 * inv_sigma2 / (chi * second);
				for (u32 k = 0; k < _enumerate; ++k) sigma2_strict[k] += (noise[k][i] - u[i]) * (noise[k][i] - u[i]) + 1 / chi;
			}
			lambda_grad /= _n * _enumerate;
			alpha_grad /= _n * _enumerate;
			gammma2_grad /= _n * _enumerate;
			for (auto& sv : sigma2_strict) sv /= _n;

			// Update parameters
			this->lambda += this->lambda_rate * lambda_grad;
			this->alpha += this->alpha_rate * alpha_grad;
			this->gammma2 += this->gammma2_rate * gammma2_grad;
			this->vec_sigma2 = sigma2_strict;
		}

		// ================================================================================
		// eigenvalue of graph laplacian for path graph
		vec calc_eigen(const vec& signal) noexcept {
			const auto _n = signal.size();
			vec _eigen(_n);
			for (u32 i = 0; i < _n; ++i) _eigen[i] = 4 * std::pow(std::sin(0.5 * M_PI * i / _n), 2);
			return _eigen;
		}

		// centerling for 2d signals
		mat centerling2D(mat signals) noexcept {
			const Type _expect = this->expect;
			for (auto& v : signals) for (auto& p : v) p -= _expect;
			return signals;
		}

		// centerling for signal
		vec centerling(vec signal) noexcept {
			const Type _expect = this->expect;
			for (auto& p : signal) p -= _expect;
			return signal;
		}

		// centerling cancellation
		vec decenterling(vec signal) noexcept {
			const Type _expect = this->expect;
			for (auto& p : signal) p += _expect;
			return signal;
		}

		// expected value
		Type calc_expect(const mat& signals) noexcept {
			Type _expect = 0.0;
			for (const auto& v : signals) for (const auto& p : v) _expect += p;
			return _expect / this->n * this->enumerate;
		}

		// signals averaging
		vec averaging(mat signals) noexcept {
			const auto _n = signals[0].size();
			const auto _enumerate = signals.size();

			vec _avg_signal(_n);
			for (const auto& signal : signals) {
				for (u32 i = 0; i < _n; ++i) {
					_avg_signal[i] += signal[i] / _enumerate;
				}
			}
			return _avg_signal;
		}

		// calculate x^T * varLambda * y 
		Type smooth_term(const vec& x, const vec& y) noexcept {
			const auto _n = x.size();
			Type tmp = 0.0;
			for (u32 i = 0; i < _n; ++i) if (i + 1 < _n) tmp += (x[i] - y[i + 1]) * (x[i] - y[i + 1]);
			return tmp;
		}
	};
}