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
	class ivhgmrf_od
	{
		using vec = std::vector<Type>;
		using mat = std::vector<std::vector<Type>>;
		using mat3D = std::vector<std::vector<std::vector<Type>>>;

	public:
		// ============================================================
		ivhgmrf_od() {
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

		~ivhgmrf_od() {}

		// ============================================================
		// Processing denoising for IVGMRF
		vec denoising(const mat& noise)
		{
			// Image setting
			this->enumerate = static_cast<u32>(noise.size());
			this->n = static_cast<u32>(noise.at(0).size());
			this->eigen = calc_eigen(noise[0]);
			this->expect = calc_expect(noise);
			this->avg_signal = averaging(noise);

			// Varizces for image denoising-
			const auto _eps = this->eps;
			const auto _mepoch = this->maxepoch;
			const vec _noise = std::move(centerling(this->avg_signal));
			const mat _noise2D = std::move(centerling2D(noise));
			vec u = _noise;
			vec v = _noise;
			vec w = _noise;

			Type error = 0.0;
			gauss_seidel_method(_noise, u, v, w, error);
			pred_parameters(_noise2D, u, v, w);

			// denoise algorithm
			u32 _epoch;
			for (_epoch = 0; _epoch < _mepoch; ++_epoch) {
				gauss_seidel_method(_noise, u, v, w, error);
#ifdef DEBUG_MODE 
				std::cout << "error : " << error << std::endl;
				std::cout << "alpha : " << this->alpha << std::endl;
				std::cout << "lambda : " << this->lambda << std::endl;
				std::cout << "gammma2 : " << this->gammma2 << std::endl;
				std::cout << "sigma2 : " << this->sigma2 << std::endl;
				std::cout << "==============" << std::endl;
#endif
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
		Type get_sigma2() { return this->sigma2; }
		u32 get_epoch() { return this->epoch; }
		vec get_avg_signal() { return this->avg_signal; }
		vec get_v() { return this->v_final; }
		vec get_w() { return this->w_final; }

		// ============================================================
	private:
		// HGMRF parameters
		Type lambda;
		Type alpha;
		Type gammma2;
		Type sigma2;

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
		vec avg_signal;
		mat centered_signals;
		vec v_final;
		vec w_final;

		// ================================================================================
		// MAP estimation
		void gauss_seidel_method(const vec& noise, vec& u, vec& v, vec& w, Type& err) noexcept
		{
			const auto _enumerate = this->enumerate;
			const auto _n = static_cast<i32>(this->n);
			const auto _lambda = this->lambda;
			const auto _alpha = this->alpha;
			const auto _sigma2 = this->sigma2;
			const auto _gammma2 = this->gammma2;

			const Type inv_sigma2 = _enumerate / _sigma2;

			// Update u and v
			for (u32 iter = 0; iter < GSITER; ++iter) {
				err = 0;
				for (i32 i = 0; i < _n; ++i) {
					// Numerator for u and v
					Type u_numerator = inv_sigma2 * noise[i] + _gammma2 * v[i];
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

					// Square grid graph structure
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
			const auto _sigma2 = this->sigma2;
			const auto _gammma2 = this->gammma2;
			const auto _eigen = this->eigen;

			const Type inv_sigma2 = _enumerate / _sigma2;

			// variances for gradient
			Type lambda_grad = 0.0;
			Type alpha_grad = (0.5 * _gammma2 * _gammma2 * smooth_term(w, w) - 0.5 * smooth_term(u, u)) / _enumerate;
			Type gammma2_grad = 0.0;
			Type sigma2_strict = 0.0;

			// parameters gradient estimation
			for (u32 i = 0; i < _n; ++i) {
				const Type first = _lambda + _alpha * _eigen[i];
				const Type second = _gammma2 + first;

				const Type psi = first * first / second;
				const Type chi = inv_sigma2 + psi;

				lambda_grad += (0.5 * _gammma2 * _gammma2 * w[i] * w[i] - 0.5 * u[i] * u[i]) / _enumerate + 0.5 * (2 / first - 1 / second) / (chi * _sigma2);
				alpha_grad += 0.5 * (2 / first - 1 / second) * _eigen[i] / (chi * _sigma2);
				gammma2_grad += 0.5 * v[i] * v[i] / _enumerate - 0.5 / (chi * _sigma2 * second);
				for (const auto& img : noise) sigma2_strict += (img[i] - u[i]) * (img[i] - u[i]) / _enumerate + 1 / chi;
			}
			lambda_grad /= _n;
			alpha_grad /= _n;
			gammma2_grad /= _n;
			sigma2_strict /= _n;

			// Update parameters
			this->lambda += this->lambda_rate * lambda_grad;
			this->alpha += this->alpha_rate * alpha_grad;
			this->gammma2 += this->gammma2_rate * gammma2_grad;
			this->sigma2 = sigma2_strict;
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