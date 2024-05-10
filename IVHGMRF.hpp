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
	class ivhgmrf
	{
		using vec = std::vector<Type>;
		using mat = std::vector<std::vector<Type>>;
		using mat3D = std::vector<std::vector<std::vector<Type>>>;

	public :
		// ============================================================
		ivhgmrf() {
			this->lambda = 1e-7;
			this->alpha = 1e-4;
			this->gammma2 = 1e-3;
			this->sigma2 = 1e+3;
			this->maxepoch = 1000;
			this->eps = 1e-3;
			this->lambda_rate = 1e-12;
			this->alpha_rate = 5e-8;
			this->gammma2_rate = 5e-8;
		}

		~ivhgmrf() {}

		// ============================================================
		// Processing denoising for IVGMRF
		mat denoising(const mat3D& noise)
		{
			// Image setting
			this->enumerate = static_cast<u32>(noise.size());
			this->rows = static_cast<u32>(noise.at(0).at(0).size());
			this->cols = static_cast<u32>(noise.at(0).size());
			this->n = this->rows * this->cols;
			this->eigen = calc_eigen();
			this->expect = calc_expect(noise);
			this->avg_img = averaging(noise);

			// Varizces for image denoising-
			const auto _eps = this->eps;
			const auto _mepoch = this->maxepoch;
			const mat _noise = std::move(centerling(this->avg_img));
			const mat3D _noise3D = std::move(centerling3D(noise));
			mat u = _noise;
			mat v = _noise;
			mat w = _noise;

			// denoise algorithm
			u32 _epoch;
			for (_epoch = 0; _epoch < _mepoch; ++_epoch) {
				Type error = 0.0;
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
				pred_parameters(_noise3D, u, v, w);
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
		mat get_avg_img() { return this->avg_img; }
		mat get_v() { return this->v_final; }
		mat get_w() { return this->w_final; }

		// ============================================================
	private :
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
		u32 rows;
		u32 cols;
		u32 n;

		// Variable for holding image data
		Type expect;
		mat eigen;
		mat avg_img;
		mat3D centered_imgs;
		mat v_final;
		mat w_final;

		// ================================================================================
		// MAP estimation
		void gauss_seidel_method(const mat& noise, mat& u, mat& v, mat& w, Type& err) noexcept
		{
			const auto _rows = static_cast<i32>(this->rows);
			const auto _cols = static_cast<i32>(this->cols);
			const auto _enumerate = this->enumerate;
			const auto _n = this->n;
			const auto _lambda = this->lambda;
			const auto _alpha = this->alpha;
			const auto _sigma2 = this->sigma2;
			const auto _gammma2 = this->gammma2;

			const Type inv_sigma2 = _enumerate / _sigma2;

			// Update u and v
			for (u32 iter = 0; iter < GSITER; ++iter) {
				err = 0;
				for (i32 i = 0; i < _cols; ++i) {
					for (i32 j = 0; j < _rows; ++j) {
						// Numerator for u and v
						Type u_numerator = inv_sigma2 * noise[i][j] + _gammma2 * v[i][j];
						Type v_numerator = _lambda * u[i][j];
						// Denominator
						Type denominator = _lambda;

						// Square grid graph structure
						if (i + 1 < _cols) {
							u_numerator += _alpha * u[i+1][j];
							v_numerator += _alpha * (u[i][j] + v[i+1][j] - u[i+1][j]);
							denominator += _alpha;
						}
						if (i - 1 >= 0) {
							u_numerator += _alpha * u[i-1][j];
							v_numerator += _alpha * (u[i][j] + v[i-1][j] - u[i-1][j]);
							denominator += _alpha;
						}
						if (j + 1 < _rows) {
							u_numerator += _alpha * u[i][j+1];
							v_numerator += _alpha * (u[i][j] + v[i][j+1] - u[i][j+1]);
							denominator += _alpha;
						}
						if (j - 1 >= 0) {
							u_numerator += _alpha * u[i][j-1];
							v_numerator += _alpha * (u[i][j] + v[i][j-1] - u[i][j-1]);
							denominator += _alpha;
						}
						
						const Type u_current = u_numerator / (denominator + inv_sigma2);
						const Type v_current = v_numerator / (denominator + _gammma2);
						err += std::abs(u[i][j] - u_current);
						u[i][j] = u_current;
						v[i][j] = v_current;
					}
				}
				err /= _n;
			}

			// update w
			for (u32 iter = 0; iter < GSITER; ++iter) {
				for (i32 i = 0; i < _cols; ++i) {
					for (i32 j = 0; j < _rows; ++j) {
						// Numerator and denominator for w
						Type w_numerator = v[i][j];
						Type w_denominator = _lambda;

						// Square grid graph structure
						if (i + 1 < _cols) {
							w_numerator += _alpha * w[i+1][j];
							w_denominator += _alpha;
						}
						if (i - 1 >= 0) {
							w_numerator += _alpha * w[i-1][j];
							w_denominator += _alpha;
						}
						if (j + 1 < _rows) {
							w_numerator += _alpha * w[i][j+1];
							w_denominator += _alpha;
						}
						if (j - 1 >= 0) {
							w_numerator += _alpha * w[i][j-1];
							w_denominator += _alpha;
						}
						w[i][j] = w_numerator / w_denominator;
					}
				}
			}
		}

		// parameters estimation
		void pred_parameters(const mat3D& noise, const mat& u, const mat& v, const mat& w) noexcept
		{
			const auto _rows = this->rows;
			const auto _cols = this->cols;
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
			for (u32 i = 0; i < _cols; ++i) {
				for (u32 j = 0; j < _rows; ++j) {
					const Type first = _lambda + _alpha * _eigen[i][j];
					const Type second = _gammma2 + first;

					const Type psi = first * first / second;
					const Type chi = inv_sigma2 + psi;

					lambda_grad += (0.5 * _gammma2 * _gammma2 * w[i][j] * w[i][j] - 0.5 * u[i][j] * u[i][j]) / _enumerate + 0.5 * (2 / first - 1 / second) / (chi * _sigma2);
					alpha_grad += 0.5 * (2 / first - 1 / second) * _eigen[i][j] / (chi * _sigma2);
					gammma2_grad += 0.5 * v[i][j] * v[i][j] / _enumerate - 0.5 / (chi * _sigma2 * second);
					for (const auto& img : noise) sigma2_strict += (img[i][j] - u[i][j]) * (img[i][j] - u[i][j]) / _enumerate + 1 / chi;
				}
			}
			lambda_grad /= _n;
			alpha_grad /=  _n;
			gammma2_grad /= _n;
			sigma2_strict /= _n;

			// Update parameters
			this->lambda += this->lambda_rate * lambda_grad;
			this->alpha += this->alpha_rate * alpha_grad;
			this->gammma2 += this->gammma2_rate * gammma2_grad;
			this->sigma2 = sigma2_strict;
		}

		// ================================================================================
		// eigenvalue of graph laplacian
		mat calc_eigen() noexcept {
			const auto _rows = this->rows;
			const auto _cols = this->cols;

			mat _eigen(_cols, vec(_rows));
			for (u32 i = 0; i < _cols; ++i) {
				for (u32 j = 0; j < _rows; ++j) {
					_eigen[i][j] = 4 * std::pow(std::sin(0.5 * M_PI * i / _cols), 2) 
						+ 4 * std::pow(std::sin(0.5 * M_PI * j / _rows), 2);
				}
			}
			return _eigen;
		}

		// centerling for 3d images
		mat3D centerling3D(mat3D imgs) noexcept {
			const Type _expect = this->expect;
			for (auto& img : imgs) for (auto& v : img) for (auto& p : v) p -= _expect;
			return imgs;
		}

		// centerling for 2d image
		mat centerling(mat img) noexcept {
			const Type _expect = this->expect;
			for (auto& v : img) for (auto& p : v) p -= _expect;
			return img;
		}

		// centerling cancellation
		mat decenterling(mat img) noexcept {
			const Type _expect = this->expect;
			for (auto& v : img) for (auto& p : v) p += _expect;
			return img;
		}

		// expected value
		Type calc_expect(const mat3D& imgs) noexcept {
			Type _expect = 0.0;
			for (const auto& img : imgs) for (const auto& v : img) for (const auto& p : v) _expect += p;
			return _expect / this->n * this->enumerate;
		}

		// 3D image averaging
		mat averaging(const mat3D& imgs) noexcept {
			const auto _rows = this->rows;
			const auto _cols = this->cols;
			const auto _enumerate = this->enumerate;

			mat _avg_img(_cols, vec(_rows));
			for (auto& img : imgs) {
				for (u32 i = 0; i < _cols; ++i) {
					for (u32 j = 0; j < _rows; ++j) {
						_avg_img[i][j] += img[i][j] / _enumerate;
					}
				}
			}
			return _avg_img;
		}

		// calculate x^T * varLambda * y 
		Type smooth_term(const mat& x, const mat& y) noexcept {
			const auto _rows = this->rows;
			const auto _cols = this->cols;

			Type tmp = 0.0;
			for (u32 i = 0; i < _cols; ++i) {
				for (u32 j = 0; j < _rows; ++j) {
					if (i + 1 < _cols) tmp += (x[i][j] - y[i + 1][j]) * (x[i][j] - y[i + 1][j]);
					if (j + 1 < _rows) tmp += (x[i][j] - y[i][j + 1]) * (x[i][j] - y[i][j + 1]);
				}
			}
			return tmp;
		}
	};
}