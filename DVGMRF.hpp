#pragma once

#include "include.hpp"

#define GSITER 2
// #define DEBUG_MODE
// #define SAVE_PARAM_CSV

#ifdef SAVE_PARAM_CSV
	#include "mylib.hpp"
#endif // SAVE_PARAM_CSV


// 異分散の複数枚画像に対するGMRF
namespace DVGMRF
{
	template<typename Type>
	class gmrf
	{
		typedef std::vector<Type> vec;
		typedef std::vector<std::vector<Type>> matrix;
		typedef std::vector<std::vector<std::vector<Type>>> mat3D;

		// =================================================================
	public:
		gmrf()
		{
			setLambda(1e-7);
			setAlpha(1e-4);
			setSigma2(1000);
			setMaxEpoch(10000);
			setEps(1e-3);
			setLambdaRate(1e-12);
			setAlphaRate(5e-07);
		}

		~gmrf()
		{}

		// =================================================================
		// 画像パラメータの初期設定
		void setImg(const mat3D& img)
		{
			this->enumerate = static_cast<int>(img.size());
			this->rows = static_cast<int>(img.at(0).at(0).size());
			this->cols = static_cast<int>(img.at(0).size());
			this->n = rows * cols;
			this->sigma2 = vecInit(this->initdev);
			this->eigenvalue = calcEigenVal();
			this->imgs = centerize3D(img);
			this->avgImg = averaged(this->imgs);
		}

		matrix processBlock(const mat3D& noise)
		{
			setImg(noise);
			if (this->enumerate > 1) setAlphaRate(1e-08);
			// std::cout << this->imgExpect << std::endl;
			matrix mean = avgImg;

#ifdef SAVE_PARAM_CSV
			std::vector<std::string> table = { "epoch","lambda","alpha","sigma","error" };
			mylib::imgio::makeCSV("SVDVCOMP", "DVGMRF_parameter", table);
#endif

			for (epoch = 0; epoch < maxepoch; epoch++) {
				this->gs(imgs, mean);
				if (error < eps) {
					break;
				}
#ifdef SAVE_PARAM_CSV
				std::vector<std::string> info_table = {
					std::to_string(this->epoch),
					mylib::utillity::float_to_string(this->lambda, 11),
					std::to_string(this->alpha),
					std::to_string(this->sigma2[0]),
					std::to_string(this->error)
			};
				mylib::imgio::saveCSV("SVDVCOMP", "DVGMRF_parameter", info_table);
#endif
				this->predParam(imgs, mean);
			}
#ifdef DEBUG_MODE
		for (int k = 0; k < sigma2.size(); k++) {
			std::cout << "sigma2[" << k << "] : " << sigma2[k] << std::endl;
		}
#endif	// DEBUG_MODE

			return decenterize(mean);
		}

		// =================================================================
		// accsessor
		Type getLambda() const
		{
			return lambda;
		}
		void setLambda(const Type _lambda)
		{
			this->lambda = static_cast<Type>(_lambda);
		}

		Type getAlpha() const
		{
			return alpha;
		}
		void setAlpha(Type _alpha)
		{
			this->alpha = static_cast<Type>(_alpha);
		}

		void setSigma2(const Type _Sigma2)
		{
			this->initdev = static_cast<Type>(_Sigma2);
		}
		vec getSigma2() const
		{
			return sigma2;
		}

		void setMaxEpoch(int _maxepoch)
		{
			this->maxepoch = _maxepoch;
		}

		void setEps(const Type _eps)
		{
			this->eps = _eps;
		}

		Type getError() const
		{
			return static_cast<Type>(this->error);
		}

		void setLambdaRate(const Type _lambdaRate)
		{
			this->lambdaRate = static_cast<Type>(_lambdaRate);
		}

		void setAlphaRate(const Type _alphaRate)
		{
			this->alphaRate = static_cast<Type>(_alphaRate);
		}

		matrix getAvgImg()
		{
			return decenterize(this->avgImg);
		}

		int getEpoch()
		{
			return epoch;
		}

		// =================================================================
	private:
		// GMRFパラメータ
		Type lambda;
		Type alpha;
		Type bias;
		Type initdev;
		vec sigma2;

		// 画像推定パラメータ
		int maxepoch;
		Type eps;
		Type error = 0.0;
		Type lambdaRate;
		Type alphaRate;
		Type imgExpect;

		// 画像サイズ
		int enumerate;
		int rows;
		int cols;
		int n;

		// 画像データ保持変数
		matrix eigenvalue;
		matrix avgImg;
		mat3D imgs;

		// 実験用データ
		uint16_t epoch;

		// ========================================================
		// denoising for GMRF!!!
		void gs(const mat3D& noise, matrix& mean)
		{
			Type inv_sigma2 = 0;
			for (int k = 0; k < enumerate; k++) {
				inv_sigma2 += 1 / sigma2[k];
			}

			// ガウス・ザイデル法
			for (int iter = 0; iter < GSITER; iter++) {
				error = 0.0;

				for (int i = 0; i < rows; i++) {
					for (int j = 0; j < cols; j++) {
						// 分母
						Type denominator = inv_sigma2 + lambda;
						// 分子
						Type numerator = 0.0;
						for (int k = 0; k < enumerate; k++) {
							numerator += noise[k][i][j] / sigma2[k];
						}

						// 正方格子のGMRF
						if (i + 1 < rows) {
							denominator += alpha;
							numerator += alpha * mean[i + 1][j];
						}
						if (i - 1 >= 0) {
							denominator += alpha;
							numerator += alpha * mean[i - 1][j];
						}
						if (j + 1 < cols) {
							denominator += alpha;
							numerator += alpha * mean[i][j + 1];
						}
						if (j - 1 >= 0) {
							denominator += alpha;
							numerator += alpha * mean[i][j - 1];
						}

						const Type currentVal = numerator / denominator;
						error += std::abs(mean[i][j] - currentVal);
						mean[i][j] = currentVal;
					}
				}
				error /= n;
			}
		}

		void predParam(const mat3D& noise, const matrix& mean)
		{
			Type inv_sigma2 = 0;
			for (int i = 0; i < enumerate; i++) {
				inv_sigma2 += 1 / sigma2[i];
			}

			Type lambdaGrad = 0.0;
			Type alphaGrad = -(0.5) * smooth_term(mean, mean);

			for (int k = 0; k < enumerate; k++) {
				for (int i = 0; i < this->rows; i++) {
					for (int j = 0; j < this->cols; j++) {
						const Type psi = lambda + alpha * eigenvalue[i][j];
						const Type chi = inv_sigma2 + psi;

						lambdaGrad += -(0.5) * mean[i][j] * mean[i][j] + (0.5) * inv_sigma2 / (chi * psi);
						alphaGrad += (0.5) * eigenvalue[i][j] * inv_sigma2 / (chi * psi);
						sigma2[k] += (noise[k][i][j] - mean[i][j]) * (noise[k][i][j] - mean[i][j]) + 1 / chi;
					}
				}
				sigma2[k] /= n;
			}
			lambdaGrad /= n * enumerate;
			alphaGrad /= n * enumerate;

			// パラメータの更新
			this->lambda += lambdaRate * lambdaGrad;
			this->alpha += alphaRate * alphaGrad;
		}

		// =================================================================
		// グラフラプラシアンの固有値
		// this->(matrix)eigenvalue
		matrix calcEigenVal()
		{
			matrix _eigenvalue = vecInit_square(0.0);
			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					_eigenvalue[i][j] = 4 * std::pow(std::sin(0.5 * M_PI * i / rows), 2)
						+ 4 * std::pow(std::sin(0.5 * M_PI * j / cols), 2);
				}
			}
			return _eigenvalue;
		}

		// 複数枚画像の平均化
		// this->(matrix)corrupted
		matrix averaged(const mat3D& img)
		{
			matrix _corrupted = vecInit_square(0.0);
			for (int k = 0; k < this->enumerate; k++) {
				for (int i = 0; i < this->rows; i++) {
					for (int j = 0; j < this->cols; j++) {
						_corrupted[i][j] += img[k][i][j] / enumerate;
					}
				}
			}

			return _corrupted;
		}

		// ノイズ画像の期待値計算
		// this->(Type)imgExpect
		Type calcExpect(const matrix& img)
		{
			Type _imgExpect = 0.0;
			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					_imgExpect += img[i][j];
				}
			}
			return _imgExpect /= n;
		}

		// ノイズ画像の中心化
		// this->(matrix)denoised
		matrix centerize(matrix& img)
		{
			matrix _centered = vecInit_square(0.0);
			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					_centered[i][j] = img[i][j] - imgExpect;
				}
			}
			return _centered;
		}

		// 中心化
		// this->(mat3D)imgs
		mat3D centerize3D(mat3D _imgs)
		{
			// 平均化画像の期待値を計算
			this->imgExpect = calcExpect(averaged(_imgs));

			// 平均化画像の期待値でそれぞれの画像の中心化
			for (int k = 0; k < this->enumerate; k++) {
				_imgs.at(k) = centerize(_imgs.at(k));
			}
			return _imgs;
		}

		// 中心化解除
		matrix decenterize(matrix img)
		{
			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					img[i][j] += imgExpect;
				}
			}
			return img;
		}

		// x^T * Λ * y の計算
		Type smooth_term(const matrix& x, const matrix& y)
		{
			Type tmp = 0.0;
			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					// リンク
					if (i + 1 < rows) {
						tmp += (x[i][j] - y[i + 1][j]) * (x[i][j] - y[i + 1][j]);
					}
					if (j + 1 < cols) {
						tmp += (x[i][j] - y[i][j + 1]) * (x[i][j] - y[i][j + 1]);
					}
				}
			}
			return tmp;
		}

		// 1次元vectorの初期化
		vec vecInit(Type initValue)
		{
			vec tmp(enumerate, initValue);
			return tmp;
		}

		// 2次元vectorをinitValueで初期化
		matrix vecInit_square(Type initValue) {
			matrix initVec(this->cols, std::vector<Type>(this->rows, initValue));
			return initVec;
		}

		// 3次元vectorをinitValueで初期化
		mat3D vecInit_cube(Type initValue) {
			mat3D initVec(this->enumerate, std::vector<std::vector<Type>>(this->cols, std::vector<Type>(this->rows, initValue)));
			return initVec;
		}
	};
}