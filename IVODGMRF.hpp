#pragma once

#include "include.hpp"

#define GSITER 2
// #define DEBUG_MODE

// 1次元データ(One-Dimentional data)に対するSVGMRF
namespace IVGMRF
{
	template<typename Type>
	class odgmrf
	{
		typedef std::vector<Type> vec;
		typedef std::vector<std::vector<Type>> matrix;

	public:
		// =================================================================
		odgmrf() {
			setLambda(1e-11);
			setAlpha(1e-8);
			setSigma2(5e-01);
			setMaxEpoch(1000);
			setEps(1e-9);
			setLambdaRate(1e-13);
			setAlphaRate(5e-7);
		}

		~odgmrf() {}

		// =================================================================
		void setData(const matrix& _data)
		{
			this->enumerate = static_cast<int>(_data.size());
			this->dataSize = static_cast<int>(_data.at(0).size());
			this->eigenvalue = calcEigenVal();
			this->data = centerize2D(_data);
		}

		void gs(const vec& noise, vec& mean)
		{
			const Type inv_sigma2 = enumerate / sigma2;

			// ガウス・ザイデル法
			for (uint16_t iter = 0; iter < GSITER; iter++) {
				error = 0.0;

				for (uint16_t i = 0; i < dataSize; i++) {
					// 分母
					Type denominator = inv_sigma2 + lambda;
					// 分子
					Type numerator = inv_sigma2 * noise[i];

					// 直線状ののGMRF
					if (i + 1 < dataSize) {
						denominator += alpha;
						numerator += alpha * mean[i + 1];
					}
					if (i - 1 >= 0) {
						denominator += alpha;
						numerator += alpha * mean[i - 1];
					}

					const Type currentVal = numerator / denominator;
					error += std::abs(mean[i] - currentVal);
					mean[i] = currentVal;
				}
				error /= dataSize;
			}
		}

		void predParam(const matrix& noise, vec& mean)
		{
			const Type inv_sigma2 = enumerate / sigma2;

			Type lambdaGrad = 0.0;
			Type alphaGrad = -(0.5) * smooth_term(mean, mean);
			sigma2 = 0.0;

			for (uint16_t i = 0; i < dataSize; i++) {
				const Type psi = lambda + alpha * eigenvalue[i];
				const Type chi = inv_sigma2 + psi;

				lambdaGrad += -(0.5) * mean[i] * mean[i] + (0.5) * inv_sigma2 / (chi * psi);
				alphaGrad += (0.5) * eigenvalue[i] * inv_sigma2 / (chi * psi);
				for (int k = 0; k < enumerate; k++) {
					sigma2 += std::pow(noise[k][i] - mean[i], 2) + 1 / chi;
				}
			}

			lambdaGrad /= dataSize * enumerate;
			alphaGrad /= dataSize * enumerate;

			// パラメータの更新
			this->sigma2 /= dataSize * enumerate;
			this->lambda += lambdaRate * lambdaGrad;
			this->alpha += alphaRate * alphaGrad;
		}

		vec processBlock(const matrix& noise)
		{
			setData(noise);
			this->avgData = averaged(noise);
			vec mean = avgData;

			this->gs(avgData, mean);
			this->predParam(data, mean);

			for (epoch = 0; epoch < maxepoch; epoch++) {
				this->gs(avgData, mean);
				if (error < eps) {
					break;
				}
				this->predParam(data, mean);
			}
#ifdef DEBUG_MODE
			std::cout << "lambda : " << lambda << std::endl;
			std::cout << "alpha : " << alpha << std::endl;
			std::cout << "sigma2 : " << sigma2 << std::endl;
			std::cout << "error/ ep=" << epoch << " :" << error << std::endl;
#endif	// DEBUG_MODE

			return decenterize(mean);
		}

		// =================================================================
		// accessor
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

		Type getSigma2() const
		{
			return sigma2;
		}
		void setSigma2(const Type _sigma2)
		{
			this->sigma2 = static_cast<Type>(_sigma2);
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

		int getEpoch()
		{
			return epoch;
		}

		vec getAvgData()
		{
			return avgData;
		}

	private:
		// =================================================================
			// GMRFパラメータ
		int enumerate;
		Type lambda;
		Type alpha;
		Type sigma2;

		// 推定パラメータ
		int maxepoch;
		Type eps;
		Type error = 0.0;
		Type lambdaRate;
		Type alphaRate;
		Type dataExpect;

		// データサイズ
		int dataSize;

		// データ保持変数
		vec eigenvalue;
		vec avgData;
		matrix data;

		// 実験データ
		int epoch;

		// =================================================================
			// グラフラプラシアンの固有値
		vec calcEigenVal()
		{
			vec _eigenvalue = vecInit();
			for (uint16_t i = 0; i < dataSize; i++)
			{
				_eigenvalue[i] = 4 * std::pow(std::sin(0.5 * M_PI * i / dataSize), 2);
			}
			return _eigenvalue;
		}

		// 複数枚劣化画像の平均化
		// this->corrupted(vec)
		vec averaged(const matrix& _data)
		{
			vec _avgData = vecInit();
			for (uint16_t k = 0; k < enumerate; k++) {
				for (uint16_t i = 0; i < dataSize; i++) {
					_avgData[i] += _data[k][i] / enumerate;
				}
			}
			return _avgData;
		}

		// 期待値の計算
		Type calcExpect(const vec& _data)
		{
			Type _dataExpect = 0.0;
			for (uint16_t i = 0; i < dataSize; i++)
			{
				_dataExpect += _data[i];
			}
			return dataExpect /= dataSize;
		}

		// データの中心化
		// this->denoised(vec)
		vec centerize(vec &_data)
		{
			vec centered = vecInit();
			for (uint16_t i = 0; i < dataSize; i++)
			{
				centered[i] = _data[i] - dataExpect;
			}
			return centered;
		}

		// 複数データの中心化
		// this->data(matrix)
		matrix centerize2D(matrix _data)
		{
			// データの期待値を計算
			this->dataExpect = calcExpect(averaged(_data));

			// 各データを中心化処理
			for (int k = 0; k < enumerate; k++) {
				_data.at(k) = centerize(_data.at(k));
			}
			return _data;
		}

		// 中心化解除
		vec decenterize(vec _data)
		{
			for (uint16_t i = 0; i < dataSize; i++) {
				_data[i] += dataExpect;
			}
			return _data;
		}

		// x^T * Λ * y の計算
		Type smooth_term(const vec& x, const vec& y)
		{
			Type tmp = 0.0;
			for (uint16_t i = 0; i < this->dataSize; i++) {
				if (i + 1 < dataSize) {
					tmp += std::pow(x[i] + y[i + 1], 2);
				}
			}
			return tmp;
		}

		// 1次元vectorの初期化
		vec vecInit()
		{
			vec tmp(dataSize);
			return tmp;
		}
	};
}
