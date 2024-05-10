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


// 同一分散の複数枚画像に対するGMRF
namespace IVGMRF
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
			setMaxEpoch(1000);
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
			this->eigenvalue = calcEigenVal();
			this->imgs = centerize3D(img);
			this->avgImg = averaged(this->imgs);
		}

		// denoising for GMRF!!!
		void gs(const matrix& noise, matrix& mean)
		{
			const Type inv_sigma2 = enumerate / sigma2;

			// ガウス・ザイデル法
			for (int iter = 0; iter < GSITER; iter++) {
				error = 0.0;

				for (int i = 0; i < rows; i++) {
					for (int j = 0; j < cols; j++) {
						// 分母
						Type denominator = inv_sigma2 + lambda;
						// 分子
						Type numerator = inv_sigma2 * noise[i][j];

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
			const Type inv_sigma2 = enumerate / sigma2;

			Type lambdaGrad = 0.0;
			Type alphaGrad = -(0.5) * smooth_term(mean, mean);
			sigma2 = 0.0;

			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					const Type psi = lambda + alpha * eigenvalue[i][j];
					const Type chi = inv_sigma2 + psi;

					lambdaGrad += -(0.5) * mean[i][j] * mean[i][j] + (0.5) * inv_sigma2 / (chi * psi);
					alphaGrad += (0.5) * eigenvalue[i][j] * inv_sigma2 / (chi * psi);
					for (int k = 0; k < enumerate; k++) {
						sigma2 += (noise[k][i][j] - mean[i][j]) * (noise[k][i][j] - mean[i][j]) + 1 / chi;
					}
				}
			}

			lambdaGrad /= n * enumerate;
			alphaGrad /= n * enumerate;

			// パラメータの更新
			this->sigma2 /= n * enumerate;
			this->lambda += lambdaRate * lambdaGrad;
			this->alpha += alphaRate * alphaGrad;
		}

		matrix processBlock(const mat3D& noise)
		{
			setImg(noise);
			// std::cout << this->imgExpect << std::endl;
			matrix mean = avgImg;

#ifdef SAVE_PARAM_CSV
			std::vector<std::string> table = { "epoch","lambda","alpha","sigma","error"};
			mylib::imgio::makeCSV("SVDVCOMP", "SVGMRF_parameter", table);
#endif

			for (epoch = 0; epoch < maxepoch; epoch++) {
				this->gs(avgImg, mean);
				if (error < eps) {
					break;
				}
#ifdef SAVE_PARAM_CSV
				std::vector<std::string> info_table = {
					std::to_string(this->epoch),
					mylib::utillity::float_to_string(this->lambda, 11),
					std::to_string(this->alpha),
					std::to_string(this->sigma2),
					std::to_string(this->error)
				};
				mylib::imgio::saveCSV("SVDVCOMP", "SVGMRF_parameter", info_table);
#endif
				this->predParam(imgs, mean);
			}
#ifdef DEBUG_MODE
			std::cout << "lambda : " << lambda << std::endl;
			std::cout << "alpha : " << alpha << std::endl;
			std::cout << "sigma2 : " << sigma2 << std::endl;
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
		Type sigma2;

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
		int epoch;

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
			for (int k = 0; k < enumerate; k++) {
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

	template<typename Type>
	class upScale
	{
		typedef std::vector<std::vector<Type>> matrix;
		typedef std::vector<std::vector<std::vector<Type>>> mat3D;
		// =================================================================
	public:
		upScale()
		{
			// 画像劣化変数
			setEnumerate(1);
			setM(0.0);
			setStddev(15);

			// 画像修復変数
			setLambda(1e-7);
			setAlpha(1e-4);
			setSigma2(std::pow(this->stddev, 2));
			setMaxEpoch(1000);
			setEps(1e-04);
			setLambdaRate(1e-12);
			setAlphaRate(5e-07);
			setEps_USGS(1e-4);
		}

		~upScale() {}

		virtual void setImg(const matrix& img)
		{
			this->rows = static_cast<int>(img.at(0).size());
			this->cols = static_cast<int>(img.size());
			this->n = rows * cols;
			this->imgs = this->corruption(img);
			this->originImg = img;
			// rows*2 x cols*2の配列を用意
			this->upScalled = vecInit_us();
			// corruptedに複数枚画像を平均化した画像で初期化
			this->corrupted = averaged();
			// 劣化画像に対するグラフラプラシアンの固有値をeigeValに取得
			calcEigenVal();
			// denoisedを中心化した画像で初期化
			this->denoised = centerize(corrupted);
			// 複数枚画像それぞれを中心化
			centerize3D();
		}

		void gs(const matrix& noise, matrix& mean)
		{
			const Type inv_sigma2 = enumerate / sigma2;

			// ガウス・ザイデル法
			for (uint16_t iter = 0; iter < GSITER; iter++) {
				error = 0.0;

				for (int i = 0; i < rows; i++) {
					for (uint16_t j = 0; j < cols; j++) {
						// 分母
						Type denominator = inv_sigma2 + lambda;
						// 分子
						Type numerator = inv_sigma2 * noise[i][j];

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
			const Type inv_sigma2 = enumerate / sigma2;

			Type lambdaGrad = 0.0;
			Type alphaGrad = -(0.5) * smooth_term(mean, mean);

			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					const Type psi = lambda + alpha * eigenvalue[i][j];
					const Type chi = inv_sigma2 + psi;

					lambdaGrad += -(0.5) * mean[i][j] * mean[i][j] + (0.5) * inv_sigma2 / (chi * psi);
					alphaGrad += (0.5) * eigenvalue[i][j] * inv_sigma2 / (chi * psi);
				}
			}

			lambdaGrad /= n * enumerate;
			alphaGrad /= n * enumerate;

			// パラメータの更新
			this->lambda += lambdaRate * lambdaGrad;
			this->alpha += alphaRate * alphaGrad;
		}

		void upScalingGS()
		{
			error = 0.0;
			// 1階層目・・・行方向、列方向にのみ隠れ層を推定
			for (int iter = 0; iter < GSITER; iter++) {
				int count = 0;
				for (int i = 0; i < rows * 2; i++) {
					for (int j = 0; j < cols * 2; j++) {
						Type denominator = lambda;	// 分母
						Type numerator = 0;		// 分子
						// 行方向の推定
						if (i % 2 == 1 && j % 2 == 0) {
							if (i + 1 < rows) {
								denominator += alpha;
								numerator += alpha * upScalled[i + 1][j];
							}
							if (i - 1 >= 0) {
								denominator += alpha;
								numerator += alpha * upScalled[i - 1][j];
							}
							const Type currentVal = numerator / denominator;
							error += std::abs(upScalled[i][j] - currentVal);
							upScalled[i][j] = currentVal;
							count++;
						}
						// 列方向の推定
						if (i % 2 == 0 && j % 2 == 1) {
							if (j + 1 < cols) {
								denominator += alpha;
								numerator += alpha * upScalled[i][j + 1];
							}
							if (j - 1 >= 0) {
								denominator += alpha;
								numerator += alpha * upScalled[i][j - 1];
							}
							const Type currentVal = numerator / denominator;
							error += std::abs(upScalled[i][j] - currentVal);
							upScalled[i][j] = currentVal;
							count++;
						}
					}
				}
				error /= count;
			}

			// 2階層目・・・残りの空白部分を推定
			for (int iter = 0; iter < GSITER; iter++) {
				int count = 0;
				for (int i = 0; i < rows * 2; i++) {
					for (int j = 0; j < cols * 2; j++) {
						Type denominator = lambda;	// 分母
						Type numerator = 0;		// 分子
						if (i % 2 == 1 && j % 2 == 1) {
							if (i + 1 < rows) {
								denominator += alpha;
								numerator += alpha * upScalled[i + 1][j];
							}
							if (i - 1 >= 0) {
								denominator += alpha;
								numerator += alpha * upScalled[i - 1][j];
							}
							if (j + 1 < cols) {
								denominator += alpha;
								numerator += alpha * upScalled[i][j + 1];
							}
							if (j - 1 >= 0) {
								denominator += alpha;
								numerator += alpha * upScalled[i][j - 1];
							}
							const Type currentVal = numerator / denominator;
							error += std::abs(upScalled[i][j] - currentVal);
							upScalled[i][j] = currentVal;
							count++;
						}
					}
				}
				error /= count;
			}

			error = 0;
			// 補完ピクセルの全学習
			while (error > eps_USGS) {
				for (int iter = 0; iter < GSITER; iter++) {
					int count = 0;
					for (int i = 0; i < rows * 2; i++) {
						for (int j = 0; j < cols * 2; j++) {
							Type denominator = lambda;	// 分母
							Type numerator = 0;		// 分子

							if (i % 2 == 0 && j % 2 == 0) {
								continue;
							}

							if (i + 1 < rows) {
								denominator += alpha;
								numerator += alpha * upScalled[i + 1][j];
							}
							if (i - 1 >= 0) {
								denominator += alpha;
								numerator += alpha * upScalled[i - 1][j];
							}
							if (j + 1 < cols) {
								denominator += alpha;
								numerator += alpha * upScalled[i][j + 1];
							}
							if (j - 1 >= 0) {
								denominator += alpha;
								numerator += alpha * upScalled[i][j - 1];
							}
							const Type currentVal = numerator / denominator;
							error += std::abs(upScalled[i][j] - currentVal);
							upScalled[i][j] = currentVal;
							count++;
						}
					}
					error /= count;
				}
			}
		}

		matrix processBlock(matrix& img)
		{
			setImg(img);
			matrix mean = denoised;

			for (epoch = 0; epoch < maxepoch; epoch++) {
				this->gs(denoised, mean);
				if (error < eps) {
					break;
				}
				this->predParam(imgs, mean);
			}
#ifdef DEBUG_MODE
			std::cout << "epoch :" << epoch << std::endl;
			std::cout << "lambda : " << lambda << std::endl;
			std::cout << "alpha : " << alpha << std::endl;
			std::cout << "sigma2 : " << sigma2 << std::endl;
#endif	// DEBUG_MODE

			upScalingGS();

			return decenterize(upScalled, calcExpect(img));
		}

		// =====================================================================
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

		matrix getCorruptedImg() const
		{
			return corrupted;
		}

		int getEpoch()	const
		{
			return epoch;
		}

		void setEnumerate(const int _enumerate)
		{
			this->enumerate = _enumerate;
		}
		void setM(const Type _m)
		{
			this->m = _m;
		}
		void setStddev(const Type _stddev)
		{
			this->stddev = _stddev;
		}
		void setEps_USGS(const Type _eps_USGS)
		{
			this->eps_USGS = _eps_USGS;
		}

		// =========================================================================
	private:
		// GMRFパラメータ
		Type lambda;
		Type alpha;
		Type sigma2;

		// 画像推定パラメータ
		int maxepoch;
		Type eps;
		Type eps_USGS;
		Type error = 0.0;
		Type lambdaRate;
		Type alphaRate;
		Type imgExpect;

		// 画像サイズ
		int rows;
		int cols;
		int n;

		// 画像データ保持変数
		matrix originImg;
		matrix eigenvalue;
		matrix corrupted;
		matrix denoised;
		mat3D imgs;
		matrix upScalled;

		// ノイズ付加用変数
		int enumerate;
		Type m;
		Type stddev;

		// 実験用変数
		int epoch;

		// =================================================================
		// グラフラプラシアンの固有値
		// this->(matrix)eigenvalue
		void calcEigenVal()
		{
			eigenvalue = vecInit_square(0.0);
			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					eigenvalue[i][j] = 4 * std::pow(std::sin(0.5 * M_PI * i / rows), 2)
						+ 4 * std::pow(std::sin(0.5 * M_PI * j / cols), 2);
				}
			}
		}

		// 複数枚画像の平均化
		// this->(matrix)corrupted
		matrix averaged()
		{
			matrix img = vecInit_square(0.0);
			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					for (int k = 0; k < this->enumerate; k++) {
						img[i][j] += imgs[k][i][j] / enumerate;
					}
				}
			}
			return img;
		}

		// 画素の期待値計算
		double calcExpect(matrix& img)
		{
			int imgRow = static_cast<int>(img.at(0).size());
			int imgCol = static_cast<int>(img.size());
			int pixs = imgRow * imgCol;
			double imgEx = 0.0;
			for (int i = 0; i < imgRow; i++) {
				for (int j = 0; j < imgCol; j++) {
					imgEx += img[i][j];
				}
			}
			return imgEx /= pixs;
		}

		// ノイズ画像の中心化
		// this->(matrix)denoised
		matrix centerize(matrix& img)
		{
			matrix centeredImg = vecInit_square(0.0);
			double imgEx = calcExpect(img);
			for (int i = 0; i < this->rows; i++) {
				for (int j = 0; j < this->cols; j++) {
					centeredImg[i][j] = img[i][j] - imgEx;
				}
			}
			return centeredImg;
		}

		// 複数枚ノイズ画像のそれぞれについて中心化
		// this->(mat3D)imgs
		void centerize3D()
		{
			for (int k = 0; k < enumerate; k++) {
				Type E = 0.0;
				for (int i = 0; i < this->rows; i++) {
					for (int j = 0; j < this->cols; j++) {
						E += imgs[k][i][j];
					}
				}
				E /= n;

				for (int i = 0; i < this->rows; i++) {
					for (int j = 0; j < this->cols; j++) {
						imgs[k][i][j] -= E;
					}
				}
			}
		}

		// 中心化解除
		matrix decenterize(matrix& img, double expect)
		{
			int alt_rows = static_cast<int>(img.at(0).size());
			int alt_cols = static_cast<int>(img.size());
			for (int i = 0; i < alt_rows; i++) {
				for (int j = 0; j < alt_cols; j++) {
					img[i][j] += expect;
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

		// 2次元vectorをinitValueで初期化
		matrix vecInit_square(Type initValue)
		{
			matrix initVec(this->rows, std::vector<Type>(this->cols, initValue));
			return initVec;
		}

		// 3次元vectorをinitValueで初期化
		mat3D vecInit_cube(Type initValue)
		{
			mat3D initVec(this->enumerate, std::vector<std::vector<Type>>(this->cols, std::vector<Type>(this->rows, initValue)));
			return initVec;
		}

		// 拡大画像を中心化した原画像で2pix飛ばしで初期化
		matrix vecInit_us()
		{
			matrix us(this->rows * 2, std::vector<Type>(this->cols * 2, 0.0));
			double imgEx = calcExpect(originImg);

			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					us[i * 2][j * 2] = originImg[i][j] - imgEx;
				}
			}
			return us;
		}

		//画像の劣化
		mat3D corruption(const matrix& img)
		{
			// ノイズ画像を格納するvectorを作成
			mat3D noise(enumerate, std::vector<std::vector<double>>(cols, std::vector<double>(rows, 0.0)));

			// ガウス乱数を生成してvectorに格納する
			std::random_device seed_gen;	//シードの初期化
			std::default_random_engine gen(seed_gen());	//乱数生成
			std::normal_distribution<double> dist(this->m, this->stddev);	//ガウス分布

			for (int k = 0; k < enumerate; k++) {
				for (int i = 0; i < rows; i++) {
					for (int j = 0; j < cols; j++) {
						const double tmp = dist(gen);
						noise[k][i][j] = img[i][j] + tmp;
						if (noise[k][i][j] > 255) {
							noise[k][i][j] = 255;
						}
						if (noise[k][i][j] < 0) {
							noise[k][i][j] = 0;
						}
					}
				}
			}
			return noise;
		}
	};
	

}