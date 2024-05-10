#include "mylib.hpp"

using namespace std::filesystem;
using std::system;

namespace
{
	struct parameters
	{
		// =====================================================================
		mylib::imgio::ImgIndex imgName = mylib::imgio::ImgIndex::LENNA;		// 入力画像

		// =====================================================================
		bool imgVisualize = false;				// 画像の表示 ? (画像)
		bool infoVisualize = true;				// 情報の表示 ? (画像, 1次元)
		bool saveImg = true;					// 画像の保存 ?	(画像)
		bool saveCsv = false;					// csvの保存 ? (画像, 1次元)
		bool savePlot = false;					// (1次元データの)プロットを保存 ? (1次元)
		u32 iter = 1;							// 処理の実行回数 (画像, 1次元)

		// =====================================================================
		// SVGMRF //
		u32 sv_enumerate = 3;				// 劣化画像の枚数
		double sv_m = 0.0;						// ガウスノイズの平均
		double sv_stddev = 30.0;				// ガウスノイズの標準偏差

		// DVGMRF //
		u32 dv_enumerate = 4;				// 劣化画像の枚数
		double dv_m = 0.0;				// ガウスノイズの平均
		vec dv_stddev = {30, 30, 30};			// ガウスノイズの標準偏差
		double dv_minStddev = 20;				// ガウスノイズの最小標準偏差(isRandom == true)
		double dv_maxStddev = 40;				// ガウスノイズの最大標準偏差(isRadnom == true)
		bool dv_isRandom = true;			// 劣化画像の分散をランダムに生成(一様分布) ?

		// ======================================================================
		// WAVE DATA //
		u32 sampleRate = 512;											// サンプルレート
		double dataTime = 1.0;											// 生成信号の時間
		u32 dataSize = static_cast<u32>(sampleRate * dataTime);			// 生成信号の総サンプル数
		u32 bufSize = dataSize;											// バッファサイズ
		vec freq = {2,3,5,7,11,13,17};										// 周波数
		double amplitude = 1.0;											// 振幅
		u32 plotSize = bufSize;											// プロットサイズ

		// ======================================================================
		// SVODGMRF //
		u32 vec_enumerate = 1;											// データの個数
		double vec_m = 0.0;												// ガウスノイズの平均
		double vec_stddev = 0.1;										// ガウスノイズの標準偏差

		// DVODGMRF //
		u32 dvod_enumerate = 3;											// ランダムデータの個数
		double dvod_m = 0.0;											// ノイズの平均
		vec dvod_stddev = {0.2};							// ノイズの標準偏差
		double dvod_minStddev = 0.1;									// ランダムノイズの最小標準偏差
		double dvod_maxStddev = 0.3;									// ランダムノイズの最大標準偏差
		bool dvod_isRandom = true;										// ノイズの標準偏差をランダムに生成(一様分布) ?

		// ======================================================================
		// imgio //
		std::string inPath = "D:\\yasuda-lab\\B4\\img\\";					// 入力画像のディレクトリ (abs path)
		std::string outPath = "D:\\yasuda-lab\\M1\\save_img\\";				// 出力画像のディレクトリ (abs path)
		std::string csvName = "200iter-" + std::to_string(dv_enumerate);	// ↑と同じ階層に出力

		// audio //
		std::string wavePlotPath = "D:\\yasuda-lab\\B4\\save_img\\ODGMRF";	// 1次元データのプロット出力先 (abs path)

		//=======================================================================
		u32 progress = 0;
	};

	parameters params;
}

// ===================================================================
// ImgIndexから画像の名前を取得
std::string mylib::imgio::getImgName(ImgIndex imgName)
{
	switch (imgName) {
	case ImgIndex::BOAT:
		return "BOAT";
	case ImgIndex::Clock:
		return "Clock";
	case ImgIndex::girl:
		return "girl";
	case ImgIndex::Jellybeans:
		return "Jelly beans";
	case ImgIndex::LENNA:
		return "LENNA";
	case ImgIndex::Lighthouse:
		return "Lighthouse";
	case ImgIndex::Man:
		return "Man";
	case ImgIndex::barbara512:
		return "barbara512";
	case ImgIndex::APC:
		return "APC";
	case ImgIndex::Stream_and_bridge:
		return "Stream and bridge";
	case ImgIndex::Tank:
		return "Tank";
	case ImgIndex::House:
		return "House1";
	case ImgIndex::Moon_surface:
		return "Moon surface";
	case ImgIndex::Tree:
		return "Tree";
	default:
		std::cout << "not found" << std::endl;
		return "";
	}
}

// 画像をvectorで読み込み
matrix mylib::imgio::loadImg()
{
	cv::Mat img = cv::imread(params.inPath + getImgName(params.imgName) + ".bmp", 0);
	matrix loadv = mylib::img::MtoV(img);
	return loadv;
}

// 画像を保存
void mylib::imgio::saveImg(std::string dir, std::string title, const matrix& img)
{
	// ディレクトリの階層構築
	std::string imgPath = params.outPath + "\\" + mylib::imgio::getImgName(params.imgName) + "\\" + dir;

	// ディレクトリの存在判定と階層の作成
	if (std::filesystem::exists(imgPath) == false) {
		std::filesystem::create_directories(imgPath);
	}
	cv::imwrite(imgPath + "\\" + title + ".jpg", mylib::img::VtoM(img));
}

// 新しいCSVを作成
void mylib::imgio::makeCSV(std::string dir, std::string title, const std::vector<std::string>& table)
{
	if (title == "") return;

	// .csvを保存するディレクトリの作成(初回のみ)
	std::string csvPath = params.outPath + "\\" + mylib::imgio::getImgName(params.imgName) + "\\" + dir;
	std::filesystem::create_directories(csvPath);

	// .csvのテーブルを作成
	std::string _table = "";
	for (u32 i = 0; i < table.size(); i++) {
		_table += table[i] + ",";
	}

	// .csvに書き込み
	std::ofstream outputfile(csvPath + "\\" + title + ".csv", std::ios_base::trunc);
	outputfile << _table << std::endl;
}

// CSVを保存
void mylib::imgio::saveCSV(std::string dir, std::string title, const std::vector<std::string>& table)
{
	if (title == "") return;

	// .csvを保存するディレクトリ
	std::string csvPath = params.outPath + "\\" + mylib::imgio::getImgName(params.imgName) + "\\" + dir;

	// ファイルが存在(すでにmakeCSV)していればCSVにtableの情報を書き込む
	if (std::filesystem::exists(csvPath) == true) {

		// .csvのテーブルを作成
		std::string _table = "";
		for (u32 i = 0; i < table.size(); i++) {
			_table += table[i] + ",";
		}

		// .csvに書き込み
		std::ofstream outputfile(csvPath + "\\" + title + ".csv", std::ios_base::app);
		outputfile << _table << std::endl;
	}
}

// ===================================================================
// cv::Mat -> std::vector
matrix mylib::img::MtoV(const cv::Mat& cvmat) {
	matrix vecs(cvmat.rows, std::vector<double>(cvmat.cols));
	//Mat -> vector
	for (int i = 0; i < cvmat.rows; i++) {
		for (int j = 0; j < cvmat.cols; j++) {
			auto x = double(cvmat.at<unsigned char>(i, j));
			vecs[i][j] = x;
		}
	}
	return vecs;
}

// std::vector -> cv::Mat
cv::Mat mylib::img::VtoM(const matrix& vecs) {
	// vecの大きさを取得
	int vec_rows = static_cast<int>(vecs.at(0).size());	// 行
	int vec_cols = static_cast<int>(vecs.size());	// 列
	cv::Mat cvmat(cv::Size(vec_cols, vec_rows), CV_64F);
	cvmat.convertTo(cvmat, CV_64F);
	for (int i = 0; i < cvmat.rows; i++) {
		for (int j = 0; j < cvmat.cols; j++) {
			cvmat.at<double>(i, j) = vecs[i][j];
		}
	}

	cvmat.convertTo(cvmat, CV_8U);
	return cvmat;
}

// resize ( size, size )
matrix mylib::img::resize(const matrix& img, const double size) 
{
	cv::Mat target = mylib::img::VtoM(img);
	cv::resize(target, target, cv::Size(), size, size);
	return mylib::img::MtoV(target);
}

// 複数の画像に同じ分散のガウスノイズを付与
mat3D mylib::img::SVcorruption(matrix& src)
{
	u32 rows = static_cast<u32>(src.at(0).size());	// 行
	u32 cols = static_cast<u32>(src.size());	// 列

	// ノイズ画像を格納するvectorを作成
	mat3D noise(params.sv_enumerate, matrix(cols, vec(rows)));

	// ガウス乱数を生成してvectorに格納する
	std::random_device seed_gen;	// シードの初期化
	std::default_random_engine gen(seed_gen());	// 乱数生成
	std::normal_distribution<double> dist(params.sv_m, params.sv_stddev);	// ガウス分布

	// クリッピング
	for(u32 k=0; k < params.sv_enumerate; k++) {
		for (u32 i = 0; i < rows; i++) {
			for (u32 j = 0; j < cols; j++) {
				const double tmp = dist(gen);
				noise[k][i][j] = src[i][j] + tmp;
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

// 複数の画像に異なる分散のガウスノイズを付与
mat3D mylib::img::DVcorruption(matrix& src)
{
	u32 rows = static_cast<u32>(src.at(0).size());	// 行
	u32 cols = static_cast<u32>(src.size());	// 列

	if (!params.dv_isRandom)
	{
		mat3D noise(params.dv_stddev.size(), std::vector<std::vector<double>>(cols, std::vector<double>(rows, 0.0)));
		for (u32 k = 0; k < params.dv_stddev.size(); k++) {
			// ガウス乱数を生成してvectorに格納する
			std::random_device seed_gen;	// シードの初期化
			std::default_random_engine gen(seed_gen());	// 乱数生成
			std::normal_distribution<double> dist(params.dv_m, params.dv_stddev[k]);	// ガウス分布

			for (u32 i = 0; i < rows; i++) {
				for (u32 j = 0; j < cols; j++) {
					const double tmp = dist(gen);
					noise[k][i][j] = src[i][j] + tmp;
				}
			}
		}
		return noise;
	}
	else
	{
		mat3D noise(params.dv_enumerate, std::vector<std::vector<double>>(cols, std::vector<double>(rows, 0.0)));
		for (u32 k = 0; k < params.dv_enumerate; k++) {
			// ガウス乱数を生成してvectorに格納する
			std::random_device seed_gen;	// シードの初期化
			std::default_random_engine gen(seed_gen());	// 乱数生成
			std::normal_distribution<double> dist(params.dv_m, mylib::img::getRand(params.dv_minStddev, params.dv_maxStddev));	// ガウス分布

			for (u32 i = 0; i < rows; i++) {
				for (u32 j = 0; j < cols; j++) {
					const double tmp = dist(gen);
					noise[k][i][j] = src[i][j] + tmp;
				}
			}
		}
		return noise;
	}
}

// PSNRを算出
double mylib::img::psnr(const matrix& img1, const matrix& img2) 
{
	u32 rows = static_cast<u32>(img1.at(0).size());	// 行
	u32 cols = static_cast<u32>(img1.size());	// 列
	u32 n = rows * cols;

	// 行列のサイズが違う場合に例外送出
	if (rows != img2.at(0).size() || cols != img2.size()) {
		std::cout << "vector size error : [" << rows << "," << cols << "] to["
			<< img2.at(0).size() << "," << img2.size() << std::endl;
		return 0;
	}

	double sum = 0.0;
	for (u32 i = 0; i < rows; i++) {
		for (u32 j = 0; j < cols; j++) {
			sum += std::pow(img1[i][j] - img2[i][j], 2);
		}
	}

	return log10(pow(255, 2) / (sum/n)) * 10.0;
}

// SSIMを算出
double mylib::img::ssim(const matrix& _img1, const matrix& _img2)
{

	// 行列のサイズが違う場合に出力
	if (_img1.at(0).size() != _img2.at(0).size() || _img1.size() != _img2.size()) {
		std::cout << "vector size error : [" << _img1.at(0).size() << "," << _img1.size() << "] to["
			<< _img2.at(0).size() << "," << _img2.size() << std::endl;
		return 0;
	}

	const double C1 = 6.5025, C2 = 58.5225;


	cv::Mat img1 = mylib::img::VtoM(_img1);
	cv::Mat img2 = mylib::img::VtoM(_img2);

	cv::Mat img1f, img2f;
	img1.convertTo(img1f, CV_32F);
	img2.convertTo(img2f, CV_32F);

	cv::Mat img1_sq = img1f.mul(img1f);
	cv::Mat img2_sq = img2f.mul(img2f);
	cv::Mat img1_img2 = img1f.mul(img2f);

	cv::Mat mu1, mu2;
	cv::GaussianBlur(img1f, mu1, cv::Size(11, 11), 1.5);
	cv::GaussianBlur(img2f, mu2, cv::Size(11, 11), 1.5);

	cv::Mat mu1_sq = mu1.mul(mu1);
	cv::Mat mu2_sq = mu2.mul(mu2);
	cv::Mat mu1_mu2 = mu1.mul(mu2);

	cv::Mat sigma1_sq, sigma2_sq, sigma12;

	GaussianBlur(img1_sq, sigma1_sq, cv::Size(11, 11), 1.5);
	sigma1_sq -= mu1_sq;

	GaussianBlur(img2_sq, sigma2_sq, cv::Size(11, 11), 1.5);
	sigma2_sq -= mu2_sq;

	GaussianBlur(img1_img2, sigma12, cv::Size(11, 11), 1.5);
	sigma12 -= mu1_mu2;

	cv::Mat ssim_map;
	cv::divide((2 * mu1_mu2 + C1) * (2 * sigma12 + C2), (mu1_sq + mu2_sq + C1) * (sigma1_sq + sigma2_sq + C2), ssim_map);

	cv::Scalar ssim_scalar = mean(ssim_map);

	return static_cast<double>(ssim_scalar[0]);
}

// 画像の中心化
matrix mylib::img::averaged(const mat3D& imgs)
{
	u32 enumerate = static_cast<u32>(imgs.size());	// 枚数の取得
	u32 rows = static_cast<u32>(imgs.at(0).at(0).size());
	u32 cols = static_cast<u32>(imgs.at(0).size());
	u32 n = rows * cols;

	matrix averaged(cols, std::vector<double>(rows, 0.0));

	for (u32 k = 0; k < enumerate; k++) {
		for (u32 i = 0; i < rows; i++) {
			for (u32 j = 0; j < cols; j++) {
				averaged[i][j] += imgs[k][i][j] / enumerate;
			}
		}
	}
	return averaged;
}

// 画像の平均値を算出
double mylib::img::calcMean(const matrix& img)
{
	u32 rows = static_cast<u32>(img.at(0).size());	// 行
	u32 cols = static_cast<u32>(img.size());	// 列
	u32 n = rows * cols;

	double imgEx = 0.0;
	for (u32 i = 0; i < rows; i++) {
		for (u32 j = 0; j < cols; j++) {
			imgEx += img[i][j];
		}
	}
	return imgEx / n;
}

// 一様乱数を生成
double mylib::img::getRand(double minVal, double maxVal) {
	// 乱数生成器
	static std::mt19937_64 mt64(0);

	// [min_val, max_val] の一様分布整数 (int) の分布生成器
	std::uniform_int_distribution<u32> get_rand_uni_int((double)minVal, (double)maxVal);

	// 乱数を生成
	return static_cast<double>(get_rand_uni_int(mt64));
}

// ================================================================
// 波形データの生成
vec mylib::audio::genWave()
{
	// buf = {0, 0, ... , 0}
	vec buf(params.dataSize);

	for (uint16_t t = 0; t < params.sampleRate; t++) {
		for (uint16_t i = 0; i < params.freq.size(); i++) {
			buf[t] += params.amplitude * std::sin(2 * M_PI * params.freq[i] * t / params.sampleRate) / params.freq.size();
		}
	}
	/*
	for (u32 i = 0; i < params.bufSize; i++) {
		std::cout << buf[i] << std::endl;
	}
	*/
	return buf;
}

// 1次元データのプロット
void mylib::audio::wavePlot(std::string title, const vec& wave)
{
	std::string path = params.wavePlotPath + "\\" + title + ".csv";
	if (params.wavePlotPath == "") return; 

	// CSVファイルの作成
	std::ofstream outputfile(path, std::ios_base::trunc);
	outputfile << "sample" << "," << "bit" << "," << title << std::endl;
	for (u32 i = 0; i < params.plotSize; i++) {
		outputfile << i << "," << wave[i] << std::endl;
	}
	
	// waveplot.pyを実行してグラフをプロット
	// system("D:\\yasuda-lab\\B4\\gmrf\\gmrf\\waveplot.py");
}

// 1次元データに同じ分散のノイズを付与
matrix mylib::audio::SVCorruptWave(const vec& src)
{
	matrix noise(params.vec_enumerate, std::vector<double>(params.bufSize));

	// ガウス乱数を生成してvectorに格納する
	std::random_device seed_gen;	// シードの初期化
	std::default_random_engine gen(seed_gen());	// 乱数生成
	std::normal_distribution<double> dist(params.vec_m, params.vec_stddev);	// ガウス分布

	for (u32 k = 0; k < params.vec_enumerate; k++) {
		for(u32 i = 0; i < src.size(); i++) {
			const double tmp = dist(gen);
			noise[k][i] = src[i] + tmp;
			if (noise[k][i] > params.amplitude) {
				noise[k][i] = 1;
			}
			if (noise[k][i] < -1.0 * params.amplitude) {
				noise[k][i] = -1;
			}
		}
	}
	return noise;
}

// 1次元データに異なる分散のノイズを付与
matrix mylib::audio::DVCorruptWave(const vec& src)
{
	if (params.dvod_isRandom != true)
	{
		matrix noise(params.dvod_stddev.size(), std::vector<double>(src.size()));
		for (u32 k = 0; k < params.dvod_stddev.size(); k++)
		{
			// ガウス乱数を生成してvectorに格納する
			std::random_device seed_gen;	// シードの初期化
			std::default_random_engine gen(seed_gen());	// 乱数生成
			std::normal_distribution<double> dist(params.vec_m, params.dvod_stddev[k]);	// ガウス分布

			for (u32 i = 0; i < src.size(); i++) {
				const double tmp = dist(gen);
				noise[k][i] = src[i] + tmp;
				if (noise[k][i] > params.amplitude) {
					noise[k][i] = params.amplitude;
				}
				if (noise[k][i] < -1.0 * params.amplitude) {
					noise[k][i] = -1.0 * params.amplitude;
				}
			}
		}
		return noise;
	}
	else
	{
		matrix noise(params.dvod_enumerate, std::vector<double>(src.size()));
		for (u32 k = 0; k < params.dvod_enumerate; k++)
		{
			// ガウス乱数を生成してvectorに格納する
			std::random_device seed_gen;	// シードの初期化
			std::default_random_engine gen(seed_gen());	// 乱数生成
			std::normal_distribution<double> dist(params.dvod_m, static_cast<double>(mylib::img::getRand(params.dvod_minStddev*100, params.dvod_maxStddev*100) / 100));	// ガウス分布

			for (u32 i = 0; i < src.size(); i++) {
				const double tmp = dist(gen);
				noise[k][i] = src[i] + tmp;
				if (noise[k][i] > params.amplitude) {
					noise[k][i] = params.amplitude;
				}
				if (noise[k][i] < -1.0 * params.amplitude) {
					noise[k][i] = -1.0 * params.amplitude;
				}
			}
		}
		return noise;
	}
}

// 1次元データのPSNR
double mylib::audio::psnr(const vec& data1, const vec& data2)
{
	u32 dataSize = static_cast<u32>(data1.size());

	// ベクトルサイズの例外処理
	if (dataSize != data2.size()) {
		std::cout << "vector size error : " << dataSize << "to" << data2.size() << std::endl;
	}

	double sum = 0.0;
	for (u32 i = 0; i < dataSize; i++) {
		sum += std::pow(data1[i] - data2[i], 2);
	}

	return log10(std::pow(params.amplitude * 2, 2) / (sum / dataSize)) * 10.0;
}

// 新しいCSVを作成
void mylib::audio::makeCSV(std::string dir, std::string title, const std::vector<std::string>& table)
{
	if (title == "") return;

	// .csvを保存するディレクトリの作成(初回のみ)
	std::string csvPath = params.wavePlotPath + "\\" + dir;
	std::filesystem::create_directories(csvPath);

	// .csvのテーブルを作成
	std::string _table = "";
	for (const auto& s : table) _table += s + ",";

	// .csvに書き込み
	std::ofstream outputfile(csvPath + "\\" + title + ".csv", std::ios_base::trunc);
	outputfile << _table << std::endl;
}

// CSVを保存
void mylib::audio::saveCSV(std::string dir, std::string title, const std::vector<std::string>& table)
{
	if (title == "") return;

	// .csvを保存するディレクトリ
	std::string csvPath = params.wavePlotPath + "\\" + dir;

	// ファイルが存在(すでにmakeCSV)していればCSVにtableの情報を書き込む
	if (std::filesystem::exists(csvPath) == true) {

		// .csvのテーブルを作成
		std::string _table = "";
		for (const auto& s : table) _table += s + ",";

		// .csvに書き込み
		std::ofstream outputfile(csvPath + "\\" + title + ".csv", std::ios_base::app);
		outputfile << _table << std::endl;
	}
}

// ================================================================
// accessor
bool mylib::utillity::getImgVisualize()
{
	return params.imgVisualize;
}
bool mylib::utillity::getInfoVisualize()
{
	return params.infoVisualize;
}
bool mylib::utillity::getSaveImg()
{
	return params.saveImg;
}
bool mylib::utillity::getSaveCsv()
{
	return params.saveCsv;
}
bool mylib::utillity::getSavePlot()
{
	return params.savePlot;
}
i32 mylib::utillity::getIter()
{
	return params.iter;
}
std::string mylib::utillity::getCsvName()
{
	return params.csvName;
}

// 進捗の表示
void mylib::utillity::progressBar()
{
	params.progress++;

	// 進捗%を計算
	u32 prog = static_cast<u32>((params.progress * 1.0 / params.iter) * 100.0);

	// 進捗の表示
	printf("\rprogress... %d%%", prog);
}

// より詳細な桁数の文字列変換
std::string mylib::utillity::float_to_string(double f, int digits)
{
	std::ostringstream oss;
	oss << std::setprecision(digits) << std::setiosflags(std::ios::fixed) << f;

	return oss.str();
}

// 数値の整合性チェック
bool mylib::utillity::is_integrity(double num) {
	if (!std::isfinite((float)num) || num < 0) return false;
	return true;
}