#pragma once

#include "include.hpp"

namespace mylib
{
	// 画像操作
	namespace img
	{
		// 関数
		matrix MtoV(const cv::Mat& img);
		cv::Mat VtoM(const matrix& img);
		mat3D SVcorruption(matrix& img);
		mat3D DVcorruption(matrix& src);
		double psnr(const matrix& img1, const matrix& img2);
		double ssim(const matrix& img1, const matrix& img2);
		double getRand(double minVal, double maxVal);
		matrix resize(const matrix& img, const double size);
		matrix averaged(const mat3D& imgs);
		double calcMean(const matrix& img);
	}

	// 画像入出力
	namespace imgio
	{
		// 画像の種類
		enum ImgIndex
		{
			BOAT,
			Clock,
			girl,
			Jellybeans,
			LENNA,
			Lighthouse,
			Man,
			barbara512,
			APC,
			Stream_and_bridge,
			Tank,
			House,
			Moon_surface,
			Tree
		};

		std::string getImgName(ImgIndex imgName);
		matrix loadImg();
		void saveImg(std::string dir, std::string title, const matrix& img);
		void makeCSV(std::string dir, std::string title, const std::vector<std::string>& table);
		void saveCSV(std::string dir, std::string title, const std::vector<std::string>& table);
	}

	// 1次元データに対する操作
	namespace audio
	{
		vec genWave();
		void wavePlot(std::string title, const vec& wave);
		matrix SVCorruptWave(const vec& src);
		matrix DVCorruptWave(const vec& src);
		double psnr(const vec& data1, const vec& data2);
		void makeCSV(std::string dir, std::string title, const std::vector<std::string>& table);
		void saveCSV(std::string dir, std::string title, const std::vector<std::string>& table);
	}

	// accessor
	namespace utillity
	{
		bool getImgVisualize();
		bool getInfoVisualize();
		bool getSaveImg();
		bool getSaveCsv();
		bool getSavePlot();
		i32 getIter();
		std::string getCsvName();
		void progressBar();
		std::string float_to_string(double f, int digits);
		bool is_integrity(double num);
	}
}