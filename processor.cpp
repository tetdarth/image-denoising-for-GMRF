#include "processor.hpp"

// ============================== SVGMRF ====================================
void SVGMRF::processor::process()
{
	matrix originImg = mylib::imgio::loadImg();
	mat3D corrutedImg = mylib::img::SVcorruption(originImg);

	// 実行時間計測 start //
	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	// GMRFの呼び出し
	IVGMRF::gmrf<double> GMRF;
	matrix denoisedImg = GMRF.processBlock(corrutedImg);

	end = std::chrono::system_clock::now();
	double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
	printf("time %lf[ms]\n", time);
	// 実行時間計測 end //

	std::cout << "epoch : " << GMRF.getEpoch() << std::endl;
	std::cout << "averaged PSNR : " << mylib::img::psnr(GMRF.getAvgImg(), originImg) << std::endl;
	std::cout << "denoised PSNR : " << mylib::img::psnr(denoisedImg, originImg) << std::endl;

	// mylib::imgio::saveImg("SVGMRF", "Avg", GMRF.getAvgImg());
	// mylib::imgio::saveImg("SVGMRF", "denoised", denoisedImg);

	cv::imshow("originImg", mylib::img::VtoM(originImg));
	cv::imshow("corruptedImg", mylib::img::VtoM(corrutedImg.at(0)));
	cv::imshow("averagedImg", mylib::img::VtoM(GMRF.getAvgImg()));
	cv::imshow("denoisedImg",mylib::img::VtoM(denoisedImg));
	cv::waitKey(0);
}

// ============================== GMRFUS ====================================
void GMRFUpScale::processor::process()
{
	matrix originImg = mylib::imgio::loadImg();
	matrix resizeImg = mylib::img::resize(originImg, 0.5);

	// 実行時間計測 start //
	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	IVGMRF::upScale<double> gmrfus;
	matrix upScalled = gmrfus.processBlock(resizeImg);

	end = std::chrono::system_clock::now();
	double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
	printf("time %lf[ms]\n", time);
	// 実行時間計測 end //

	std::cout << "cv - psnr :" << mylib::img::psnr(mylib::img::resize(resizeImg, 2.0), originImg) << std::endl;
	std::cout << "GMRFpsnr :" << mylib::img::psnr(upScalled, originImg) << std::endl;

	mylib::imgio::saveImg("GMRFUS", "ocvUS", (mylib::img::resize(resizeImg, 2.0)));
	mylib::imgio::saveImg("GMRFUS", "GMRFUS", upScalled);

	// cv::imshow("origin", mylib::img::VtoM(originImg));
	cv::imshow("opencv", mylib::img::VtoM(mylib::img::resize(resizeImg, 2.0)));
	cv::imshow("upScale",mylib::img::VtoM(upScalled));
	cv::waitKey(0);
}

// ============================== DVGMRF ====================================
void DVGMRF::processor::process()
{
	matrix originImg = mylib::imgio::loadImg();
	mat3D corrutedImg = mylib::img::DVcorruption(originImg);

	// 実行時間計測 start //
	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	// GMRFの呼び出し
	DVGMRF::gmrf<double> GMRF;
	matrix denoisedImg = GMRF.processBlock(corrutedImg);

	end = std::chrono::system_clock::now();
	double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
	printf("time %lf[ms]\n", time);
	// 実行時間計測 end //

	std::cout << "lambda : " << GMRF.getLambda() << std::endl;
	std::cout << "alpha : " << GMRF.getAlpha() << std::endl;

	std::cout << "epoch : " << GMRF.getEpoch() << std::endl;
	std::cout << "averaged PSNR : " << mylib::img::psnr(GMRF.getAvgImg(), originImg) << std::endl;
	std::cout << "denoised PSNR : " << mylib::img::psnr(denoisedImg, originImg) << std::endl;
	std::cout << "denoised SSIM : " << mylib::img::ssim(denoisedImg, originImg) << std::endl;
	std::cout << "enumerate : " << corrutedImg.size() << std::endl;

	mylib::imgio::saveImg("DVGMRF", "Avg", GMRF.getAvgImg());
	mylib::imgio::saveImg("DVGMRF", "denoised", denoisedImg);

	cv::imshow("originImg", mylib::img::VtoM(originImg));
	cv::imshow("corruptedImg", mylib::img::VtoM(corrutedImg.at(0)));
	cv::imshow("averagedImg", mylib::img::VtoM(mylib::img::averaged(corrutedImg)));
	cv::imshow("denoisedImg", mylib::img::VtoM(denoisedImg));
	cv::waitKey(0);
}

// ============================== DVSVCOMP ==================================
void DVSVCOMP::processor::process()
{
	std::string dirName = "SVDVCOMP";

	// CSVの作成
	if (mylib::utillity::getSaveCsv() == true) {
		std::vector<std::string> table = { 
			"avgPSNR", 
			"avgSSIM", 
			"svPSNR", 
			"svSSIM", 
			"dvPSNR", 
			"dvSSIM", 
			"svTime", 
			"dvTime",
			"",
			// 平均算出用テーブル
			"=AVERAGE(A:A)",
			"=AVERAGE(B:B)",
			"=AVERAGE(C:C)",
			"=AVERAGE(D:D)",
			"=AVERAGE(E:E)",
			"=AVERAGE(F:F)",
			"=AVERAGE(G:G)",
			"=AVERAGE(H:H)"
		};
		mylib::imgio::makeCSV(dirName, mylib::utillity::getCsvName(), table);
	}

	// 並列で処理を実行
#pragma omp parallel for num_threads(THREADS)
	for (i32 i = 0; i < mylib::utillity::getIter(); i++) {
		// 画像の呼び出しと劣化
		matrix originImg = mylib::imgio::loadImg();
		mat3D corrutedImg = mylib::img::DVcorruption(originImg);

		// 実行時間計測 start //
		std::chrono::system_clock::time_point start, end;
		start = std::chrono::system_clock::now();

		// GMRFの呼び出し
		IVGMRF::gmrf<double> ivgmrf;
		matrix SVdenoising = ivgmrf.processBlock(corrutedImg);

		end = std::chrono::system_clock::now();
		double svTime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
		// 実行時間計測 end //


		// 実行時間計測 start //
		start = std::chrono::system_clock::now();

		// GMRFの呼び出し
		DVGMRF::gmrf<double> dvgmrf;
		matrix DVdenoising = dvgmrf.processBlock(corrutedImg);

		end = std::chrono::system_clock::now();
		double dvTime = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
		// 実行時間計測 end //

		// info
		if (mylib::utillity::getInfoVisualize() == true) {
			std::cout << "averaged PSNR : " << mylib::img::psnr(dvgmrf.getAvgImg(), originImg) << std::endl;
			std::cout << "averaged SSIM : " << mylib::img::ssim(dvgmrf.getAvgImg(), originImg) << std::endl;
			std::cout << std::endl;
			std::cout << "============ SVGMRF ============= " << std::endl;
			std::cout << "processing time : " << svTime << std::endl;
			std::cout << "epoch : " << ivgmrf.getEpoch() << std::endl;
			std::cout << "lambda : " << ivgmrf.getLambda() << std::endl;
			std::cout << "alpha : " << ivgmrf.getAlpha() << std::endl;
			std::cout << "denoised PSNR : " << mylib::img::psnr(SVdenoising, originImg) << std::endl;
			std::cout << "denoised SSIM : " << mylib::img::ssim(SVdenoising, originImg) << std::endl;
			std::cout << std::endl;
			std::cout << "============ DVGMRF ============= " << std::endl;
			std::cout << "processing time : " << dvTime << std::endl;
			std::cout << "epoch : " << dvgmrf.getEpoch() << std::endl;
			std::cout << "lambda : " << dvgmrf.getLambda() << std::endl;
			std::cout << "alpha : " << dvgmrf.getAlpha() << std::endl;
			std::cout << "denoised PSNR : " << mylib::img::psnr(DVdenoising, originImg) << std::endl;
			std::cout << "denoised SSIM : " << mylib::img::ssim(DVdenoising, originImg) << std::endl;
		}

		// CSVに書き込み
		if (mylib::utillity::getSaveCsv() == true) {
			std::vector<std::string> info_table = { std::to_string(mylib::img::psnr(ivgmrf.getAvgImg(), originImg)),
													std::to_string(mylib::img::ssim(ivgmrf.getAvgImg(), originImg)),
													std::to_string(mylib::img::psnr(SVdenoising, originImg)),
													std::to_string(mylib::img::ssim(SVdenoising, originImg)),
													std::to_string(mylib::img::psnr(DVdenoising, originImg)),
													std::to_string(mylib::img::ssim(DVdenoising, originImg)),
													std::to_string(svTime),
													std::to_string(dvTime)
			};
			mylib::imgio::saveCSV(dirName, mylib::utillity::getCsvName(), info_table);
		}
		
		// 画像の保存
		if (mylib::utillity::getSaveImg() == true)
		{
			// mylib::imgio::saveImg(dirName, "Original", originImg);
			mylib::imgio::saveImg(dirName, "Avg", ivgmrf.getAvgImg());
			mylib::imgio::saveImg(dirName, "SVGMRF", SVdenoising);
			mylib::imgio::saveImg(dirName, "DVGMRF", DVdenoising);
		}

		// 画像の表示
		if (mylib::utillity::getImgVisualize() == true) {
			cv::imshow("originImg", mylib::img::VtoM(originImg));
			cv::imshow("averagedImg", mylib::img::VtoM(ivgmrf.getAvgImg()));
			cv::imshow("denoised for SVGMRF", mylib::img::VtoM(SVdenoising));
			cv::imshow("denoised for DVGMRF", mylib::img::VtoM(DVdenoising));
			cv::waitKey(0);
		}

		// 進捗の表示
		mylib::utillity::progressBar();
	}
}

// ============================== ODGMRF ====================================
void SVODGMRF::processor::process()
{
	vec originData = mylib::audio::genWave();
	mylib::audio::wavePlot("origindata", originData);
	matrix corruptedData = mylib::audio::SVCorruptWave(originData);
	mylib::audio::wavePlot("corrupted", corruptedData[0]);

	// 実行時間計測 start //
	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	IVGMRF::odgmrf<double> GMRF;
	vec denoisedData = GMRF.processBlock(corruptedData);

	end = std::chrono::system_clock::now();
	double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
	printf("time %lf[ms]\n", time);
	// 実行時間計測 end //

	mylib::audio::wavePlot("denoised", denoisedData);

	std::cout << "epoch : " << GMRF.getEpoch() << std::endl;
	std::cout << "averaged PSNR : " << mylib::audio::psnr(GMRF.getAvgData(), originData) << std::endl;
	std::cout << "denoised PSNR : " << mylib::audio::psnr(denoisedData, originData) << std::endl;
}

// ============================== DVODGMRF ==================================
void DVODGMRF::processor::process()
{
	vec originData = mylib::audio::genWave();
	mylib::audio::wavePlot("origindata", originData);
	matrix corruptedData = mylib::audio::DVCorruptWave(originData);
	/*
	mylib::audio::wavePlot("corrupt[0]", corruptedData[0]);
	mylib::audio::wavePlot("corrupt[1]", corruptedData[1]);
	mylib::audio::wavePlot("corrupt[2]", corruptedData[2]);
	*/

	// 実行時間計測 start //
	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	DVGMRF::odgmrf<double> GMRF;
	vec denoisedData = GMRF.processBlock(corruptedData);

	end = std::chrono::system_clock::now();
	double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
	printf("time %lf[ms]\n", time);
	// 実行時間計測 end //
	
	vec averagedData = GMRF.getAvgData();
	mylib::audio::wavePlot("averaged", averagedData);
	mylib::audio::wavePlot("denoised", denoisedData);

	std::cout << "epoch : " << GMRF.getEpoch() << std::endl;
	std::cout << "averaged PSNR : " << mylib::audio::psnr(averagedData, originData) << std::endl;
	std::cout << "denoised PSNR : " << mylib::audio::psnr(denoisedData, originData) << std::endl;
}

// ============================== DVSVODCOMP ================================
void DVSVODCOMP::processor::process()
{
	std::string dirName = "SVDVODCOMP";

	// CSVのテーブルを作成
	if (mylib::utillity::getSaveCsv() == true) {
		std::vector<std::string> table = {
			"avgPSNR",
			"svPSNR",
			"dvPSNR",
			"svTime",
			"dvTime",
			"",
			// 平均算出用テーブル
			"=AVERAGE(A:A)",
			"=AVERAGE(B:B)",
			"=AVERAGE(C:C)",
			"=AVERAGE(D:D)",
			"=AVERAGE(E:E)"
		};
		mylib::audio::makeCSV(dirName, mylib::utillity::getCsvName(), table);
	}

	// 並列で処理を実行
#pragma omp parallel for num_threads(THREADS)
	for (i32 i = 0; i < mylib::utillity::getIter(); i++) {
		// 波形の生成と劣化
		vec originData = mylib::audio::genWave();
		mylib::audio::wavePlot("origindata", originData);
		matrix corruptedData = mylib::audio::DVCorruptWave(originData);

		/*
		mylib::audio::wavePlot("corrupt[0]", corruptedData[0]);
		mylib::audio::wavePlot("corrupt[1]", corruptedData[1]);
		mylib::audio::wavePlot("corrupt[2]", corruptedData[2]);
		*/

		std::chrono::system_clock::time_point start, end;

		// 実行時間計測 start //
		start = std::chrono::system_clock::now();

		IVGMRF::odgmrf<double> SVGMRF;
		vec SVdenoisedData = SVGMRF.processBlock(corruptedData);

		end = std::chrono::system_clock::now();
		double sv_time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
		// 実行時間計測 end //

		// 実行時間計測 start //
		start = std::chrono::system_clock::now();

		DVGMRF::odgmrf<double> DVGMRF;
		vec DVdenoisedData = DVGMRF.processBlock(corruptedData);

		end = std::chrono::system_clock::now();
		double dv_time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
		// 実行時間計測 end //

		// CSVに書き込み
		if (mylib::utillity::getSaveCsv() == true) {
			std::vector<std::string> table = {
				std::to_string(mylib::audio::psnr(DVGMRF.getAvgData(), originData)),
				std::to_string(mylib::audio::psnr(SVdenoisedData, originData)),
				std::to_string(mylib::audio::psnr(DVdenoisedData, originData)),
				std::to_string(sv_time),
				std::to_string(dv_time)
			};
			mylib::audio::saveCSV(dirName, mylib::utillity::getCsvName(), table);
		}
		
		// plot
		if (mylib::utillity::getSavePlot() == true) {
			mylib::audio::wavePlot("corrupted0", corruptedData.at(0));
			mylib::audio::wavePlot("corrupted1", corruptedData.at(1));
			mylib::audio::wavePlot("corrupted2", corruptedData.at(2));
			mylib::audio::wavePlot("averaged", DVGMRF.getAvgData());
			mylib::audio::wavePlot("denoised_for_SVGMRF", SVdenoisedData);
			mylib::audio::wavePlot("denoised_for_DVGMRF", DVdenoisedData);
		}

		// info
		if (mylib::utillity::getInfoVisualize() == true) {
			std::cout << "averaged PSNR : " << mylib::audio::psnr(DVGMRF.getAvgData(), originData) << std::endl;
			std::cout << std::endl;
			std::cout << "========== SVGMRF =========" << std::endl;
			std::cout << "execute time : " << sv_time << "[ms]" << std::endl;
			std::cout << "epoch : " << SVGMRF.getEpoch() << std::endl;
			std::cout << "alpha for SVGMRF : " << SVGMRF.getAlpha() << std::endl;
			std::cout << "lambda for SVGMRF : " << SVGMRF.getLambda() << std::endl;
			std::cout << "denoised PSNR for SVGMRF : " << mylib::audio::psnr(SVdenoisedData, originData) << std::endl;
			std::cout << std::endl;
			std::cout << "========== DVGMRF =========" << std::endl;
			std::cout << "execute time : " << dv_time << "[ms]" << std::endl;
			std::cout << "epoch : " << DVGMRF.getEpoch() << std::endl;
			std::cout << "alpha for DVGMRF : " << DVGMRF.getAlpha() << std::endl;
			std::cout << "lambda for DVGMRF : " << DVGMRF.getLambda() << std::endl;
			std::cout << "denoised PSNR for DVGMRF : " << mylib::audio::psnr(DVdenoisedData, originData) << std::endl;
		}

		// 進捗の表示
		mylib::utillity::progressBar();
	}
}

// ============================== IVHGMRF ===================================
void IVHGMRF::processor::process() 
{
	matrix originImg = mylib::imgio::loadImg();
	mat3D corrutedImg = mylib::img::DVcorruption(originImg);

	// 実行時間計測 start //
	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	// GMRFの呼び出し
	HGMRF::ivhgmrf<double> GMRF;
	matrix denoisedImg = GMRF.denoising(corrutedImg);

	end = std::chrono::system_clock::now();
	double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
	printf("time %lf[ms]\n", time);
	// 実行時間計測 end //

	std::cout << "epoch : " << GMRF.get_epoch() << std::endl;
	std::cout << "averaged PSNR : " << mylib::img::psnr(GMRF.get_avg_img(), originImg) << std::endl;
	std::cout << "denoised PSNR : " << mylib::img::psnr(denoisedImg, originImg) << std::endl;
	std::cout << "alpha : " << GMRF.get_alpha() << std::endl;
	std::cout << "lambda ; " << GMRF.get_lambda() << std::endl;
	std::cout << "gammma2 : " << GMRF.get_gammma2() << std::endl;
	std::cout << "sigma2 : " << GMRF.get_sigma2() << std::endl;

	cv::imshow("originImg", mylib::img::VtoM(originImg));
	cv::imshow("corruptedImg", mylib::img::VtoM(corrutedImg.at(0)));
	cv::imshow("averagedImg", mylib::img::VtoM(GMRF.get_avg_img()));
	cv::imshow("denoisedImg", mylib::img::VtoM(denoisedImg));
	cv::waitKey(0);
}


void HGMRFCOMP::processor::process()
{
	matrix originImg = mylib::imgio::loadImg();
	mat3D corrutedImg = mylib::img::SVcorruption(originImg);

	std::string dirName = "HGMRFCOMP";

	// 実行時間計測 start //
	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	// GMRFの呼び出し
	IVGMRF::gmrf<double> GMRF;
	matrix denoised_gmrf = GMRF.processBlock(corrutedImg);

	end = std::chrono::system_clock::now();
	double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
	printf("time %lf[ms]\n", time);
	// 実行時間計測 end //

	std::cout << "epoch : " << GMRF.getEpoch() << std::endl;
	// std::cout << "averaged PSNR : " << mylib::img::psnr(GMRF.getAvgImg(), originImg) << std::endl;
	std::cout << "denoised PSNR : " << mylib::img::psnr(denoised_gmrf, originImg) << std::endl;
	std::cout << "denoised SSIM : " << mylib::img::ssim(denoised_gmrf, originImg) << std::endl;
	std::cout << "alpha : " << GMRF.getAlpha() << std::endl;
	std::cout << "lammbda : " << GMRF.getLambda() << std::endl;
	std::cout << "sigma2 : " << GMRF.getSigma2() << std::endl;
	std::cout << "===========================================================" << std::endl;

	// mylib::imgio::saveImg("SVGMRF", "Avg", GMRF.getAvgImg());
	// mylib::imgio::saveImg("SVGMRF", "denoised", denoisedImg);

	// 実行時間計測 start //
	start = std::chrono::system_clock::now();

	// HGMRFの呼び出し
	HGMRF::ivhgmrf<double> HGMRF;
	matrix denoised_hgmrf = HGMRF.denoising(corrutedImg);

	end = std::chrono::system_clock::now();
	time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
	printf("time %lf[ms]\n", time);
	// 実行時間計測 end //

	std::cout << "epoch : " << HGMRF.get_epoch() << std::endl;
	// std::cout << "averaged PSNR : " << mylib::img::psnr(HGMRF.get_avg_img(), originImg) << std::endl;
	std::cout << "denoised PSNR : " << mylib::img::psnr(denoised_hgmrf, originImg) << std::endl;
	std::cout << "denoised SSIM : " << mylib::img::ssim(denoised_hgmrf, originImg) << std::endl;
	std::cout << "alpha : " << HGMRF.get_alpha() << std::endl;
	std::cout << "lambda ; " << HGMRF.get_lambda() << std::endl;
	std::cout << "gammma2 : " << HGMRF.get_gammma2() << std::endl;
	std::cout << "sigma2 : " << HGMRF.get_sigma2() << std::endl;
	std::cout << "===========================================================" << std::endl;

	if (mylib::utillity::getSaveImg) {
		mylib::imgio::saveImg(dirName, "GMRF", denoised_gmrf);
		mylib::imgio::saveImg(dirName, "HGMRF", denoised_hgmrf);
	}

	cv::imshow("originImg", mylib::img::VtoM(originImg));
	cv::imshow("corruptedImg", mylib::img::VtoM(corrutedImg.at(0)));
	cv::imshow("averagedImg", mylib::img::VtoM(HGMRF.get_avg_img()));
	cv::imshow("denoised_gmrf", mylib::img::VtoM(denoised_gmrf));
	cv::imshow("denoised_hgmrf", mylib::img::VtoM(denoised_hgmrf));
	cv::waitKey(0);
}


void DVHGMRF::processor::process() {
	matrix originImg = mylib::imgio::loadImg();
	mat3D corrutedImg = mylib::img::DVcorruption(originImg);

	// 実行時間計測 start //
	std::chrono::system_clock::time_point start, end;
	start = std::chrono::system_clock::now();

	// GMRFの呼び出し
	HGMRF::dvhgmrf<double> GMRF;
	matrix denoisedImg = GMRF.denoising(corrutedImg);

	end = std::chrono::system_clock::now();
	double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);
	printf("time %lf[ms]\n", time);
	// 実行時間計測 end //

	std::cout << "epoch : " << GMRF.get_epoch() << std::endl;
	std::cout << "averaged PSNR : " << mylib::img::psnr(GMRF.get_avg_img(), originImg) << std::endl;
	std::cout << "denoised PSNR : " << mylib::img::psnr(denoisedImg, originImg) << std::endl;
	std::cout << "alpha : " << GMRF.get_alpha() << std::endl;
	std::cout << "lambda ; " << GMRF.get_lambda() << std::endl;
	std::cout << "gammma2 : " << GMRF.get_gammma2() << std::endl;
	vec vec_sigma2 = GMRF.get_vec_sigma2();
	for (int k = 0; k < vec_sigma2.size(); ++k) std::cout << "sigma2[" << k << "] : " << vec_sigma2[k] << std::endl;

	cv::imshow("originImg", mylib::img::VtoM(originImg));
	cv::imshow("corruptedImg", mylib::img::VtoM(corrutedImg.at(0)));
	cv::imshow("averagedImg", mylib::img::VtoM(GMRF.get_avg_img()));
	cv::imshow("denoisedImg", mylib::img::VtoM(denoisedImg));
	cv::waitKey(0);
}

void GMRFCOMP::processor::process() 
{
	std::string dirName = "GMRFCOMP";

	// CSVのテーブルを作成
	if (mylib::utillity::getSaveCsv()) {
		std::vector<std::string> table = {
			"IVGMRF_PSNR",
			"IVGMRF_SSIM",
			"DVGMRF_PSNR",
			"DVGMRF_SSIM",
			"IVHGMRF_PSNR",
			"IVHGMRF_SSIM",
			"DVHGMRF_PSNR",
			"DVHGMRF_SSIM",
			"",
			"",
			// 平均算出用テーブル
			"=AVERAGE(A:A)",
			"=AVERAGE(B:B)",
			"=AVERAGE(C:C)",
			"=AVERAGE(D:D)",
			"=AVERAGE(E:E)",
			"=AVERAGE(F:F)",
			"=AVERAGE(G:G)",
			"=AVERAGE(H:H)",
		};
		mylib::imgio::makeCSV(dirName, mylib::utillity::getCsvName(), table);
	}

#pragma omp parallel for num_threads(THREADS)
	for (i32 i = 0; i < mylib::utillity::getIter(); ++i) 
	{
		matrix originImg = mylib::imgio::loadImg();
		mat3D corrutedImg = mylib::img::DVcorruption(originImg);

		std::chrono::system_clock::time_point start, end;
		
		// ======== denoising for IVGMRF ========
		start = std::chrono::system_clock::now();

		IVGMRF::gmrf<double> ivgmrf;
		matrix denoised_ivgmrf = ivgmrf.processBlock(corrutedImg);

		end = std::chrono::system_clock::now();
		double time_ivgmrf = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);

		// ======== denoising for DVGMRF ========
		start = std::chrono::system_clock::now();

		DVGMRF::gmrf<double> dvgmrf;
		matrix denoised_dvgmrf = dvgmrf.processBlock(corrutedImg);

		end = std::chrono::system_clock::now();
		double time_dvgmrf = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);

		// ======== denoising for IVHGMRF ========
		start = std::chrono::system_clock::now();

		HGMRF::ivhgmrf<double> ivhgmrf;
		matrix denoised_ivhgmrf = ivhgmrf.denoising(corrutedImg);

		end = std::chrono::system_clock::now();
		double time_ivhgmrf = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);


		// ======= denoising for DVHGMRF ========
		start = std::chrono::system_clock::now();

		HGMRF::dvhgmrf<double> dvhgmrf;
		matrix denoised_dvhgmrf = dvhgmrf.denoising(corrutedImg);

		end = std::chrono::system_clock::now();
		double time_dvhgmrf = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0);

		// =====================================
		// 評価値の確保
		double ivgmrf_psnr = mylib::img::psnr(denoised_ivgmrf, originImg);
		double ivgmrf_ssim = mylib::img::ssim(denoised_ivgmrf, originImg);
		double dvgmrf_psnr = mylib::img::psnr(denoised_dvgmrf, originImg);
		double dvgmrf_ssim = mylib::img::ssim(denoised_dvgmrf, originImg);
		double ivhgmrf_psnr = mylib::img::psnr(denoised_ivhgmrf, originImg);
		double ivhgmrf_ssim = mylib::img::ssim(denoised_ivhgmrf, originImg);
		double dvhgmrf_psnr = mylib::img::psnr(denoised_dvhgmrf, originImg);
		double dvhgmrf_ssim = mylib::img::ssim(denoised_dvhgmrf, originImg);

		// 値不正合判定
		if (!mylib::utillity::is_integrity(ivgmrf_psnr) || !mylib::utillity::is_integrity(ivgmrf_ssim)
			|| !mylib::utillity::is_integrity(dvgmrf_psnr) || !mylib::utillity::is_integrity(dvgmrf_ssim)
			|| !mylib::utillity::is_integrity(ivhgmrf_psnr) || !mylib::utillity::is_integrity(ivhgmrf_ssim)
			|| !mylib::utillity::is_integrity(dvhgmrf_psnr) || !mylib::utillity::is_integrity(dvhgmrf_ssim))
		{ 
			mylib::utillity::progressBar();
			continue; 
		}

		// =====================================
		// 情報表示
		if (mylib::utillity::getInfoVisualize()) {
			std::cout << "IVGMRF_PSNR : " << ivgmrf_psnr << std::endl;
			std::cout << "DVGMRF_PSNR : " << dvgmrf_psnr << std::endl;
			std::cout << "IVHGMRF_PSNR : " << ivhgmrf_psnr << std::endl;
			std::cout << "DVHGMRF_PSNR : " << dvhgmrf_psnr << std::endl;
			std::cout << " " << std::endl;
			std::cout << "IVGMRF_SSIM : " << ivgmrf_ssim << std::endl;
			std::cout << "DVGMRF_SSIM : " << dvgmrf_ssim << std::endl;
			std::cout << "IVHGMRF_SSIM : " << ivhgmrf_ssim << std::endl;
			std::cout << "DVHGMRF_SSIM : " << dvhgmrf_ssim << std::endl;
			std::cout << " " << std::endl;
			std::cout << "IVGMRF_time : " << time_ivgmrf << std::endl;
			std::cout << "DVGMRF_time : " << time_dvgmrf << std::endl;
			std::cout << "IVHGMRF_time : " << time_ivhgmrf << std::endl;
			std::cout << "DVHGMRF_time : " << time_dvhgmrf << std::endl;
		}

		// CSVに書き込み
		if (mylib::utillity::getSaveCsv()) {
			std::vector<std::string> table = {
				std::to_string(ivgmrf_psnr),
				std::to_string(ivgmrf_ssim),
				std::to_string(dvgmrf_psnr),
				std::to_string(dvgmrf_ssim),
				std::to_string(ivhgmrf_psnr),
				std::to_string(ivhgmrf_ssim),
				std::to_string(dvhgmrf_psnr),
				std::to_string(dvhgmrf_ssim),
			};
			mylib::imgio::saveCSV(dirName, mylib::utillity::getCsvName(), table);
		}

		// 画像の保存
		if (mylib::utillity::getSaveImg()) {
			mylib::imgio::saveImg(dirName, "IVGMRF", denoised_ivgmrf);
			mylib::imgio::saveImg(dirName, "DVGMRF", denoised_dvgmrf);
			mylib::imgio::saveImg(dirName, "IVHGMRF", denoised_ivhgmrf);
			mylib::imgio::saveImg(dirName, "DVHGMRF", denoised_dvhgmrf);
		}

		// 画像の表示
		if (mylib::utillity::getImgVisualize()) {
			cv::imshow("IVGMRF", mylib::img::VtoM(denoised_ivgmrf));
			cv::imshow("DVGMRF", mylib::img::VtoM(denoised_dvgmrf));
			cv::imshow("IVHGMRF", mylib::img::VtoM(denoised_ivhgmrf));
			cv::imshow("DVHGMRF", mylib::img::VtoM(denoised_dvhgmrf));
			cv::waitKey(0);
		}

		// 進捗の表示
		mylib::utillity::progressBar();
	}
}