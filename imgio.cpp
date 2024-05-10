#include "imgio.hpp"

// Mat -> vector
matrix imgio::MtoV(const cv::Mat& cvmat) {
	matrix vecs(cvmat.rows, std::vector<double>(cvmat.cols));
	for (int i = 0; i < cvmat.rows; i++) {
		for (int j = 0; j < cvmat.cols; j++) {
			auto x = double(cvmat.at<unsigned char>(i, j));
			vecs[i][j] = x;
		}
	}
	return vecs;
}

// vector -> Mat
cv::Mat imgio::VtoM(const matrix& vecs) {
	//vec�̑傫�����擾
	uint16_t vec_rows = static_cast<uint16_t>(vecs.at(0).size());	//�s
	uint16_t vec_cols = static_cast<uint16_t>(vecs.size());	//��
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

// �K�E�X�m�C�Y�̕t�^
void imgio::corruption(matrix& src)
{
	uint16_t rows = static_cast<int>(src.at(0).size());	//�s
	uint16_t cols = static_cast<int>(src.size());	//��

	std::random_device seed_gen;	//�V�[�h�̏�����
	std::default_random_engine gen(seed_gen());	//��������
	std::normal_distribution<double> dist(this->m, this->stddev);	//�K�E�X���z

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			const double tmp = dist(gen);
			src[i][j] += tmp;
		}
	}
	imgio::limitation(src);
}

void imgio::limitation(matrix& img) 
{
	uint16_t rows = static_cast<int>(img.at(0).size());	//�s
	uint16_t cols = static_cast<int>(img.size());	//��

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (img[i][j] <= 0) {
				img[i][j] = 0;
			}
			else if (img[i][j] >= 255) {
				img[i][j] = 255;
			}
		}
	}
}

double imgio::MSE(const matrix& img1, const matrix& img2)
{
	uint16_t rows = static_cast<int>(img1.at(0).size());	//�s
	uint16_t cols = static_cast<int>(img1.size());	//��
	uint32_t n = rows * cols;
	double sum = 0.0;

	// �s��̃T�C�Y���Ⴄ�ꍇ�̗�O����
	if (rows != img2.at(0).size() || cols != img2.at(1).size())
	{
		std::cout << "vector size error : [" << rows << "," << cols << "] to["
			<< img2.at(0).size() << "," << img2.at(1).size() << std::endl;
		return 0;
	}

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			sum += std::pow(img1[i][j] - img2[i][j], 2);
		}
	}
	return sum / n * 1.0;
}

double imgio::PSNR(const matrix& origin, const matrix& corrupt)
{
	return static_cast<double>(log10(pow(255, 2) / MSE(origin, corrupt)) * 10.0);
}

void imgio::setEnumerate(uint32_t _enumerate)
{
	this->enumerate = _enumerate;
}

void imgio::setM(double _m)
{
	this->m = _m;
}

void imgio::setStddev(double _stddev)
{
	this->stddev = _stddev;
}