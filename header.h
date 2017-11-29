#ifndef HEADER_H
#define HEADER_H

#include <vector>
#include <string>

typedef long double ld;

class Mat{
public:
	Mat() = default;
	Mat(int d);

	void resize(int d);

	std::vector<ld>& operator[](int row);
	std::vector<ld> operator[](int row) const;
	Mat operator*(Mat& rhs) const;

	void fillFromFile(std::string str);

	// fields
	int dim = 0;
	std::vector<ld> svobChleny;

	std::vector<ld> ownNumbers;
	std::vector<std::vector<ld>> ownVectors;

// private:
	std::vector<std::vector<ld>> mat;
};

Mat methodGaussJordan(const Mat m);
Mat multRow(const Mat& m, const int ind, const ld value);
Mat divRow(const Mat& m, const int ind, const ld value);
Mat swapRows(const Mat& m, const int ind1, const int ind2);

void print(const Mat m);
void print(const std::vector<ld> vec);

ld getDeterminant(const Mat old_matr);
Mat getMWithoutIJ(const Mat mas, const int i, const int j);
Mat getMWithIColChangedToSvobChlen(const Mat matr, const int index_stolba);
std::vector<ld> kramer(const Mat matr);
std::vector<ld> relaxation(const Mat A, const ld parameter);
void LUdecompos(const Mat A, Mat &L, Mat &U);

// std::pair<ld, std::vector<ld>> powMethod(const Mat& A, const std::vector<ld> svobChleny);
std::vector<std::pair<ld, std::vector<ld>>> powMethod(const Mat& A, const std::vector<ld> svobChleny);
std::vector<ld> unmMatNaVec(const Mat& A, const std::vector<ld>& vec);
std::vector<ld> methodGauss(const Mat& A0, const ld l);

Mat umnMat(const Mat & m1, const Mat & m2);

std::vector<ld> normal(const std::vector<ld> y);
#endif
