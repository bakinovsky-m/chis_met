#include "header.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

typedef long double ld;

Mat::Mat(int d){
	this->dim = d;
	this->mat.resize(dim);
	for (int i = 0; i < this->dim; ++i)
	{
		this->mat[i] = vector<ld>(this->dim, 0);
	}
}

vector<ld>& Mat::operator[](int row) {
	return mat[row];
}

vector<ld> Mat::operator[](int row) const{
	return mat[row];	
}

Mat Mat::operator*(Mat& rhs) const {
	Mat res = rhs;
	for (int i = 0; i < rhs.dim; ++i) {
		for(int j = 0; j < rhs.dim; ++j) {
			res[i][j] = 0;
			for(int k = 0; k < rhs.dim; ++k){
				res[i][j] += this->mat[i][k] * rhs[k][j];
			}
		}
	}
	return res;
}

void Mat::resize(int d){
	this->dim = d;
	mat.resize(dim);
	for (int i = 0; i < dim; ++i)
	{
		mat[i].resize(dim);
	}
}

void Mat::fillFromFile(std::string sstr){
	int dim = 123;
	ifstream fin(sstr, ifstream::in);
	fin >> dim;
	this->resize(dim);
	this->dim = dim;
	for(int i = 0; i < dim; ++i){
		this->mat[i].resize(dim);
		for (int j = 0; j < dim; ++j)
		{
			fin >> this->mat[i][j];
		}
	}
	if(fin.good()){
		for(int i = 0; i < dim; ++i){
			this->svobChleny.resize(dim);
			fin >> this->svobChleny[i];
		}
	}
	fin.close();
}


Mat methodGaussJordan(const Mat old_m){
	Mat m = old_m;
	Mat res = old_m;
	for(int i = 0; i < old_m.dim; ++i){
		for(int j = 0; j < old_m.dim; ++j){
			if(i == j){
				res[i][j] = 1;
			} else {
				res[i][j] = 0;
			}
		}
	}

	for(int str_num = 0; str_num < m.dim; ++str_num){
		ld del = m[str_num][str_num];
		m = divRow(m, str_num, del);
		res = divRow(res, str_num, del);

		for(int n = str_num + 1; n < m.dim; ++n){
			ld u_m = m[n][str_num];

			for(int el = 0; el < m.dim; ++el){
				m[n][el] -= m[str_num][el] * u_m;
			}

			for(int el = 0; el < res.dim; ++el){
				res[n][el] -= res[str_num][el] * u_m;
			}
		}
	}

	for(int str_num = m.dim - 1; str_num >= 0; --str_num){
		for(int n = str_num - 1; n >= 0; --n){
			ld u_m = m[n][str_num];

			for(int el = 0; el < m.dim; ++el){
				m[n][el] -= m[str_num][el] * u_m;
			}

			for(int el = 0; el < res.dim; ++el){
				res[n][el] -= res[str_num][el] * u_m;
			}
		}
	}
	return res;
}

Mat swapRows(const Mat& m, const int ind1, const int ind2){
	Mat res = m;

	vector<ld> temp_row = res[ind1];
	res[ind1] = res[ind2];
	res[ind2] = temp_row;

	return res;
}

Mat multRow(const Mat& m, const int ind, const ld value){
	Mat res = m;
	for(int i = 0; i < res.dim; ++i){
		res[ind][i] *= value;
	}
	return res;
}

Mat divRow(const Mat& m, const int ind, const ld value){
	if(value == 0) throw invalid_argument("div by 0");

	Mat res = m;
	for(int i = 0; i < res.dim; ++i){
		res[ind][i] /= value;
	}
	return res;
}

void print(const Mat m){
	for(int i = 0; i < m.dim; ++i){
		for(int j = 0; j < m.dim; ++j){
			cout << m[i][j] << " ";
			if(m.svobChleny.size() > 0 && j == m.dim - 1) {
				cout << "|" << m.svobChleny[i];
			}
		}
		cout << endl;
	}
}

void print(const vector<ld> vec){
	for (unsigned int i = 0; i < vec.size(); ++i){
		cout << vec[i] << endl;
	}
}

vector<ld> kramer(const Mat matr){
	vector<ld> result;
	result.resize(matr.dim);
	for(int i = 0; i < matr.dim; i++) {
		Mat p = getMWithIColChangedToSvobChlen(matr, i);
		if(getDeterminant(matr) != 0){
			result[i] =  1.0f/getDeterminant(matr) * getDeterminant(p);
		} else {
			result[i] = -1;
		}
	}
	return result;
}


ld getDeterminant(const Mat old_matr) {
	int dim = old_matr.dim;
    int k = 0;
    Mat new_matr = Mat();
    new_matr.resize(dim);
    for (int i = 0; i < dim; i++)
      new_matr[i].resize(dim);
    int det = 0;
    k = 1; 
    if (dim < 1) cout << "No determinant for you!";
    if (dim == 1) {
      det = old_matr[0][0];
      return(det);
    } else if (dim == 2) {
      det = old_matr[0][0] * old_matr[1][1] - (old_matr[1][0] * old_matr[0][1]);
      return(det);
    } else if (dim>2) {
      for (int i = 0; i<dim; i++) {
        new_matr = getMWithoutIJ(old_matr, i, 0);
        det = det + k * old_matr[i][0] * getDeterminant(new_matr);
        k = -k;
      }
    }
    return det;
}

Mat getMWithIColChangedToSvobChlen(const Mat matr, const int index_stolba){
	Mat p;
	int dimension = matr.dim;
	p.resize(dimension);
	for(int i = 0; i < dimension; ++i){
		p[i].resize(dimension);
	}
	for(int i = 0; i < dimension; i++){
		for(int j = 0; j < dimension; j++){
			if(j == index_stolba){
				p[i][j] = matr.svobChleny[i];
			} else {
				p[i][j] = matr[i][j];
			}
		}
	}

	return p;
}

Mat getMWithoutIJ(const Mat matr, const int i_ind, const int j_ind) {
	int iOffset, jOffset;
    iOffset = 0;
    int dimension = matr.dim - 1;

    Mat p;
	p.resize(dimension);

    for (int ki = 0; ki < dimension; ki++) {
        if (ki == i_ind) {
        	iOffset = 1;
        }
        jOffset = 0;
        for (int kj = 0; kj < dimension; kj++) {
            if (kj == j_ind) {
            	jOffset = 1;
            }
            p[ki][kj] = matr[ki + iOffset][kj + jOffset];
	    }
	}

	return p;
}

vector<ld> relaxation(const Mat A, const ld parameter) {
    vector<ld> cur = A.svobChleny;
    int n = A.dim;
    ld checker = 1e-5;
    // int iterStep = 0;

    while (true) {
    	// iterStep++;
        vector<ld> old = cur;
        for (int i = 0;i < n; ++i) {
            ld sum = 0;
            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                sum -= (A[i][j] * cur[j]);
            }
            sum += A.svobChleny[i];
            sum *= (parameter / A[i][i]);
            cur[i] = sum + (1 - parameter) * cur[i];
        }
        ld mx = abs(cur[0] - old[0]);
        for (int i = 0; i < n; ++i) {
            mx = max(abs(cur[i] - old[i]), mx);
        }
        if (mx < checker) break;
    }
    // cout << iterStep << " steps" << endl;
    return cur;
}

void LUdecompos(const Mat A, Mat &L, Mat &U){
	int n = A.dim;
	U = A;

	for(int i = 0; i < n; i++)
		for(int j = i; j < n; j++)
			L.mat[j][i] = U.mat[j][i] / U.mat[i][i];
	
	for(int k = 1; k < n; k++)
	{
		for(int i = k-1; i < n; i++)
			for(int j = i; j < n; j++)
				L.mat[j][i] = U.mat[j][i] / U.mat[i][i];

		for(int i = k; i < n; i++)
			for(int j = k-1; j < n; j++)
				U.mat[i][j] = U.mat[i][j] - L.mat[i][k-1] * U.mat[k-1][j];
	}

}

vector<pair<ld, vector<ld>>> powMethod(const Mat& A0, const vector<ld> svobChleny){
	vector<pair<ld, vector<ld>>> res(0);
	pair<ld, vector<ld>> res1(0, vector<ld>(0));
	pair<ld, vector<ld>> res2(0, vector<ld>(0));
	ld eps = 1e-5;

	Mat A = A0;
	vector<ld> y0 = svobChleny;
	vector<ld> y1(y0.size(),0);

	bool done = false;

	ld y_k_prev_iter = 0;
	ld y_k = 0;

	y_k_prev_iter = 0;

	while(!done){
		y1 = unmMatNaVec(A, y0);

		y_k = y1[0] / y0[0];

		if(fabs(y_k - y_k_prev_iter) < eps){
			done = true;
		}

		y_k_prev_iter = y_k;

		y0 = y1;
	}

	for(unsigned int i = 0; i < y1.size(); ++i){
		y1[i] = y1[i] / y1[y1.size() - 1];
	}

	res1.first = y_k;
	res1.second = y1;

	res.push_back(res1);

	return res;
}

vector<ld> unmMatNaVec(const Mat& A, const vector<ld>& vec){
	vector<ld> res(0);

	for(int i = 0; i < A.dim; ++i){
		ld t = 0;
		for(int k = 0; k < A.dim; ++k){
			t += A[i][k] * vec[k];
		}
		res.push_back(t);
	}

	return res;
}

vector<ld> methodGauss(const Mat& A0, const ld l){
	vector<ld> res(A0.dim, 0);
	res[res.size() - 1] = 1;
	Mat A = A0;

	for (int i = 0; i < A.dim; ++i){
		for(int j = 0; j < A.dim; ++j){
			if(i == j){
				A[i][j] -= l;
			}
		}
	}

	for(int str_num = 0; str_num < A.dim; ++str_num){
		ld del = A[str_num][str_num];
		A = divRow(A, str_num, del);

		for(int n = str_num + 1; n < A.dim; ++n){
			ld u_m = A[n][str_num];

			for(int el = 0; el < A.dim; ++el){
				A[n][el] -= A[str_num][el] * u_m;
			}
		}
	}


	for(int i = A.dim - 1; i >= 0; --i){
		for(int j = i + 1; j < A.dim; ++j){
			res[i] -= res[j] * A[i][j];
		}
		res[i] /= A[i][i];
	}

	return res;
}

Mat umnMat(const Mat & m1, const Mat & m2){
	Mat temp = m1;
	for(int i = 0; i < m1.dim; ++i){
		for(int j = 0; j < m1.dim; ++j){
			temp.mat[i][j] = 0;
			for(int k = 0; k < m1.dim; ++k){
				temp.mat[i][j] += m1.mat[i][k] * m2.mat[k][j];
			}
		}
	}
	return temp;
}