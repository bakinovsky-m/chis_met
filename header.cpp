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
	const char * str = sstr.c_str();
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
	for (int i = 0; i < vec.size(); ++i){
		cout << vec[i] << " ";
	}
	cout << endl;
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
    int k, n;
    Mat new_matr = Mat();
    new_matr.resize(dim);
    for (int i = 0; i < dim; i++)
      new_matr[i].resize(dim);
    int det = 0;
    k = 1; //(-1) в степени i
    n = dim - 1;
    if (dim < 1) cout << "Определитель вычислить невозможно!";
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

    for (int ki = 0; ki < dimension; ki++) { // проверка индекса строки
        if (ki == i_ind) {
        	iOffset = 1;
        }
        jOffset = 0;
        for (int kj = 0; kj < dimension; kj++) { // проверка индекса столбца
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

void print_v(vector<ld> v) //Распечатка векторов
{
	for(int i = 0; i < v.size(); i++)
	{
		std::cout << v[i] << " ";
	}
	std::cout << std::endl;
}

ld norma(vector<ld> iter1, vector<ld> iter2) //разность векторов
{
	ld eps = 0.0; 
	for(int i = 0; i < iter1.size(); i++)
	{
		if(!isnan(iter1[i]) && !isnan(iter2[i]))
			eps += std::fabs(iter1[i] - iter2[i]);
	}
	return eps;
}

vector<ld> normal(const vector<ld> y) //нормирование вектора
{
	ld norma = 0.0;
	for(int i = 0; i < y.size(); i++)
	{
		norma += y[i] * y[i];
	}
	norma = std::pow(norma, 0.5);
	std::vector<ld> v;
	for(int i = 0; i < y.size(); i++)
	{
		v.push_back(y[i] / norma);
	}
	/*cout << "normal_v\n";
	for(int i = 0; i < y.size();i++)
		cout << y[i] << " ";
	cout << endl;*/
	return v;
}

vector<ld> div(const vector<ld> y, const vector<ld> x) //деление векторов
{
	vector<ld> res;
	for(int i = 0; i < y.size(); i++)
	{
		if( x[i] != 0)
			res.push_back(y[i] / x[i]);
		else
			res.push_back(0.0 / 0.0);
	}
	return res;	
}




pair<ld, vector<ld>> pmAlgorythm(const Mat& A0, const vector<ld> svobChleny) //собсна алгоритм
{
	vector <ld> lambda_prev(svobChleny.size());
	auto y0 = svobChleny;
	auto x0 = normal(y0);
	auto y1 = unmMatNaVec(A0, x0);
	auto lambda_current = div(y1, x0);
	while(norma(lambda_current, lambda_prev) > 0.000001)
	{

		y0 = y1;
		x0 = normal(y0);
		y1 = unmMatNaVec(A0, x0);
		lambda_prev = lambda_current;
		lambda_current = div(y1, x0);
	}
	x0 = normal(y1);
	ld res = 0.0;
	int count = 0;
	for (int i = 0; i < lambda_current.size(); i++)
	{
		if (lambda_current[i] != 0.0 / 0.0)
		{
			res += lambda_current[i];
			count++;
		}
	}
	if(count == 0)
	{
		std::cout << "Shit:(" << std::endl;
	}
	else 
	{
		std::cout << "eigenval: ";
		std::cout << res / count << std::endl;
		std::cout << "eigenvec: \n";
		for(int i = 0; i < x0.size(); i++)
			std::cout << x0[i] << " ";
		std::cout << std::endl<<endl;
		auto lll =  std::make_pair(res / count, x0);
		std::vector<ld> v = svobChleny;
		for(int i = 0; i < 10; i++)
		{
			v = unmMatNaVec(A0, v);
		}
		auto v1 = unmMatNaVec(A0, v);
		auto v2 = unmMatNaVec(A0, v1);
		vector<ld> y2(v2.size());
		for(int i = 0; i < y2.size(); i++)
		{
		 	y2[i] = (v2[i] - (res / count) * v1[i]) / (v1[i] - (res / count) * v[i]); 
		}

		return lll;
		// auto y2 = unmMatNaVec(A0, x0);
		// vector<ld> v(y2.size());
		// std::cout << "vectors" << std::endl;
		// print_v(y0);
		// print_v(y1);
		// print_v(y2);
		// std::cout << "here!!!\n";
		// for(int i = 0; i < y2.size(); i++)
		// {
		// 	v[i] = (y2[i] - (res / count) * y1[i]) / (y1[i] - (res / count) * y0[i]); 
		// }
		// res = 0.0;
		// count = v.size();
		// for(int i = 0; i < v.size(); i++)
		// {
		// 	cout << v[i];
		// 	res += v[i]; 
		// }
		// std::cout << "Second" << res << std::endl;
	}
	pair<ld, vector<ld>> l;
	return l;
} 

// vector<pair<ld, vector<ld>>> powMethod(const Mat& A0, const vector<ld> svobChleny){
// 	vector<pair<ld, vector<ld>>> res(0);
// 	pair<ld, vector<ld>> res1(0, vector<ld>(0));
// 	pair<ld, vector<ld>> res2(0, vector<ld>(0));
// 	ld eps = 1e-5;

// 	Mat A = A0;
// 	vector<ld> y0 = svobChleny;
// 	vector<ld> y1(y0.size(),0);

	

// 	bool done = false;

// 	ld y_k_prev_iter = 0;
// 	ld y_k = 0;

// 	ld a_plus_1 = 1;
// 	ld a = 1;
// 	ld a_minus_1 = 1;

// 	vector<ld> temp_cont(0);

// 	int counter = 0;
// 	done = false;
// 	while(!done){
// 	// while(fabs(a - y_k * a_minus_1) > eps){
// 		// cout << "fabs: " << fabs(a - y_k * a_minus_1) << endl;
// 		counter++;
// 		if(counter > 31){
// 			done = true;
// 		}
// 		y1 = unmMatNaVec(A, y0);

// 		y_k = y1[0] / y0[0];
// 		// cout<<"y_k = "<< y_k<<endl<<endl;
// 		a_minus_1 = a;
// 		a = a_plus_1;
// 		a_plus_1 = y0[0];

// 		ld temp = (a_plus_1 - (y_k * a)) / (a - (y_k * a_minus_1));


// 		// cout << "(" << a_plus_1 << " - " << y_k << " * " << a << ") / (" << a << " - " << y_k << " * " << a_minus_1 << ")" << endl;
// 		// cout << "temp: " << temp << endl;

// 		temp_cont.push_back(temp);

// 		y0 = y1;
// 	}

// 	y0 = svobChleny;
// 	done = false;
// 	y_k_prev_iter = 0;
// 	while(!done){
// 		// cout << "iter: " << counter << endl;

// 		y1 = unmMatNaVec(A, y0);
// 				for (int i = 0; i < y1.size(); ++i)
// 		y_k = y1[0] / y0[0];

// 		// cout<<"y1 = ";
// 		// for (int i = 0; i < y1.size(); ++i)
// 		// {
// 		// 	cout<<y1[i]<<" ";		
// 		// }
// 		// cout<<endl;
// 		// cout<<"y0 = ";
// 		// for (int i = 0; i < y0.size(); ++i)
// 		// {
// 		// 	cout<<y0[i]<<" ";		
// 		// }
// 		// cout<<endl<<endl;

// 		if(fabs(y_k - y_k_prev_iter) < eps){
// 			done = true;
// 			// break;
// 		}
// 		y_k_prev_iter = y_k;		
// 		y0 = y1;
// 	}

// 	for(int i = 0; i < y1.size(); ++i){
// 		y1[i] = y1[i] / y1[y1.size() - 1];
// 	}

// 	res1.first = y_k;

// 	// y0 = svobChleny;
// 	// y1 = 
// 	// while(!done)
// 	// {

// 	// }

// 	// res1.second = y1;

// 	// res2.first = (a_plus_1 - (res1.first * a)) / (a - (res1.first * a_minus_1));
// 	ld temp = 0;
// 	for(int i = 0; i < temp_cont.size(); ++i){
// 		// cout << "temp += " << temp_cont[i] << endl;
// 		temp += temp_cont[i];
// 	}
// 	// cout << temp << " / " << temp_cont.size() << endl;
// 	res2.first = temp / temp_cont.size();
// 	res2.second = methodGauss(A, res2.first);

// 	res.push_back(res1);
// 	res.push_back(res2);

// 	return res;
// }

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

	ld a_plus_1 = 1;
	ld a = 1;
	ld a_minus_1 = 1;

	vector<ld> temp_cont(0);

	int counter = 0;
	done = false;


	y0 = svobChleny;
	done = false;
	y_k_prev_iter = 0;
	while(!done){
		// cout << "iter: " << counter << endl;

		y1 = unmMatNaVec(A, y0);
				for (int i = 0; i < y1.size(); ++i)
		y_k = y1[0] / y0[0];

		// cout<<"y1 = ";
		// for (int i = 0; i < y1.size(); ++i)
		// {
		// 	cout<<y1[i]<<" ";		
		// }
		// cout<<endl;
		// cout<<"y0 = ";
		// for (int i = 0; i < y0.size(); ++i)
		// {
		// 	cout<<y0[i]<<" ";		
		// }
		// cout<<endl<<endl;

		if(fabs(y_k - y_k_prev_iter) < eps){
			done = true;
			// break;
		}
		y_k_prev_iter = y_k;

		a_minus_1 = a;
		a = a_plus_1;
		a_plus_1 = y0[0];

		y0 = y1;
	}

	for(int i = 0; i < y1.size(); ++i){
		y1[i] = y1[i] / y1[y1.size() - 1];
	}


	res1.first = y_k;

	// y0 = svobChleny;
	// y1 = 
	// while(!done)
	// {

	// }

	res1.second = y1;


	res2.second = methodGauss(A, res2.first);

	res.push_back(res1);
	res.push_back(res2);

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