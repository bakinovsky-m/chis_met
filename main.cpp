#include <iostream>
#include <vector>
#include "header.h"

using namespace std;

int main(){
	// 1 lab
	Mat m1;
	m1.fillFromFile("tests/for_first.txt");
	vector<ld> res = kramer(m1);

	cout << "!!! 1 zadanie !!!" << endl;
	print(m1);
	cout << endl;

	cout << "Kramer:" << endl;
	print(res);
	cout << endl;

	res = relaxation(m1, 1.1);
	cout << "Relaxation: " << endl;
	print(res);
	// end 1 lab

	// 2 lab
	Mat m2;
	m2.fillFromFile("tests/for_second.txt");

	cout << endl << "!!! 2 zadanie !!!" << endl << endl;
	print(m2);
	cout << endl;
	cout << "Obratnaya: " << endl;
	print(methodGaussJordan(m2));
	cout << endl;
	auto p = powMethod(m2, m2[0]);
	cout << "sobs ch: " << p[0].first << endl;
	cout << "sobs vec: " << endl;
	print(p[0].second);
	cout << endl;
	// end 2 lab

	// 3 lab 
	Mat m3;
	m3.fillFromFile("tests/for_third.txt");

	cout << "!!! 3 zadanie !!!" << endl << endl;
	print(m3);
	Mat L = Mat(m3.dim);
	Mat U = Mat(m3.dim);

	LUdecompos(m3, L, U);

	double eps = 1e-15;

	double pogr = 100;
	Mat temp = m3;
	Mat a1 = m3;
	while(abs(pogr) > eps) {
		LUdecompos(a1, L, U);
		temp = a1;
		a1 = umnMat(U, L);

		pogr = temp.mat[0][0] - a1.mat[0][0];
	}

	cout << endl;

	for (int i = 0; i < temp.dim; ++i) {
		for (int j = 0; j < temp.dim; ++j) {
			if (i == j){
				cout << "sobs chislo: " << temp.mat[i][j] << endl;
				auto sobs_vec = methodGauss(m3, temp.mat[i][j]);
				print(normal(sobs_vec));
				cout << endl;
			}
		}
	}
	// end 3 lab
	return 0;
}