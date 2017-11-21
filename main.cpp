#include <iostream>
#include <vector>
#include "header.h"

using namespace std;

int main(){
	// 1 lab
	Mat m1;
	m1.fillFromFile("tests/for_first.txt");
	vector<ld> res = kramer(m1);

	cout << "1 zadanie" << endl;
	print(m1);
	cout << endl;

	cout << "Kramer:" << endl;
	// cout << "res: " << endl;
	for(int i = 0; i < res.size(); ++i){
		cout << res[i] << " ";
	}
	cout << endl;
	cout << endl;

	res = relaxation(m1, 1.1);
	cout << "Relaxation: " << endl;
	for(int i = 0; i < res.size(); ++i){
		cout << res[i] << " ";
	}
	cout << endl;
	// end 1 lab

	// 2 lab
	Mat m2;
	m2.fillFromFile("tests/for_second.txt");

	cout << endl << "2 zadanie" << endl;
	print(m2);
	cout << endl;
	cout << "Obratnaya: " << endl;
	print(methodGaussJordan(m2));
	cout << endl;
	// powMethod(m2, m2.svobChleny);
	auto p = powMethod(m2, m2[0]);
	cout << "sobs ch 1: " << p[0].first << endl;
	cout << "sobs vec 1: " << endl;
	print(p[0].second);
	cout << endl;

	cout << "sobs ch 2: " << p[1].first << endl;
	cout << "sobs vec 2: " << endl;
	print(p[1].second);
	cout << endl;
	// end 2 lab

	// 3 lab 

	// end 3 lab
	return 0;
}