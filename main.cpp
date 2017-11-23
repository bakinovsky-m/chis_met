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


	auto a = pmAlgorythm(m2, m2[0]);
	Mat m_shit(m2.dim);
	vector<ld> x(m2.dim);
	for(int i = 0; i < x.size(); i++)
	{
		x[i] += 1.0 / (a.first * a.second[1]) * m2[1][i];
	}

	for(int i = 0; i < m_shit.dim; i++)
	{
		for(int j = 0; j < m_shit.dim; j++)
		{
			if(i == j)
			{
				m_shit[i][j] = m2[i][j]  - a.first * a.second[i] * x[i];
			}
			else
			{
				m_shit[i][j] = m2[i][j] - a.first * a.second[i] * x[j];
			}
		}
	}


	auto b = pmAlgorythm(m_shit, m_shit[0]);
	
	vector<long double> vec(b.second.size());
	long double ll = 0.0;
	for(int i = 0; i < b.second.size(); i++)
	{
		ll += b.second[i] * a.second[i];
	}
	for(int i = 0; i < vec.size(); i++)
	{
		vec[i] = (b.first - a.first) * b.second[i] + a.first * ll * a.second[i];
	}
	// std::cout << "vec2\n";
	// for(int i = 0; i < vec.size(); i++)
	// {
	// 	std::cout << vec[i] << " ";
	// }
	std::cout << std::endl;
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