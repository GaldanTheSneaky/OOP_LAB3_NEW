#include"inherit_vector.h"
#include<iostream>

int main()
{
	tnzr::Matrix<double> matrix;
	tnzr::Vector<double> vector1 = { 1, 2, 3 };
	tnzr::Vector<double> vector2 = { 4, 5, 6 };
	tnzr::Vector<double> vector3 = { 4, 5, 6 };
	vector3 = vector1 + vector2;
	std::cout << vector3;
	return 0;
}