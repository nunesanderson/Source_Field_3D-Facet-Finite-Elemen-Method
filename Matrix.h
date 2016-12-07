#include <string>
#include<vector>
using namespace std;
using std::string;

#ifndef _MATRIX_INTERNAL_INCLUDED_
#define _MATRIX_INTERNAL_INCLUDED_

class Matrix
{

public:
	int rows, cols;
	Matrix();
	~Matrix();
	Matrix(int,int);
	vector<vector<double>> mat;
	Matrix mutiply(Matrix &mat1, Matrix &mat2);
	void print(string matName = "");
	void resize(int rows,int cols);
	void ones();
	double Det_3x3();
	double Det_2x2();


private:

};

class Vector1D
{
public:
	vector <double> subtract(vector <double> arrA, vector <double> arrB);
	double Abs(vector <double> arr);
	vector <double> sum(vector <double> arrA, vector <double> arrB);
	vector <vector <double>> sum(vector <vector <double>> arrA, vector <vector <double>> arrB);
	vector<double> multiScal(vector<double>arr, double scal);
	vector<double> crossProduct(vector <double> u, vector <double> v);
	double distance(vector <double> u, vector <double> v);
	double dot(vector <double> u, vector <double> v);

	
};


#endif