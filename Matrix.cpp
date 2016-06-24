#include "stdafx.h"
#include <iostream>
#include<math.h>
using namespace std;
#include<vector>
#include <string>
#include <fstream>
#include <sstream>

#include "Matrix.h"

Matrix::Matrix(int d11, int d22)
{
	rows = d11;
	cols = d22;
	vector<vector<double>> p1(rows, vector<double>(cols, 0));
	mat = p1;
}

Matrix::Matrix()
{
	rows = 0;
	cols = 0;
	vector<vector<double>> p1(rows, vector<double>(cols, 0));
	mat = p1;
}

Matrix::~Matrix()
{
	vector<vector<double>>().swap(mat);
}


Matrix Matrix::mutiply(Matrix &mat1, Matrix &mat2)
{
	int r1 = mat1.rows;
	int c1 = mat1.cols;
	int r2 = mat2.rows;
	int c2 = mat2.cols;

	Matrix ans(r1, c2);

	for (int i = 0; i<r1; ++i)
		for (int j = 0; j<c2; ++j)
			for (int k = 0; k<c1; ++k)
			{
				ans.mat[i][j] += mat1.mat[i][k] * mat2.mat[k][j];
			}
	return ans;
}

void Matrix::print(string matName)
{
	cout << "=============================" << endl;
	if (matName != "")
	{
		cout << matName << endl;
	}

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}
}

void Matrix::resize(int _rows, int _cols)
{

	cols = _cols;
	rows = _rows;
	if (rows == 0 && cols == 0)
	{
		mat.resize(_cols, vector<double>(_rows));

	}
	else
	{
		mat.resize(_rows);
		for (int i = 0; i < rows; i++)
		{
			mat[i].resize(_cols);
		}

	}
}

void Matrix::ones()
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			mat[i][j] = 1.0;
		}
	}

}

double Matrix::Det_3x3()
{
	double determinant;
	switch (rows)
	{
	case 2:
		determinant = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
		break;
	case 3:
		determinant = mat[0][0] * ((mat[1][1] * mat[2][2]) - (mat[2][1] * mat[1][2])) - mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2])
			+ mat[0][2] * (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]);
		break;
	default:
		printf("Only 2X2 and 3X3 matrices ");
		break;
	}
	return determinant;
}

vector<double> Vector1D::subtract(vector<double> arrA, vector<double> arrB)
{
	int sizeA = arrA.size();
	int sizeB = arrB.size();
	vector<double> ans(sizeA);

	if (sizeA == sizeB)
	{
		for (int i = 0; i < sizeA; i++)
		{
			ans[i] = arrA[i] - arrB[i];
		}
	}
	return ans;
}

double Vector1D::Abs(vector<double> arr)
{
	double ans = 0;
	int size = arr.size();
	for (int i = 0; i < size; i++)
	{
		ans += pow(arr[i], 2.0);

	}
	ans = sqrt(ans);
	return ans;
}

vector<double> Vector1D::sum(vector<double> arrA, vector<double> arrB)
{
	int sizeA = arrA.size();
	int sizeB = arrB.size();
	vector<double> ans(sizeA);

	if (sizeA == sizeB)
	{
		for (int i = 0; i < sizeA; i++)
		{
			ans[i] = arrA[i] + arrB[i];
		}
	}
	return ans;
}

vector<double> Vector1D::multiScal(vector<double> arr, double scal)
{
	int size = arr.size();
	vector<double>ans(size);
	for (int i = 0; i < size; i++)
	{
		ans[i] = arr[i] * scal;
	}
	return ans;
}

vector<double> Vector1D::crossProduct(vector<double> u, vector<double> v)
{
	int sizeA = u.size();
	int sizeB = v.size();
	vector<double>ans(3);
	ans[0] = 0;
	ans[1] = 0;
	ans[2] = 0;
	if ((sizeA == 3) && (sizeB == 3))
	{
		ans[0] = u[1] * v[2] - u[2] * v[1];
		ans[1] = u[2] * v[0] - u[0] * v[2];
		ans[2] = u[0] * v[1] - u[1] * v[0];
	}
	return ans;
}


double Vector1D::distance(vector<double> P1, vector<double> P2)
{
	double ans = 0;
	for (int i = 0; i < P2.size(); i++)
	{
		ans += pow(P2[i] - P1[i], 2.0);

	}
	ans = sqrt(ans);

	return ans;
}
