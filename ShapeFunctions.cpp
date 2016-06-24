#include "stdafx.h"
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <limits>
using namespace std;


/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include "ShapeFunctions.h"
#include "Messages.h"
#include "Gmsh.h"


Messages msg;

ShapeFunctions::~ShapeFunctions() {}
/* ------------------------------------------------------------------------
Nodal shape functions

[N1 N2...Nn]
---------------------------------------------------------------------------*/

Matrix ShapeFunctions::GetNodalShapeFunction(int ElemType, double u, double v, double p) {
	Matrix ans(0,0);

	switch (ElemType)
	{
	case 1: //First order line
		ans = nodalLineFirstOrder(u);
		break;

	case 4: //First order tetrahedral
		ans=nodalTetrahedralFirstOrder(u, v, p);
		break;

	case 5://First order hexahedral
		ans=nodalHexahedralFirstOrder(u, v, p);
		break;

	case 8: //Second order line
		ans=nodalLineSecondOrder(u);
		break;
	case 11: // Second order tetrahedral
		ans=nodalTetrahedralSecondOrder(u, v, p);
		break;

	default:
		msg.NotImplementedElement(ElemType, "GetNodalShapeFunction");
	}
	return ans;
}

Matrix ShapeFunctions::nodalLineFirstOrder(double u) {
	Matrix ans(1, 2);
	ans.mat[0][0] = 0.5 * (1.0 - u);
	ans.mat[0][1] = 0.5 * (1.0 + u);
	return ans;
}

/*
Second order line
0    2    1
*----*----*
*/
Matrix ShapeFunctions::nodalLineSecondOrder(double u) {
	Matrix ans(1, 3);

	ans.mat[0][0] = -u*(1.0 - u) / 2.0;
	ans.mat[0][1] = u*(1.0 + u) / 2.0;
	ans.mat[0][2] = (1.0 - u*u);
	return ans;
}

Matrix ShapeFunctions::nodalTetrahedralFirstOrder(double u, double v, double p) {
	// Ref Ida&Bastos, pg 319, table 8.9
	Matrix ans(1, 4);

	ans.mat[0][0] = 1.0 - u - v - p;
	ans.mat[0][1] = u;
	ans.mat[0][2] = v;
	ans.mat[0][3] = p;

	return ans;
}
Matrix ShapeFunctions::nodalTetrahedralSecondOrder(double u, double v, double p)
{
	Matrix ans(1, 10);
	double t = 1.0 - u - v - p;

	ans.mat[0][0] = -t*(1.0 - 2.0*t);
	ans.mat[0][1] = 4.0*u*t;
	ans.mat[0][2] = -u*(1.0 - 2.0*u);
	ans.mat[0][3] = 4.0*u*v;
	ans.mat[0][4] = -v*(1.0 - 2.0*v);
	ans.mat[0][5] = 4.0*v*t;
	ans.mat[0][6] = 4.0*p*t;
	ans.mat[0][7] = 4.0*u*p;
	ans.mat[0][8] = 4.0*v*p;
	ans.mat[0][9] = -p*(1.0 - 2.0*p);


	return ans;
}
Matrix ShapeFunctions::nodalHexahedralFirstOrder(double u, double v, double p) {
	//Ref Ida&Bastos, pg 321, table 8.12
	Matrix ans( 1, 8);


	double a1 = 1.0 + u;
	double a2 = 1.0 - u;
	double b1 = 1.0 + v;
	double b2 = 1.0 - v;
	double c1 = 1.0 + p;
	double c2 = 1.0 - p;
	ans.mat[0][0] = a2*b2*c2 / 8.0;
	ans.mat[0][1] = a1*b2*c2 / 8.0;
	ans.mat[0][2] = a1*b1*c2 / 8.0;
	ans.mat[0][3] = a2*b1*c2 / 8.0;
	ans.mat[0][4] = a2*b2*c1 / 8.0;
	ans.mat[0][5] = a1*b2*c1 / 8.0;
	ans.mat[0][6] = a1*b1*c1 / 8.0;
	ans.mat[0][7] = a2*b1*c1 / 8.0;
	return ans;
}

/* ------------------------------------------------------------------------
Gradient of nodal shape functions
[dN1/du		dN2/du	... dNn/du,
dN1/dv		dN2/dv	... dNn/dv,
dN1/dp		dN2/dp	... dNn/dp]
---------------------------------------------------------------------------*/

Matrix ShapeFunctions::GetGradNodalShapeFunction(int ElemType, double u, double v, double p) {
	Matrix ans(0,0);


	switch (ElemType)
	{
	case 1: //First order line
		ans=gradNodalLineFirstOrder();
		break;

	case 4: //First order tetrahedral
		ans=gradNodalTetrahedralFirstOrder();
		break;

	case 5://First order hexahedral
		ans = gradNodalHexahedralFirstOrder(u, v, p);
		break;

	case 8://Second order line
		ans=gradNodalLineSecondOrder(u);
		break;
	case 11://Second order tetrahedral
		ans=gradNodalTetrahedralSecondOrder(u, v, p);
		break;
	default:
		msg.NotImplementedElement(ElemType, "GetGradNodalShapeFunction");

	}

	return ans;
}

Matrix ShapeFunctions::gradNodalLineFirstOrder() {
	Matrix ans( 3, 2);
	ans.mat[0][0] = -0.5;	ans.mat[0][1] = 0.5;
	ans.mat[1][0] = 0.0;	ans.mat[1][1] = 0.0;
	ans.mat[2][0] = 0.0;	ans.mat[2][1] = 0.0;

	return ans;

}

Matrix ShapeFunctions::gradNodalLineSecondOrder(double u) {
	Matrix ans (3, 3);
	ans.mat[0][0] = (-1.0 + 2 * u) / 2.0;	ans.mat[0][1] = (1.0 + 2 * u) / 2.0;	ans.mat[0][2] = -2 * u;
	ans.mat[1][0] = 0.0;					ans.mat[1][1] = 0.0;					ans.mat[1][2] = 0.0;
	ans.mat[2][0] = 0.0;					ans.mat[2][1] = 0.0;					ans.mat[2][2] = 0.0;

	return ans;

}
Matrix ShapeFunctions::gradNodalTetrahedralFirstOrder() {

	/*ref:
	Ida&Bastos, pg 319, table 8.9
	First order tetrahedral element*/

	Matrix ans(3, 4);
	ans.mat[0][0] = -1.0;	ans.mat[0][1] = 1.0;	ans.mat[0][2] = 0.0;	ans.mat[0][3] = 0.0;
	ans.mat[1][0] = -1.0;	ans.mat[1][1] = 0.0;	ans.mat[1][2] = 1.0;	ans.mat[1][3] = 0.0;
	ans.mat[2][0] = -1.0;	ans.mat[2][1] = 0.0;	ans.mat[2][2] = 0.0;	ans.mat[2][3] = 1.0;

	return ans;

}
Matrix ShapeFunctions::gradNodalTetrahedralSecondOrder(double u, double v, double p)
{
	Matrix ans(3, 10);
	double t = 1.0 - u - v - p;

	//gradNu
	ans.mat[0][0] = 1.0 - 4.0 * t;
	ans.mat[0][1] = 4.0 * (t - u);
	ans.mat[0][2] = -1.0 + 4.0 * u;
	ans.mat[0][3] = 4.0 * v;
	ans.mat[0][4] = 0.0;
	ans.mat[0][5] = -4.0 * v;
	ans.mat[0][6] = -4.0 * p;
	ans.mat[0][7] = 4.0 * p;
	ans.mat[0][8] = 0.0;
	ans.mat[0][9] = 0.0;

	//gradNv
	ans.mat[1][0] = 1.0 - 4.0 * t;
	ans.mat[1][1] = -4.0 * u;
	ans.mat[1][2] = 0.0;
	ans.mat[1][3] = 4.0 * u;
	ans.mat[1][4] = -1.0 + 4.0 * v;
	ans.mat[1][5] = 4.0 * (t - v);
	ans.mat[1][6] = -4.0 * p;
	ans.mat[1][7] = 0.0;
	ans.mat[1][8] = 4.0* p;
	ans.mat[1][9] = 0.0;

	//gradNp
	ans.mat[2][0] = 1.0 - 4.0*t;
	ans.mat[2][1] = -4.0*u;
	ans.mat[2][2] = 0.0;
	ans.mat[2][3] = 0.0;
	ans.mat[2][4] = 0.0;
	ans.mat[2][5] = -4.0 * v;
	ans.mat[2][6] = 4.0 * (t - p);
	ans.mat[2][7] = 4.0 * u;
	ans.mat[2][8] = 4.0*v;
	ans.mat[2][9] = -1.0 + 4.0*p;

	return ans;
}
Matrix ShapeFunctions::gradNodalHexahedralFirstOrder(double u, double v, double p) {

	/*ref:
	Ida&Bastos, pg 319, table 8.9*/

	Matrix ans(3, 8);

	double a1 = 1.0 + u;
	double a2 = 1.0 - u;
	double b1 = 1.0 + v;
	double b2 = 1.0 - v;
	double c1 = 1.0 + p;
	double c2 = 1.0 - p;
	double  k = 1.0 / 8.0;

	ans.mat[0][0] = -b2*c2;	ans.mat[0][1] = b2*c2;	ans.mat[0][2] = b1*c2;	ans.mat[0][3] = -b1*c2;
	ans.mat[1][0] = -a2*c2;	ans.mat[1][1] = -a1*c2;	ans.mat[1][2] = a1*c2;	ans.mat[1][3] = a2*c2;
	ans.mat[2][0] = -a2*b2;	ans.mat[2][1] = -a1*b2;	ans.mat[2][2] = -a1*b1;	ans.mat[2][3] = -a2*b1;

	ans.mat[0][4] = -b2*c1;	ans.mat[0][5] = b2*c1;	ans.mat[0][6] = b1*c1;	ans.mat[0][7] = -b1*c1;
	ans.mat[1][4] = -a2*c1;	ans.mat[1][5] = -a1*c1;	ans.mat[1][6] = a1*c1;	ans.mat[1][7] = a2*c1;
	ans.mat[2][4] = a2*b2;	ans.mat[2][5] = a1*b2;	ans.mat[2][6] = a1*b1;	ans.mat[2][7] = a2*b1;

	return ans;

}

/* ------------------------------------------------------------------------
Gauss Points
u	v	p
point 1
point 2
.
.
.
point n
---------------------------------------------------------------------------*/

GaussLegendrePoints::~GaussLegendrePoints() {}

GaussLegendrePoints::GaussLegendrePoints(int ElemType) {

	switch (ElemType)
	{
	case 1: //First order line
		lineOnePoint();
		break;

	case 4: //First order tetrahedral
		tetrahedralFourPoinst();
		break;

	case 5: //First order hexahedral
		hexahedralEightPoinstInside();
		break;

	case 8://Second order line
		lineTwoPoints();
		break;
	case 11: //First order tetrahedral
		tetrahedralFivePoinst();
		break;

	default:
		msg.NotImplementedElement(ElemType, "GaussLegendrePoints");

	}
}

void GaussLegendrePoints::lineTwoPoints() {

	double k = sqrt(1.0 / 3.0);
	pointsCoordinates.resize(2, 3);

	pointsCoordinates.mat[0][0] = -k;	pointsCoordinates.mat[0][1] = 0.0;	pointsCoordinates.mat[0][2] = 0.0;
	pointsCoordinates.mat[1][0] = k;	pointsCoordinates.mat[1][1] = 0.0;	pointsCoordinates.mat[1][2] = 0.0;

	weights.push_back(1);
	weights.push_back(1);
}

void GaussLegendrePoints::lineOnePoint() {
	pointsCoordinates.resize(1, 3);
	pointsCoordinates.mat[0][0] = 0.0;	pointsCoordinates.mat[0][1] = 0.0;	pointsCoordinates.mat[0][2] = 0.0;
	
	weights.push_back(1);
}

void GaussLegendrePoints::tetrahedralOnePoint() {

	double k = 1.0 / 4.0;
	pointsCoordinates.resize(1, 3);

	pointsCoordinates.mat[0][0] = k;	pointsCoordinates.mat[0][1] = k;	pointsCoordinates.mat[0][2] = k;

	weights.push_back(1.0 / 6.0);
}

void GaussLegendrePoints::tetrahedralFourPoinst() {
	double a = 0.1381966;
	double b = 0.5854102;
	pointsCoordinates.resize(4, 3);

	pointsCoordinates.mat[0][0] = a;	pointsCoordinates.mat[0][1] = a;	pointsCoordinates.mat[0][2] = a;
	pointsCoordinates.mat[1][0] = a;	pointsCoordinates.mat[1][1] = a;	pointsCoordinates.mat[1][2] = b;
	pointsCoordinates.mat[2][0] = a;	pointsCoordinates.mat[2][1] = b;	pointsCoordinates.mat[2][2] = a;
	pointsCoordinates.mat[3][0] = b;	pointsCoordinates.mat[3][1] = a;	pointsCoordinates.mat[3][2] = a;

	double k = 1.0 / 24.0;
	weights.push_back(k);
	weights.push_back(k);
	weights.push_back(k);
	weights.push_back(k);
}

void GaussLegendrePoints::tetrahedralFivePoinst()
{
	double a = 1.0 / 4.0;
	double b = 1.0 / 6.0;
	double c = 1.0 / 2.0;
	pointsCoordinates.resize(5, 3);

	pointsCoordinates.mat[0][0] = a;	pointsCoordinates.mat[0][1] = a;	pointsCoordinates.mat[0][2] = a;
	pointsCoordinates.mat[1][0] = b;	pointsCoordinates.mat[1][1] = b;	pointsCoordinates.mat[1][2] = b;
	pointsCoordinates.mat[2][0] = b;	pointsCoordinates.mat[2][1] = b;	pointsCoordinates.mat[2][2] = c;
	pointsCoordinates.mat[3][0] = b;	pointsCoordinates.mat[3][1] = c;	pointsCoordinates.mat[3][2] = b;
	pointsCoordinates.mat[4][0] = c;	pointsCoordinates.mat[4][1] = b;	pointsCoordinates.mat[4][2] = b;


	double k = 3.0 / 40.0;
	weights.push_back(-2.0 / 15.0);
	weights.push_back(k);
	weights.push_back(k);
	weights.push_back(k);
	weights.push_back(k);


}


void GaussLegendrePoints::hexahedralEightPoinstInside() {
	double r_3 = 1.0 / sqrt(3.0);
	pointsCoordinates.resize(8, 3);

	pointsCoordinates.mat[0][0] = r_3;	pointsCoordinates.mat[0][1] = r_3;	pointsCoordinates.mat[0][2] = r_3;
	pointsCoordinates.mat[1][0] = -r_3;	pointsCoordinates.mat[1][1] = r_3;	pointsCoordinates.mat[1][2] = r_3;
	pointsCoordinates.mat[2][0] = r_3;	pointsCoordinates.mat[2][1] = -r_3;	pointsCoordinates.mat[2][2] = r_3;
	pointsCoordinates.mat[3][0] = -r_3;	pointsCoordinates.mat[3][1] = -r_3;	pointsCoordinates.mat[3][2] = r_3;
	pointsCoordinates.mat[4][0] = r_3;	pointsCoordinates.mat[4][1] = r_3;	pointsCoordinates.mat[4][2] = -r_3;
	pointsCoordinates.mat[5][0] = -r_3;	pointsCoordinates.mat[5][1] = r_3;	pointsCoordinates.mat[5][2] = -r_3;
	pointsCoordinates.mat[6][0] = r_3;	pointsCoordinates.mat[6][1] = -r_3;	pointsCoordinates.mat[6][2] = -r_3;
	pointsCoordinates.mat[7][0] = -r_3;	pointsCoordinates.mat[7][1] = -r_3;	pointsCoordinates.mat[7][2] = -r_3;


	double k = 1.0;
	weights.push_back(k);
	weights.push_back(k);
	weights.push_back(k);
	weights.push_back(k);
	weights.push_back(k);
	weights.push_back(k);
	weights.push_back(k);

}

/* ------------------------------------------------------------------------
Shape functions operations
---------------------------------------------------------------------------*/

void Operations::getGaussPoints(vector<vector<double>> &gaussPointsCoord, vector<vector<int>> &pointsIDPerElement, GetMesh mesh, int volIDField)
{
	Messages messages;
	messages.logMessage("Calculating Gauss points");
	Operations oper;
	ShapeFunctions shape;
	vector< vector<double>> thisGauss;

	// Loop for all the elements
	int thisElemType = 150;
	int pointCounter = 0;
	for (int i = 0; i < mesh.numElements; i++)
	{
		if (mesh.physicalTags[i] == volIDField)
		{

			// Only gets the Gauss points when 
			if (thisElemType != mesh.elemTypes[i])
			{
				thisElemType = mesh.elemTypes[i];
				GaussLegendrePoints thisElemGauss(thisElemType);
				thisGauss = thisElemGauss.pointsCoordinates.mat;

			}

			// Loop for all the Gauss points
			vector<int> thisPointsID;
			for (int eachPoint = 0; eachPoint < thisGauss.size(); eachPoint++)
			{

				vector<double> ptUV;
				ptUV = thisGauss[eachPoint];

				vector<double> ptXY = oper.scalLocalToReal(thisElemType, i, mesh, ptUV);
				gaussPointsCoord.push_back(ptXY);
				thisPointsID.push_back(pointCounter);
				pointCounter++;
			}
			pointsIDPerElement.push_back(thisPointsID);
		}
	}
	messages.logMessage("Calculating Gauss points: Done");

}



Matrix Operations::getCoordJac(int elemID, GetMesh mesh) {

	size_t numNodes = mesh.elemNodes[elemID].size();
	Matrix coordJac(numNodes, 3);

	for (int i = 0; i < numNodes; i++)
	{
		int thisNode = mesh.elemNodes[elemID][i];
		coordJac.mat[i][0] = mesh.nodesCoordinates[thisNode][0];
		coordJac.mat[i][1] = mesh.nodesCoordinates[thisNode][1];
		coordJac.mat[i][2] = mesh.nodesCoordinates[thisNode][2];
	}
	return coordJac;

}

vector<double>  Operations::scalLocalToReal(int elemType, int elemID, GetMesh mesh, vector<double> point) {

	ShapeFunctions shape;
	Operations oper;
	GmshNodesNumbering gmshNumbering;

	vector<int>order = gmshNumbering.GetGmshNodesNumbering(elemType);
	Matrix shapeFunction = shape.GetNodalShapeFunction(elemType, point[0], point[1], point[2]);
	//shapeFunction.print("shapeFunction 1");
	shapeFunction= gmshNumbering.ConvertToGmshNumbering(shapeFunction, order);
	//shapeFunction.print("shapeFunction 2");

	Matrix coordJac = oper.getCoordJac(elemID, mesh);
	Matrix xy = coordJac.mutiply(shapeFunction, coordJac);

	if (false)
	{
		cout << endl;
		cout << "=================================" << endl;
		cout << "Ploting scalLocalToReal" << endl;
		cout << "u=" << point[0] << " v=" << point[1] << " p=" << point[2] << endl;
		shapeFunction.print("shapeFunction");
		coordJac.print("coordJac");
		xy.print("xy");
	}

	vector<double> ans = { xy.mat[0][0],xy.mat[0][1],xy.mat[0][2]};

	return ans;
}

Matrix Operations::Jacobian(int elemType, int elemID, GetMesh mesh, vector<double> point) {

	ShapeFunctions shape;
	Operations oper;
	GmshNodesNumbering gmshNumbering;
	vector<int>order = gmshNumbering.GetGmshNodesNumbering(elemType);

	Matrix gradN = shape.GetGradNodalShapeFunction(elemType, point[0], point[1], point[2]);
	gradN = gmshNumbering.ConvertToGmshNumbering(gradN, order);

	Matrix coordJac = oper.getCoordJac(elemID, mesh);
	Matrix jac = jac.mutiply(gradN, coordJac);

	/*
	cout << "=================================" << endl;
	cout << "=================================" << endl;
	gradN.print("gradn");
	cout << "=================================" << endl;
	coordJac.print("coordjac");
	cout << "=================================" << endl;
	jac.print("jac");
	cout << "=================================" << endl;
	cout << "det: " << jac.Det_3x3();
	*/

	return jac;

}

vector<int> GmshNodesNumbering::GetGmshNodesNumbering(int ElemType)
{
	vector<int> ans;
	switch (ElemType)
	{
	case 1: //First order line
		ans = GmshLineFirstOrder();
		break;

	case 4: //First order tetrahedral
		ans = GmshTetrahedralFirstOrder();
		break;

	case 11://Second order tetrahedral
		ans = GmshTetrahedralSecondtOrder();
		break;

	default:
		msg.NotImplementedElement(ElemType, "GetGmshNodesNumbering");
	}
	return ans;

}

Matrix GmshNodesNumbering::ConvertToGmshNumbering(Matrix mat, vector<int> order)
{
	Matrix mat2(mat.rows, mat.cols);
	int col_counter = 0;
	for each (int col in order)
	{
		for (int row = 0; row < mat.rows; row++)
		{
			mat2.mat[row][col_counter] = mat.mat[row][col];
		}
		col_counter++;
	}
	return mat2;
}

vector<int> GmshNodesNumbering::GmshLineFirstOrder()
{
	vector<int> ans = { 0,1 };
	return ans;
}

vector<int> GmshNodesNumbering::GmshTetrahedralFirstOrder()
{
	vector<int> ans = { 0,2,3,1 };
	return ans;
}

vector<int> GmshNodesNumbering::GmshTetrahedralSecondtOrder()
{
	vector<int> ans = { 0,4,9,2,5,8,6,1,7,3 };
	return ans;

}