#include <string>
#include<vector>
using namespace std;
using std::string;

/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include "Matrix.h"
#include "Gmsh.h"

#ifndef _SHAPE_FUNCTIONS_INCLUDED_
#define _SHAPE_FUNCTIONS_INCLUDED_

class GmshNodesNumbering
{
public:
	vector<int> GetGmshNodesNumbering(int ElemType);
	Matrix ConvertToGmshNumbering(Matrix mat, vector<int> order);

private:
	vector<int> GmshLineFirstOrder();
	vector<int> GmshLineSecondOrder();

	vector<int> GmshTriangleFirstOrder();
	vector<int> GmshTriangleSecondOrder();

	vector<int> GmshTetrahedralFirstOrder();
	vector<int> GmshTetrahedralSecondtOrder();

};

class ShapeFunctions
{
public:
	/**
	Nodal shape functions

	@param ElemType: element type
	@param u,v,p: local coordinates
	@return 2D Matrix
	[N1 N2...Nn]
	*/
	Matrix GetNodalShapeFunction(int ElemType, double u, double v, double p);

	/**
	Gradient of nodal shape functions

	@param ElemType: element type
	@param u,v,p: local coordinates
	@return 2D Matrix
	[dN1/du		dN2/du	... dNn/du,
	dN1/dv		dN2/dv	... dNn/dv,
	dN1/dp		dN2/dp	... dNn/dp]
	*/
	Matrix GetGradNodalShapeFunction(int ElemType, double u = 0, double v = 0, double p = 0);

	~ShapeFunctions();

private:
	// Nodal shape functions
	Matrix  nodalLineFirstOrder(double u);
	Matrix  nodalLineSecondOrder(double u);

	Matrix  nodalTriangleFirstOrder(double u,double v);
	Matrix  nodalTriangleSecondOrder(double u, double v);


	Matrix  nodalTetrahedralFirstOrder(double u, double v, double p);
	Matrix  nodalTetrahedralSecondOrder(double u, double v, double p);

	Matrix  nodalHexahedralFirstOrder(double u, double v, double p);

	// Gradient of nodal shape functions
	Matrix  gradNodalLineFirstOrder();
	Matrix  gradNodalLineSecondOrder(double u);
	
	Matrix  gradNodalTriangleFirstOrder();
	Matrix  gradNodalTriangleSecondOrder(double u, double v);
	
	Matrix  gradNodalTetrahedralFirstOrder();
	Matrix  gradNodalTetrahedralSecondOrder(double u, double v, double p);

	Matrix  gradNodalHexahedralFirstOrder(double u, double v, double p);
};

class GaussLegendrePoints
{
public:
	/**
	Gauss points

	@param ElemType: element type
	*/
	GaussLegendrePoints(int ElemType);

	/**
	Local coordinates

	@return 2D Matrix
	[u1	v1 p1,
	u2	v2 p2,
	un	vn pn]
	*/
	Matrix pointsCoordinates;

	/**
	Weight for each point

	@return 2D Matrix
	[W1 w2 ... Wn]
	*/
	vector<double> weights;
	~GaussLegendrePoints();

private:
	void lineTwoPoints();
	void lineOnePoint();

	void triangleOnePointsInside();
	void triangleThreePointsInside();
	void triangleFourPointsInside();
	void triangleSevenPointsInside();

	void tetrahedralOnePoint();
	void tetrahedralFourPoinst();
	void tetrahedralFivePoinst();
	
	void hexahedralEightPoinstInside();

};


class Operations {
public:

	int getElemDimension(int ElemType);
	void getGaussPoints(vector<int> &elem_ID_list, vector<vector<double>> &gaussPointsCoord, vector<vector<int>> &pointsIDPerElement, GetMesh mesh, vector<int> volIDField);
	void getGaussPoints(vector<vector<double>> &gaussPointsCoord, vector<vector<int>> &pointsIDPerElement, GetMesh mesh, vector<int> volIDField);
	void getGaussPointsVol(vector<int> &elem_ID_list,vector<vector<double>> &gaussPointsCoord, vector<vector<int>> &pointsIDPerElement, GetMesh mesh, vector<int> volIDField);
	void getGaussPointsAdaptive(vector<int> &elem_ID_list, int &point_counter, vector<vector<double>> &gaussPointsCoord, vector<vector<int>> &pointsIDPerElement, GetMesh mesh, vector<int> GaussPointID,vector<int> forceFourPoints);

	//(GetMesh mesh, int phySurfaceFilter, vector<int> phyVolumeFilter, vector<vector<double>> &normalVectors, vector<double> &area
	
	/**
	Coordinates of the points of one element

	@param elemID: mesh ID of the element
	@mesh mesh information
	@return 2D Matrix
	[x1	y1 z1,
	x2	y2 z2,
	xn	yn zn]
	*/
	Matrix getCoordJac(int elemID, GetMesh mesh);

	/**
	Convert an scalar from local to global

	@param elemType: type of the element
	@param elemID: ID of the element in the mesh arrays
	@param mesh/ mesh data
	@param point: 2D array with the coordinates u,v,p
	@return 2D Matrix
	[x	y z]

	*/
	vector<double> scalLocalToReal(int elemType, int elemID, GetMesh mesh, vector<double> point);
	Matrix Jacobian(int elemType, int elemID, GetMesh mesh, vector<double> point);
	double getDetJac1D(Matrix mat);

private:
	void getGaussPoints_private(vector<int> &elem_ID_list,vector<vector<double>> &gaussPointsCoord, vector<vector<int>> &pointsIDPerElement, GetMesh mesh);


};


#endif