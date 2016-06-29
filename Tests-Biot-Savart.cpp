#include "stdafx.h"
#include <iostream>
#include<math.h>
using namespace std;
#include<vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include "Biot-Savart.h"
#include "ShapeFunctions.h"
#include "Matrix.h"
#include "Gmsh.h"
#include "Messages.h"
#include "Tests-Biot-Savart.h"
#include "AnalyticalSolutions.h"


void TestsBiotSavart::volumeDomain()
{

////////////////////////////////////////////////////////////////
//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\23_Biot_Savart_3D\\Solid winding\\");
	string meshFile_first("model_first.msh");
	string meshFile_second("model_second.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh_first(path + meshFile_first);
	GetMesh mesh_second(path + meshFile_second);

////////////////////////////////////////////////////////////////
//Points
	vector<vector<double>>pointsCoordinates;
	vector<vector<int>>pointsID;
	Operations points;
	points.getGaussPoints(pointsCoordinates, pointsID, mesh_first, 1002);

////////////////////////////////////////////////////////////////
// Analytical
	BiotSavartAnalyt analy;
	vector<double> analyPoints;
	vector<vector<double>> plotAnaly;
	for each (vector<double> thisPoint in pointsCoordinates)
	{
		analyPoints.push_back(thisPoint[2]);
	}
	
	vector<double> Hanaly = analy.getHcoil((45.0+60.0)/2000.0, analyPoints,0.5,1.0,100.0);
	int counter = 0;
	for each (double thisData in Hanaly)
	{
		vector<double> line;
		line.push_back(analyPoints[counter]);
		line.push_back(thisData);
		plotAnaly.push_back(line);
		counter++;
	}

////////////////////////////////////////////////////////////////
// Biot Savart inetgration
	BiotSavart biotSavart;
	double current_Density = 13333.33333;
	vector<double>centerPosition = { 0.0,0.0,0.0 };
	int volIDCurrente = 1000;
	int alongAxis = 0;
	vector<vector<double>> field_first = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_first, pointsCoordinates, path);
	vector<vector<double>> field_second = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_second, pointsCoordinates, path);

	vector<vector<double>> BiotSavart_first;
	vector<vector<double>> BiotSavart_second;


	counter = 0;
	for each (vector<double> thisData in field_first)
	{
		vector<double> line_first;
		line_first.push_back(pointsCoordinates[counter][0]);
		double mag = sqrt(pow(field_first[counter][0], 2) + pow(field_first[counter][2], 2));
		line_first.push_back(mag);

		vector<double> line_second;
		line_second.push_back(pointsCoordinates[counter][0]);
		mag = sqrt(pow(field_second[counter][0], 2) + pow(field_second[counter][2], 2));
		line_second.push_back(mag);

		BiotSavart_first.push_back(line_first);
		BiotSavart_second.push_back(line_second);

		counter++;
	}

////////////////////////////////////////////////////////////////
// Post processing
	PostProcessing post;
	post.writeVectorField(pointsCoordinates, field_first, "H", path + "\\results\\Gmsh_H_vector_first.txt");
	post.writeVectorField(pointsCoordinates, field_second, "H", path + "\\results\\Gmsh_H_vector_second.txt");

	post.writeDataResults(plotAnaly, path, "analyt");
	post.writeDataResults(BiotSavart_first, path, "Biot_first");
	post.writeDataResults(BiotSavart_second, path, "Biot_sec");
}

void TestsBiotSavart::currentLoop()
{

////////////////////////////////////////////////////////////////
//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\23_Biot_Savart_3D\\Line\\");
	string meshFile_first("model_first.msh");
	string meshFile_second("model_second.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh_first(path + meshFile_first);
	GetMesh mesh_second(path + meshFile_second);


////////////////////////////////////////////////////////////////
//Points
	vector<vector<double>>pointsCoordinatesZ;
	vector<vector<int>>pointsIDZ;
	vector<vector<double>>pointsCoordinatesX;
	vector<vector<int>>pointsIDX;
	Operations points;
	points.getGaussPoints(pointsCoordinatesZ, pointsIDZ, mesh_first, 1001);
	points.getGaussPoints(pointsCoordinatesX, pointsIDX, mesh_first, 1002);

////////////////////////////////////////////////////////////////
// Analytical
	BiotSavartAnalyt analy;
	vector<double> analyPoints;
	double current = 100.0;
	vector<vector<double>> plotAnaly;
	for each (vector<double> thisPoint in pointsCoordinatesZ)
	{
		analyPoints.push_back(thisPoint[2]);
	}

	vector<double> Hanaly = analy.getHCurrentLoop(150.0/1000.0,analyPoints, current);
	int counter = 0;
	for each (double thisData in Hanaly)
	{
		vector<double> line;
		line.push_back(analyPoints[counter]);
		line.push_back(thisData);
		plotAnaly.push_back(line);
		counter++;
	}

////////////////////////////////////////////////////////////////
// Biot Savart inetgration
	BiotSavart biotSavart;
	int volID = 1000;
	vector<vector<double>> field_first_Z = biotSavart.integrateLine(current,volID,mesh_first, pointsCoordinatesZ, path);
	vector<vector<double>> field_second_Z = biotSavart.integrateLine(current, volID, mesh_second, pointsCoordinatesZ, path);
	vector<vector<double>> field_first_X = biotSavart.integrateLine(current, volID, mesh_first, pointsCoordinatesX, path);
	vector<vector<double>> field_second_X = biotSavart.integrateLine(current, volID, mesh_second, pointsCoordinatesX, path);

	vector<vector<double>> BiotSavart_first_Z;
	vector<vector<double>> BiotSavart_second_Z;
	vector<vector<double>> BiotSavart_first_X;
	vector<vector<double>> BiotSavart_second_X;

// Line along Z
	counter = 0;
	for each (vector<double> thisData in pointsCoordinatesZ)
	{
		vector<double> line_first;
		line_first.push_back(pointsCoordinatesZ[counter][2]);
		line_first.push_back(abs(field_first_Z[counter][2]));

		vector<double> line_second;
		line_second.push_back(pointsCoordinatesZ[counter][2]);
		line_second.push_back(abs(field_second_Z[counter][2]));

		BiotSavart_first_Z.push_back(line_first);
		BiotSavart_second_Z.push_back(line_second);

		counter++;
	}

// Line along X
	counter = 0;
	for each (vector<double> thisData in pointsCoordinatesX)
	{
		vector<double> line_first;
		line_first.push_back(pointsCoordinatesX[counter][0]);
		double mag = sqrt(pow(field_first_X[counter][0], 2) + pow(field_first_X[counter][2], 2));
		line_first.push_back(mag);

		vector<double> line_second;
		line_second.push_back(pointsCoordinatesX[counter][0]);
		mag = sqrt(pow(field_second_X[counter][0], 2) + pow(field_second_X[counter][2], 2));

		line_second.push_back(mag);

		BiotSavart_first_X.push_back(line_first);
		BiotSavart_second_X.push_back(line_second);

		counter++;
	}


////////////////////////////////////////////////////////////////
// Post processing
	PostProcessing post;

	post.writeDataResults(plotAnaly, path, "analyt");
	post.writeDataResults(BiotSavart_first_Z, path, "Biot_first_z");
	post.writeDataResults(BiotSavart_second_Z, path, "Biot_second_z");
	post.writeDataResults(BiotSavart_first_X, path, "Biot_first_x");
	post.writeDataResults(BiotSavart_second_X, path, "Biot_second_x");
	
}

void TestsBiotSavart::TwoD()
{

////////////////////////////////////////////////////////////////
//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\23_Biot_Savart_3D\\2D\\BS\\");
	string meshFile_first("model_first.msh");
	string meshFile_second("model_second.msh");
	GetMesh mesh_first(path + meshFile_first);
	GetMesh mesh_second(path + meshFile_second);

////////////////////////////////////////////////////////////////
//Points
	vector<vector<double>>pointsCoordinatesX;
	vector<vector<int>>pointsIDX;
	vector<vector<double>>pointsCoordinatesY;
	vector<vector<int>>pointsIDY;
	Operations points;
	points.getGaussPoints(pointsCoordinatesX, pointsIDX, mesh_first, 1001);
	points.getGaussPoints(pointsCoordinatesY, pointsIDY, mesh_first, 1000);

	vector<vector<double>>teste = { {0,0,0} };
	////////////////////////////////////////////////////////////////
	// Biot Savart inetgration
	BiotSavart biotSavart;
	int volID = 1003;
	double current_Density = 6250000;
	int numberGaussPoints = teste.size();
	vector<vector<double>> Hresults(numberGaussPoints, vector<double>(3, 0));

	biotSavart.integrateTwoD(Hresults, volID, current_Density, mesh_second, teste, path);

}

