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
	vector<vector<double>>pointsCoordinatesZ;
	vector<vector<int>>pointsIDZ;
	vector<vector<double>>pointsCoordinatesX;
	vector<vector<int>>pointsIDX;
	Operations points;
	vector<int> volIDZ = {1001};
	vector<int> volIDX = { 1002 };
	points.getGaussPoints(pointsCoordinatesZ, pointsIDZ, mesh_first, volIDZ);
	points.getGaussPoints(pointsCoordinatesX, pointsIDX, mesh_first, volIDX);


////////////////////////////////////////////////////////////////
// Analytical
	BiotSavartAnalyt analy;
	vector<double> analyPoints;
	vector<vector<double>> plotAnaly;
	for each (vector<double> thisPoint in pointsCoordinatesZ)
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
	vector<vector<double>> field_first_along_z = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_first, pointsCoordinatesZ, path);
	vector<vector<double>> field_second_along_z = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_second, pointsCoordinatesZ, path);
	vector<vector<double>> field_first_along_x = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_first, pointsCoordinatesX, path);
	vector<vector<double>> field_second_along_x = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_second, pointsCoordinatesX, path);

	vector<vector<double>> BiotSavart_first_z;
	vector<vector<double>> BiotSavart_second_z;
	vector<vector<double>> BiotSavart_first_x;
	vector<vector<double>> BiotSavart_second_x;

// line along Z
	counter = 0;
	for each (vector<double> thisData in field_first_along_z)
	{
		vector<double> line_first;
		line_first.push_back(pointsCoordinatesZ[counter][2]);
		line_first.push_back(field_first_along_z[counter][2]);

		vector<double> line_second;
		line_second.push_back(pointsCoordinatesZ[counter][2]);
		line_second.push_back(field_second_along_z[counter][2]);

		BiotSavart_first_z.push_back(line_first);
		BiotSavart_second_z.push_back(line_second);

		counter++;
	}

	// line along X
	counter = 0;
	for each (vector<double> thisData in field_first_along_x)
	{
		vector<double> line_first;
		line_first.push_back(pointsCoordinatesX[counter][0]);
		double mag = sqrt(pow(field_first_along_x[counter][0], 2) + pow(field_first_along_x[counter][2], 2));
		line_first.push_back(mag);

		vector<double> line_second;
		line_second.push_back(pointsCoordinatesX[counter][0]);
		mag = sqrt(pow(field_second_along_x[counter][0], 2) + pow(field_second_along_x[counter][2], 2));
		line_second.push_back(mag);

		BiotSavart_first_x.push_back(line_first);
		BiotSavart_second_x.push_back(line_second);

		counter++;
	}
////////////////////////////////////////////////////////////////
// Post processing
	PostProcessing post;
		post.writeDataResults(plotAnaly, path, "analyt");
	post.writeDataResults(BiotSavart_first_x, path, "Biot_first_x");
	post.writeDataResults(BiotSavart_second_x, path, "Biot_sec_x");
	post.writeDataResults(BiotSavart_first_z, path, "Biot_first_z");
	post.writeDataResults(BiotSavart_second_z, path, "Biot_sec_z");
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
	vector<int> volIDz = { 1001 };
	vector<int> volIDx = { 1002 };
	points.getGaussPoints(pointsCoordinatesZ, pointsIDZ, mesh_first, volIDz);
	points.getGaussPoints(pointsCoordinatesX, pointsIDX, mesh_first, volIDx);

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
	vector<int> volIDz = { 1001 };
	vector<int> volIDy = { 1000 };
	points.getGaussPoints(pointsCoordinatesX, pointsIDX, mesh_first, volIDz);
	points.getGaussPoints(pointsCoordinatesY, pointsIDY, mesh_first, volIDy);

////////////////////////////////////////////////////////////////
// Biot Savart inetgration
	BiotSavart biotSavart;
	int volID1 = 1003;
	int volID2 = 1004;
	double current_Density = 13333.33333;

	// Line along Y 
	int numberGaussPointsY = pointsCoordinatesY.size();
	vector<vector<double>> HresultsYFirst(numberGaussPointsY, vector<double>(3, 0));
	biotSavart.integrateTwoD(HresultsYFirst, volID1, -current_Density, mesh_first, pointsCoordinatesY, path);
	biotSavart.integrateTwoD(HresultsYFirst, volID2, current_Density, mesh_first, pointsCoordinatesY, path);

	vector<vector<double>> HresultsYSecond(numberGaussPointsY, vector<double>(3, 0));
	biotSavart.integrateTwoD(HresultsYSecond, volID1, -current_Density, mesh_second, pointsCoordinatesY, path);
	biotSavart.integrateTwoD(HresultsYSecond, volID2, current_Density, mesh_second, pointsCoordinatesY, path);

	// Line along X 
	int numberGaussPointsX = pointsCoordinatesX.size();
	vector<vector<double>> HresultsXFirst(numberGaussPointsX, vector<double>(3, 0));
	biotSavart.integrateTwoD(HresultsXFirst, volID1, -current_Density, mesh_first, pointsCoordinatesX, path);
	biotSavart.integrateTwoD(HresultsXFirst, volID2, current_Density, mesh_first, pointsCoordinatesX, path);

	vector<vector<double>> HresultsXSecond(numberGaussPointsX, vector<double>(3, 0));
	biotSavart.integrateTwoD(HresultsXSecond, volID1, -current_Density, mesh_second, pointsCoordinatesX, path);
	biotSavart.integrateTwoD(HresultsXSecond, volID2, current_Density, mesh_second, pointsCoordinatesX, path);


////////////////////////////////////////////////////////////////
// Post processing

	vector<vector<double>> BiotSavartXFirst;
	vector<vector<double>> BiotSavartXSecond;
	vector<vector<double>> BiotSavartYFirst;
	vector<vector<double>> BiotSavartYSecond;


	int counter = 0;
	for each (vector<double> thisData in pointsCoordinatesY)
	{
		vector<double> lineYFirst;
		lineYFirst.push_back(pointsCoordinatesY[counter][1]);
		lineYFirst.push_back(abs(HresultsYFirst[counter][0]));
		BiotSavartYFirst.push_back(lineYFirst);

		vector<double> lineYSecond;
		lineYSecond.push_back(pointsCoordinatesY[counter][1]);
		lineYSecond.push_back(abs(HresultsYSecond[counter][0]));
		BiotSavartYSecond.push_back(lineYFirst);
		
		counter++;
	}


	counter = 0;
	for each (vector<double> thisData in pointsCoordinatesX)
	{
		vector<double> lineXFirst;
		lineXFirst.push_back(pointsCoordinatesX[counter][0]);
		lineXFirst.push_back(HresultsXFirst[counter][0]);
		BiotSavartXFirst.push_back(lineXFirst);

		vector<double> lineXSecond;
		lineXSecond.push_back(pointsCoordinatesX[counter][0]);
		lineXSecond.push_back(HresultsXSecond[counter][0
		]);
		BiotSavartXSecond.push_back(lineXSecond);

		counter++;
	}
	PostProcessing post;
	post.writeDataResults(BiotSavartXFirst, path, "BiotSavartXFirst");
	post.writeDataResults(BiotSavartXSecond, path, "BiotSavartXSecond");
	post.writeDataResults(BiotSavartYFirst, path, "BiotSavartYFirst");
	post.writeDataResults(BiotSavartYSecond, path, "BiotSavartYSecond");


}

void TestsBiotSavart::Rele3D()
{

	////////////////////////////////////////////////////////////////
	//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\25_Rele_3D\\03_FFFM_complete\\");
	string meshFile("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh(path + meshFile);

	////////////////////////////////////////////////////////////////
	//Points
	vector<vector<double>>pointsCoordinates;
	vector<vector<int>>pointsID;
	Operations points;
	vector<int>volIDy = {2000,2001,2002,2003,2004,2005};
	points.getGaussPoints(pointsCoordinates, pointsID, mesh, volIDy);


	////////////////////////////////////////////////////////////////
	// Biot Savart inetgration
	BiotSavart biotSavart;
	double current_Density = 41666666.67;
	vector<double>centerPosition = { 0.1,0.0,0.0 };
	int volIDCurrente = 1000;
	int alongAxis = 0;
	/*vector<vector<double>> field_first = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_first, pointsCoordinates, path);
	vector<vector<double>> field_second = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_second, pointsCoordinates, path);
*/
	//vector<vector<double>> BiotSavart_first;
	//vector<vector<double>> BiotSavart_second;


	//counter = 0;
	//for each (vector<double> thisData in field_first)
	//{
	//	vector<double> line_first;
	//	line_first.push_back(pointsCoordinates[counter][0]);
	//	double mag = sqrt(pow(field_first[counter][0], 2) + pow(field_first[counter][2], 2));
	//	line_first.push_back(mag);

	//	vector<double> line_second;
	//	line_second.push_back(pointsCoordinates[counter][0]);
	//	mag = sqrt(pow(field_second[counter][0], 2) + pow(field_second[counter][2], 2));
	//	line_second.push_back(mag);

	//	BiotSavart_first.push_back(line_first);
	//	BiotSavart_second.push_back(line_second);

	//	counter++;
	//}

	//////////////////////////////////////////////////////////////////
	//// Post processing
	//PostProcessing post;
	//post.writeVectorField(pointsCoordinates, field_first, "H", path + "\\results\\Gmsh_H_vector_first.txt");
	//post.writeVectorField(pointsCoordinates, field_second, "H", path + "\\results\\Gmsh_H_vector_second.txt");

	//post.writeDataResults(plotAnaly, path, "analyt");
	//post.writeDataResults(BiotSavart_first, path, "Biot_first");
	//post.writeDataResults(BiotSavart_second, path, "Biot_sec");



}

