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
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\25_Rele_3D\\06_Source_Field\\");
	string meshFile_points("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh_points(path + meshFile_points);

	string meshFile_integration("model_wind_seg.msh");
	GetMesh mesh_integration(path + meshFile_integration);

	////////////////////////////////////////////////////////////////
	//Points
	vector<vector<double>>pointsCoordinates;
	vector<vector<int>>pointsID;
	Operations points;
	//vector<int>volID = {2000,2001,2002,2003,2004,2005};
	vector<int>volID = { 2000,2001,2002};

	vector<int> Elem_IDs;
	points.getGaussPointsVol(Elem_IDs,pointsCoordinates, pointsID, mesh_points, volID);

		
	////////////////////////////////////////////////////////////////
	// Biot Savart inetgration
	BiotSavart biotSavart;
	double current_Density = pow(10,7);
	vector<double>centerPosition = { 0.01,0.00,0.01 };
	int volIDCurrente = 2005;
	int alongAxis = 0;
	vector<vector<double>> field = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_integration, pointsCoordinates, path);

	//////////////////////////////////////////////////////////////////
	//// Post processing
	PostProcessing post;


	//writes the field solution to Gmsh
	post.writeVectorField(pointsCoordinates, field, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	
	//Writes the Gauss points information
	//IDs_Gauss_Points
	//Coordinates_Gauss_Points

	post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);

	//Writes the magnetic field
	post.writeDataResults(field, path, "Gauss_Points_H_field");


	//////////////////////////////////////////////////////////////////
	//vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	//// Biot Savart inetgration
	//vector<vector<double>> field = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_integration, points_list, path);

	//

	//////Writes the magnetic field
	//post.writeDataResults(field, path, "Gauss_Points_H_field_line");

	//////writes the field solution to Gmsh
	//post.writeVectorField(points_list, field, "H", path + "\\results\\Gauss_Points_H_field_line_Gmsh.txt");



}

void TestsBiotSavart::Teste3D()
{

	////////////////////////////////////////////////////////////////
	//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\27_Teste_Campo_Fonte_3D\\");
	string meshFile("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh(path + meshFile);

	////////////////////////////////////////////////////////////////
	//Points
	vector<vector<double>>pointsCoordinates;
	vector<vector<int>>pointsID;
	Operations points;
	vector<int>volID = { 1000 };
	points.getGaussPoints(pointsCoordinates, pointsID, mesh, volID);

	////////////////////////////////////////////////////////////////
	// Biot Savart inetgration
	BiotSavart biotSavart;
	double current= 1;
	int volID_int =  1001 ;
	vector<vector<double>> field = biotSavart.integrateLine(current, volID_int,mesh,pointsCoordinates, path);

	//////////////////////////////////////////////////////////////////
	//// Post processing
	PostProcessing post;

	//writes the field solution to Gmsh
	post.writeVectorField(pointsCoordinates, field, "H", path + "\\results\\\Gauss_Points_H_field_Gmsh.txt");

	//Writes the Gauss points information
	//IDs_Gauss_Points
	//Coordinates_Gauss_Points

	//post.writeGaussPointsIDs(pointsID, pointsCoordinates, path);

	//Writes the magnetic field
	post.writeDataResults(field, path, "\Gauss_Points_H_field.txt");
}

void TestsBiotSavart::Teste3D_core()
{

	////////////////////////////////////////////////////////////////
	//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\27_Teste_Campo_Fonte_3D\\Dispositivo_simples\\");
	string meshFile("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh(path + meshFile);

	////////////////////////////////////////////////////////////////
	//Points
	vector<vector<double>>pointsCoordinates;
	vector<vector<int>>pointsID;
	Operations points;
	vector<int>volID = {1000,1001};
	points.getGaussPoints(pointsCoordinates, pointsID, mesh, volID);



	////////////////////////////////////////////////////////////////
	// Biot Savart inetgration
	BiotSavart biotSavart;
	double current_Density = 41666666.67;
	vector<double>centerPosition = { 0.00,0.00,0.00 };
	int volIDCurrente = 1001;
	int alongAxis = 0;
	////vector<vector<double>> field = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);

	//PM
	vector<vector<double>> field;
	

	//////////////////////////////////////////////////////////////////
	//// Post processing
	PostProcessing post;

	//writes the field solution to Gmsh
	post.writeVectorField(pointsCoordinates, field, "H", path + "\\results\\Gmsh_H_vector_first.txt");

	//Writes the Gauss points information
	//IDs_Gauss_Points
	//Coordinates_Gauss_Points

	//post.writeGaussPointsIDs(pointsID, pointsCoordinates, path);

	//Writes the magnetic field
	post.writeDataResults(field, path, "Hfield_Gauss_Points");
}

void TestsBiotSavart::Subdomain()
{

	//*******************************
	// Mesh layer
	//*******************************

	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\31_Subdomain_Dular_2009\\Subdomain\\");
	string meshFile("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh(path + meshFile);

	int surf = 177;
	vector<int> vol = { 174,173,175,176 };
	//vector<int> vol = { 174};

	vector<vector<double>> normalVectors;
	vector<int> Elem_IDs;
	//vector<int> Elem_IDs = mesh.defineBoundary(mesh, surf, vol, 1, normalVectors);
	//mesh.writeMesh(mesh, path + "\\mesh_layer.msh", Elem_IDs);
	//////////////////////////////////////////////////////////////
	//Gauss Points
	vector<vector<double>>pointsCoordinates;
	vector<vector<int>>pointsID;

	Operations points;

	points.getGaussPointsVol(Elem_IDs,pointsCoordinates, pointsID, mesh, vol);

	////////////////////////////////////////////////////////////////
	// Biot Savart inetgration
	BiotSavart biotSavart;
	double current_Density = pow(10, 7)*1000.0;
	vector<double>centerPosition = { 0.01,0.00,0.01 };
	int volIDCurrente = 173;
	int alongAxis = 0;
	vector<vector<double>> field1 = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	vector<vector<double>> field2 = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);

	//sum the two fields
	int rowCounter = 0;
	for each (vector<double> line in field1)
	{
		int colCounter = 0;
		for each (double col in line)
		{
			field1[rowCounter][colCounter] += field2[rowCounter][colCounter];
			colCounter ++;
		}
		rowCounter ++;
	}

	
	//////////////////////////////////////////////////////////////////
	//// Post processing
	PostProcessing post;


	//writes the field solution to Gmsh
	post.writeVectorField(pointsCoordinates, field1, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");


	//Writes the Gauss points information
	//IDs_Gauss_Points
	//Coordinates_Gauss_Points

	post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);

	//Writes the magnetic field
	post.writeDataResults(field1, path, "Gauss_Points_H_field");

//	
	//////////////////////////////////////////////////////////////////
	vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	// Biot Savart inetgration
	vector<vector<double>> field3 = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, points_list, path);
	vector<vector<double>> field4 = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, points_list, path);

	//sum the two fields
	rowCounter = 0;
	for each (vector<double> line in field3)
	{
		int colCounter = 0;
		for each (double col in line)
		{
			field3[rowCounter][colCounter] += field4[rowCounter][colCounter];
			colCounter++;
		}
		rowCounter++;
	}
	
	//Writes the magnetic field
	post.writeDataResults(field3, path, "Gauss_Points_H_field_line");

	//writes the field solution to Gmsh
	post.writeVectorField(points_list, field3, "H", path + "\\results\\Gauss_Points_H_field_line_Gmsh.txt");
}


////get component of the field
//vector<vector<double>> tangent;
//vector<vector<double>> normal;
//vector<vector<double>> newField;
//
//
//Vector1D thisMath;
//int counter = 0;
//for each (vector<double> line in field1)
//{
//
//	vector<double> thisComponent1 = thisMath.crossProduct(normalVectors[counter], line);
//	vector<double> thisTangent1 = thisMath.crossProduct(normalVectors[counter], thisComponent1);
//	vector<double> thisTangent = thisMath.multiScal(thisTangent1, -1.0);
//	vector<double> thisNormal = thisMath.subtract(line, thisTangent);
//
//
//	tangent.push_back(thisTangent);
//	normal.push_back(thisNormal);
//	counter++;
//}
////////////////////////////////////////////////////////////////////
////// Post processing
//PostProcessing post;
//
//
////writes the field solution to Gmsh
//post.writeVectorField(pointsCoordinates, field1, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
//post.writeVectorField(pointsCoordinates, tangent, "Ht", path + "\\results\\Gauss_Points_H_tangent_field_Gmsh.txt");
//post.writeVectorField(pointsCoordinates, normal, "Hn", path + "\\results\\Gauss_Points_H_normal_field_Gmsh.txt");