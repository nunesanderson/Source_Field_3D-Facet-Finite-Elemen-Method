#include "stdafx.h"
#include <iostream>
#include<math.h>
using namespace std;
#include<vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <numeric>
#include <list>
#include <vector>
#include <numeric>
#include <map>
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
	vector<int> Elem_IDsz;
	vector<int> Elem_IDsx;

	points.getGaussPointsVol(Elem_IDsz, pointsCoordinatesZ, pointsIDZ, mesh_first, volIDZ);
	points.getGaussPointsVol(Elem_IDsx, pointsCoordinatesX, pointsIDX, mesh_first, volIDX);


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
	int alongAxis = 2;
	vector<vector<double>> field_first_along_z = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_first, pointsCoordinatesZ, path);
	vector<vector<double>> field_second_along_z = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_second, pointsCoordinatesZ, path);
	/*vector<vector<double>> field_first_along_x = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_first, pointsCoordinatesX, path);
	vector<vector<double>> field_second_along_x = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_second, pointsCoordinatesX, path);
*/
//	vector<vector<double>> BiotSavart_first_z;
//	vector<vector<double>> BiotSavart_second_z;
//	vector<vector<double>> BiotSavart_first_x;
//	vector<vector<double>> BiotSavart_second_x;
//
//// line along Z
//	counter = 0;
//	for each (vector<double> thisData in field_first_along_z)
//	{
//		vector<double> line_first;
//		line_first.push_back(pointsCoordinatesZ[counter][2]);
//		line_first.push_back(field_first_along_z[counter][2]);
//
//		vector<double> line_second;
//		line_second.push_back(pointsCoordinatesZ[counter][2]);
//		line_second.push_back(field_second_along_z[counter][2]);
//
//		BiotSavart_first_z.push_back(line_first);
//		BiotSavart_second_z.push_back(line_second);
//
//		counter++;
//	}
//
//	// line along X
//	counter = 0;
//	for each (vector<double> thisData in field_first_along_x)
//	{
//		vector<double> line_first;
//		line_first.push_back(pointsCoordinatesX[counter][0]);
//		double mag = sqrt(pow(field_first_along_x[counter][0], 2) + pow(field_first_along_x[counter][2], 2));
//		mag = field_first_along_x[counter][2];
//		line_first.push_back(mag);
//
//		vector<double> line_second;
//		line_second.push_back(pointsCoordinatesX[counter][0]);
//		mag = sqrt(pow(field_second_along_x[counter][0], 2) + pow(field_second_along_x[counter][2], 2));
//		mag = field_second_along_x[counter][2];
//		line_second.push_back(mag);
//
//		BiotSavart_first_x.push_back(line_first);
//		BiotSavart_second_x.push_back(line_second);
//
//		counter++;
//	}
//////////////////////////////////////////////////////////////////
//// Post processing
//	PostProcessing post;
//		post.writeDataResults(plotAnaly, path, "analyt");
//	post.writeDataResults(BiotSavart_first_x, path, "Biot_first_x");
//	post.writeDataResults(BiotSavart_second_x, path, "Biot_sec_x");
//	post.writeDataResults(BiotSavart_first_z, path, "Biot_first_z");
//	post.writeDataResults(BiotSavart_second_z, path, "Biot_sec_z");
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
	vector<int> Elem_IDsz;
	vector<int> Elem_IDsx;

	points.getGaussPointsVol(Elem_IDsz,pointsCoordinatesZ, pointsIDZ, mesh_first, volIDz);
	points.getGaussPointsVol(Elem_IDsx,pointsCoordinatesX, pointsIDX, mesh_first, volIDx);

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
		mag = field_first_X[counter][2];
		line_first.push_back(mag);

		vector<double> line_second;
		line_second.push_back(pointsCoordinatesX[counter][0]);
		mag = sqrt(pow(field_second_X[counter][0], 2) + pow(field_second_X[counter][2], 2));
		mag = field_second_X[counter][2];

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
	vector<int> Elem_IDsx;
	vector<int> Elem_IDsy;

	points.getGaussPointsVol(Elem_IDsx,pointsCoordinatesX, pointsIDX, mesh_first, volIDz);
	points.getGaussPointsVol(Elem_IDsy,pointsCoordinatesY, pointsIDY, mesh_first, volIDy);

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
		lineYFirst.push_back((HresultsYFirst[counter][0]));
		BiotSavartYFirst.push_back(lineYFirst);

		vector<double> lineYSecond;
		lineYSecond.push_back(pointsCoordinatesY[counter][1]);
		lineYSecond.push_back((HresultsYSecond[counter][0]));
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
		lineXSecond.push_back(HresultsXSecond[counter][0]);
		BiotSavartXSecond.push_back(lineXSecond);

		counter++;
	}
	PostProcessing post;
	post.writeDataResults(BiotSavartXFirst, path, "BiotSavartXFirst");
	post.writeDataResults(BiotSavartXSecond, path, "BiotSavartXSecond");
	post.writeDataResults(BiotSavartYFirst, path, "BiotSavartYFirst");
	post.writeDataResults(BiotSavartYSecond, path, "BiotSavartYSecond");


}

void TestsBiotSavart::Atuador3D()
{
	//*********************************************************************************************
	// Total volume
	//*********************************************************************************************

	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\34_Atuador_vert\\FFEM_Complete\\");
	string meshFile("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh(path + meshFile);
	PostProcessing post;
	Vector1D thisMath;

	vector<int>vol = { 2001,2002,2003,2004,2005 };
	vector<int> Elem_IDs;
	vector<vector<double>>pointsCoordinates;
	vector<vector<int>>pointsID;
	Operations points;
	points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh, vol);

	//-------------------------------------------------------------------------------------
	// Biot Savart integration
	BiotSavart biotSavart;
	double current_Density = 3.0*pow(10, 7);
	vector<double>centerPosition_1 = { -0.04,0.00,0.00 };
	vector<double>centerPosition_2 = { 0.04,0.00,0.00 };
	int alongAxis = 2;
	vector<vector<double>> field_1 = biotSavart.integrateSolidWinding(2003, current_Density, centerPosition_1, alongAxis, mesh, pointsCoordinates, path);
	vector<vector<double>> field_2 = biotSavart.integrateSolidWinding(2004, -current_Density, centerPosition_2, alongAxis, mesh, pointsCoordinates, path);
	vector<vector<double>> fieldTotal = thisMath.sum(field_1, field_2);
	//-------------------------------------------------------------------------------------
	 //Post processing
	/*post.writeVectorField(pointsCoordinates, fieldTotal, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);
	post.writeDataResults(fieldTotal, path, "Gauss_Points_H_field");
*/

	////*********************************************************************************************
	//// VS - only
	////*********************************************************************************************
	//string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\34_Atuador_vert\\FFEM_VS\\");
	//string meshFile("model.msh");
	//string HFieldFile("results\\H_Field.txt");
	//GetMesh mesh(path + meshFile);
	//PostProcessing post;
	//	Vector1D thisMath;

	//vector<int>vol = {2002,2005 };
	//vector<int> Elem_IDs;

	//vector<vector<double>>pointsCoordinates;
	//vector<vector<int>>pointsID;
	//Operations points;
	//points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh, vol);

	////-------------------------------------------------------------------------------------
	//// Biot Savart integration
	//BiotSavart biotSavart;
	//double current_Density = 3.0*pow(10, 7);
	//vector<double>centerPosition_1 = { -0.04,0.00,0.00 };
	//vector<double>centerPosition_2 = { 0.04,0.00,0.00 };
	//int alongAxis = 2;
	//vector<vector<double>> field_1 = biotSavart.integrateSolidWinding(2003, current_Density, centerPosition_1, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> field_2 = biotSavart.integrateSolidWinding(2004, -current_Density, centerPosition_2, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> fieldTotal = thisMath.sum(field_1, field_2);
	////-------------------------------------------------------------------------------------
	//// Post processing
	//post.writeVectorField(pointsCoordinates, fieldTotal, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	//post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);
	//post.writeDataResults(fieldTotal, path, "Gauss_Points_H_field");
	//
	////-------------------------------------------------------------------------------------
	//// Line
	//vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	//// Biot Savart inetgration
	//vector<vector<double>> field1Line = biotSavart.integrateSolidWinding(2003, current_Density, centerPosition_1, alongAxis, mesh, points_list, path);
	//vector<vector<double>> field2Line = biotSavart.integrateSolidWinding(2004, -current_Density, centerPosition_2, alongAxis, mesh, points_list, path);
	//vector<vector<double>> field1LineTotal = thisMath.sum(field1Line, field2Line);

	////Writes the magnetic field
	//post.writeDataResults(field1LineTotal, path, "line_field_BS");

	////writes the field solution to Gmsh
	//post.writeVectorField(points_list, field1LineTotal, "H", path + "\\results\\line_field_BS_Gmsh.txt");

	//*********************************************************************************************
	// VS corrected
	//*********************************************************************************************
	//string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\34_Atuador_vert\\FFEM_VS_corrected\\");
	//string meshFile("model.msh");
	//string HFieldFile("results\\H_Field.txt");
	//GetMesh mesh(path + meshFile);
	//PostProcessing post;
	//BiotSavart biotSavart;
	//Vector1D thisMath;

	//vector<int> vol = {2002,2005};
	//vector<int> Elem_IDs;

	//vector<vector<double>>pointsCoordinates;
	//vector<vector<int>>pointsID;
	//Operations points;
	//points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh, vol);

	////-------------------------------------------------------------------------------------
	//// Biot Savart integration
	//double current_Density = 3.0*pow(10, 7);
	//vector<double>centerPosition_1 = { -0.04,0.00,0.00 };
	//vector<double>centerPosition_2 = { 0.04,0.00,0.00 };
	//int alongAxis = 2;
	//vector<vector<double>> field_1 = biotSavart.integrateSolidWinding(2003, current_Density, centerPosition_1, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> field_2 = biotSavart.integrateSolidWinding(2004, -current_Density, centerPosition_2, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> fieldTotal = thisMath.sum(field_1, field_2);
	////-------------------------------------------------------------------------------------
	////Post processing
	//post.writeVectorField(pointsCoordinates, fieldTotal, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	//post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);
	//post.writeDataResults(fieldTotal, path, "Gauss_Points_H_field");

	////-------------------------------------------------------------------------------------
	//// Line
	//vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	//// Biot Savart inetgration
	//vector<vector<double>> field1Line = biotSavart.integrateSolidWinding(2003, current_Density, centerPosition_1, alongAxis, mesh, points_list, path);
	//vector<vector<double>> field2Line = biotSavart.integrateSolidWinding(2004, -current_Density, centerPosition_2, alongAxis, mesh, points_list, path);
	//vector<vector<double>> field1LineTotal = thisMath.sum(field1Line, field2Line);


	////Writes the magnetic field
	//post.writeDataResults(field1LineTotal, path, "line_field_BS");

	////writes the field solution to Gmsh
	//post.writeVectorField(points_list, field1LineTotal, "H", path + "\\results\\line_field_BS_Gmsh.txt");
	////*********************************************************************************************
	////BS correction
	////*********************************************************************************************
	//vector<int> surf_ID = { 2006 };
	//vector<int> Elem_IDsSF;

	//vector<vector<double>>pointsCoordinatesSF;
	//vector<vector<int>>pointsIDSF;
	//points.getGaussPointsVol(Elem_IDsSF, pointsCoordinatesSF, pointsIDSF, mesh, surf_ID);

	////-------------------------------------------------------------------------------------
	////External surface

	//// Biot Savart inetgration
	//vector<vector<double>> field1LineBS = biotSavart.integrateSolidWinding(2003, current_Density, centerPosition_1, alongAxis, mesh, pointsCoordinatesSF, path);
	//vector<vector<double>> field2LineBS = biotSavart.integrateSolidWinding(2004, -current_Density, centerPosition_2, alongAxis, mesh, pointsCoordinatesSF, path);
	//vector<vector<double>> field1LineTotalBS = thisMath.sum(field1LineBS, field2LineBS);


	////get component of the field
	//vector<vector<double>> tangent;
	//vector<vector<double>> normal;
	//vector<double> areas;

	//post.getFieldComponents(normal, tangent, areas, field1LineTotalBS, pointsCoordinatesSF, Elem_IDsSF, pointsIDSF, "Gauss_Points_H_field_surface", path, mesh);
	//vector<int>volID = {2001,2002,2003,2004,2005,2006};


	//points.getGaussPointsVol(Elem_IDs,pointsCoordinates, pointsID, mesh, volID);

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

void TestsBiotSavart::H_discon()
{
	//*********************************************************************************************
	// Convert the mesh
	//*********************************************************************************************
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\39_H_discon\\FFEM_Complete\\");
	string meshFileTemp("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh meshTemp(path + meshFileTemp);
	PostProcessing post;
	BiotSavart biotSavart;
	Vector1D thisMath;

	// Layer of elements - Inner
	vector<int> vol = {60,61};
	int sufaceIDTemp = 62;
	vector<vector<double>>normalVectors;
	vector<double>areasTemp;
	vector<int> elemIDLayerTemp;


	//*********************************************************************************************
	// Fiels calculation
	//*********************************************************************************************
	string meshFile("model.msh");
	GetMesh mesh(path + meshFile);
	vector<int> elemIDLayer;
	vector<vector<int>>pointsID;
	vector<vector<double>>pointsCoordinates;
	Operations points;
	elemIDLayer = mesh.defineVolumeBoundary(mesh, sufaceIDTemp, vol, 2, normalVectors, areasTemp, path);
	points.getGaussPointsVol(elemIDLayer, pointsCoordinates, pointsID, mesh, vol);
	meshTemp.writeMesh(mesh, path + "\\layer_mesh.msh", elemIDLayer, vol);

	//-------------------------------------------------------------------------------------
	//Field components along the element surface
	int counter = 0;
	vector<vector<double>>tangentlField;

	for each (int elem in elemIDLayer)
	{
		vector<double>this_tangentField = { 0,0,0 };


		for (size_t i = 0; i < 4; i++)
		{
			if (pointsCoordinates[counter][0]>0)
			{

			this_tangentField = { 0,-1,0 };
			}
			else
			{
				this_tangentField = { 0,1,0 };

			}
			
			tangentlField.push_back(this_tangentField);
			counter++;
		}
	}


	
	//////-------------------------------------------------------------------------------------
	//////Post processing
	post.writeVectorField(pointsCoordinates, tangentlField, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	post.writeGaussPointsIDs(elemIDLayer, pointsID, pointsCoordinates, path);
	post.writeDataResults(tangentlField, path, "Gauss_Points_H_field");
	post.writeVectorField(pointsCoordinates, normalVectors, "H", path + "\\results\\Gauss_Points_H_field_normal_Gmsh.txt");
	post.writeVectorField(pointsCoordinates, tangentlField, "H", path + "\\results\\Gauss_Points_H_field_tangential_Gmsh.txt");

	////-------------------------------------------------------------------------------------
	//// Line
	////vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	////// Biot Savart inetgration
	////vector<vector<double>> field1Line = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, points_list, path);
	////vector<vector<double>> field2Line = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, points_list, path);
	////vector<vector<double>> field1LineTotal = thisMath.sum(field1Line, field2Line);


	//////Writes the magnetic field
	////post.writeDataResults(field1LineTotal, path, "line_field_BS");

	//////writes the field solution to Gmsh
	////post.writeVectorField(points_list, field1LineTotal, "H", path + "\\results\\line_field_BS_Gmsh.txt");

	////*********************************************************************************************
	////BS correction
	////*********************************************************************************************
	//vector<int> surf_ID = { 178,177 };
	//vector<int> Elem_IDsSF;

	//vector<vector<double>>pointsCoordinatesSF;
	//vector<vector<int>>pointsIDSF;
	//points.getGaussPointsVol(Elem_IDsSF, pointsCoordinatesSF, pointsIDSF, mesh, surf_ID);

	////-------------------------------------------------------------------------------------
	////External surface

	//// Biot Savart inetgration
	//vector<vector<double>> field1LineBS = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinatesSF, path);
	//vector<vector<double>> field2LineBS = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinatesSF, path);
	//vector<vector<double>> field1LineTotalBS = thisMath.sum(field1LineBS, field2LineBS);


	////get component of the field
	//vector<vector<double>> tangent;
	//vector<vector<double>> normal;
	//vector<double>areasSF;

	//post.getFieldComponents(normal, tangent, areasSF, field1LineTotalBS, pointsCoordinatesSF, Elem_IDsSF, pointsIDSF, "Gauss_Points_H_field_surface", path, mesh);

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

	//*********************************************************************************************
	// Total volume
	//*********************************************************************************************

	//string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\31_Subdomain_Dular_2009\\FFEM_Complete\\");
	//string meshFile("model.msh");
	//string HFieldFile("results\\H_Field.txt");
	//GetMesh mesh(path + meshFile);
	//PostProcessing post;
	//BiotSavart biotSavart;
	//Vector1D thisMath;

	//vector<int> vol = { 174,173,175,176};
	//vector<vector<double>> normalVectors;
	//vector<int> Elem_IDs;
	//

	//vector<vector<double>>pointsCoordinates;
	//vector<vector<int>>pointsID;
	//Operations points;
	//points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh, vol);

	////-------------------------------------------------------------------------------------
	//// Biot Savart integration
	//double current_Density = pow(10, 7)*1000.0;
	//vector<double>centerPosition = { 0.01,0.00,0.01 };
	//int volIDCurrente = 173;
	//int alongAxis = 0;
	//vector<vector<double>> field1 = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> field2 = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> fieldTotal = thisMath.sum(field1, field2);
	////-------------------------------------------------------------------------------------
	//// Post processing
	//post.writeVectorField(pointsCoordinates, fieldTotal, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	//post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);
	//post.writeDataResults(fieldTotal, path, "Gauss_Points_H_field");
	

	////*********************************************************************************************
	//// VS - only
	////*********************************************************************************************
	//string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\31_Subdomain_Dular_2009\\FFEM_VS\\");
	//string meshFile("model.msh");
	//string HFieldFile("results\\H_Field.txt");
	//GetMesh mesh(path + meshFile);
	//PostProcessing post;
	//BiotSavart biotSavart;
	//Vector1D thisMath;

	//vector<int> vol = { 174 };
	//vector<int> Elem_IDs;

	//vector<vector<double>>pointsCoordinates;
	//vector<vector<int>>pointsID;
	//Operations points;
	//points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh, vol);

	////-------------------------------------------------------------------------------------
	//// Biot Savart integration
	//double current_Density = pow(10, 7)*1000.0;
	//vector<double>centerPosition = { 0.01,0.00,0.01 };
	//int volIDCurrente = 173;
	//int alongAxis = 0;
	//vector<vector<double>> field1 = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> field2 = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> fieldTotal = thisMath.sum(field1, field2);
	////-------------------------------------------------------------------------------------
	//// Post processing
	//post.writeVectorField(pointsCoordinates, fieldTotal, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	//post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);
	//post.writeDataResults(fieldTotal, path, "Gauss_Points_H_field");
	//
	////-------------------------------------------------------------------------------------
	//// Line
	//vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	//// Biot Savart inetgration
	//vector<vector<double>> field1Line = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, points_list, path);
	//vector<vector<double>> field2Line = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, points_list, path);
	//vector<vector<double>> field1LineTotal = thisMath.sum(field1Line, field2Line);

	////Writes the magnetic field
	//post.writeDataResults(field1LineTotal, path, "line_field_BS");

	////writes the field solution to Gmsh
	//post.writeVectorField(points_list, field1LineTotal, "H", path + "\\results\\line_field_BS_Gmsh.txt");

	//*********************************************************************************************
	// VS corrected
	//*********************************************************************************************
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\31_Subdomain_Dular_2009\\FFEM_VS_corrected\\");
	string meshFile("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh(path + meshFile);
	PostProcessing post;
	BiotSavart biotSavart;
	Vector1D thisMath;

	vector<int> vol = { 174 };
	vector<int> Elem_IDs;

	vector<vector<double>>pointsCoordinates;
	vector<vector<int>>pointsID;
	Operations points;
	points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh, vol);

	//-------------------------------------------------------------------------------------
	// Biot Savart integration
	double current_Density = pow(10, 7)*1000.0;
	vector<double>centerPosition = { 0.01,0.00,0.01 };
	int volIDCurrente = 173;
	int alongAxis = 0;
	vector<vector<double>> field1 = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	vector<vector<double>> field2 = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	vector<vector<double>> fieldTotal = thisMath.sum(field1, field2);
	//-------------------------------------------------------------------------------------
	 //Post processing
	post.writeVectorField(pointsCoordinates, fieldTotal, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);
	post.writeDataResults(fieldTotal, path, "Gauss_Points_H_field");

	//-------------------------------------------------------------------------------------
	// Line
	vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	// Biot Savart inetgration
	vector<vector<double>> field1Line = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, points_list, path);
	vector<vector<double>> field2Line = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, points_list, path);
	vector<vector<double>> field1LineTotal = thisMath.sum(field1Line, field2Line);


	//Writes the magnetic field
	post.writeDataResults(field1LineTotal, path, "line_field_BS");

	//writes the field solution to Gmsh
	post.writeVectorField(points_list, field1LineTotal, "H", path + "\\results\\line_field_BS_Gmsh.txt");
	//*********************************************************************************************
	//BS correction
	//*********************************************************************************************
	vector<int> surf_ID = { 178 };
	vector<int> Elem_IDsSF;

	vector<vector<double>>pointsCoordinatesSF;
	vector<vector<int>>pointsIDSF;
	points.getGaussPointsVol(Elem_IDsSF, pointsCoordinatesSF, pointsIDSF, mesh, surf_ID);

	//-------------------------------------------------------------------------------------
	//External surface

	// Biot Savart inetgration
	vector<vector<double>> field1LineBS = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinatesSF, path);
	vector<vector<double>> field2LineBS = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinatesSF, path);
	vector<vector<double>> field1LineTotalBS = thisMath.sum(field1LineBS, field2LineBS);


	//get component of the field
	vector<vector<double>> tangent;
	vector<vector<double>> normal;
	vector<double> areas;

	post.getFieldComponents(normal, tangent, areas, field1LineTotalBS, pointsCoordinatesSF, Elem_IDsSF, pointsIDSF,  "Gauss_Points_H_field_surface",  path,  mesh);
	
}

void TestsBiotSavart::AdaptiveProcessWinding(double fw)
{

	////////////////////////////////////////////////////////////////
	//Files
	//string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\40_3D_Winding\\normal\\");
	string   path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\40_3D_Winding\\algo\\");
	string meshFile_points("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh_points(path + meshFile_points);
	Messages messages;

	////////////////////////////////////////////////////////////////
	//Inegration Points
	vector<vector<double>>gaussPointsCoord;
	vector<vector<int>>pointsIDPerElement;
	Operations points;
	vector<int>volID = {10000,10001,10002};
	vector<int>forceFourPoints = {10001};
	vector<int> elem_ID_list;

	int point_counter = 0;
	vector<int> gaussPoints = { 0 };
	//points.getGaussPointsVol(elem_ID_list, gaussPointsCoord, pointsIDPerElement, mesh_points, volID);
	points.getGaussPointsAdaptive(elem_ID_list, point_counter, gaussPointsCoord, pointsIDPerElement, mesh_points, gaussPoints, forceFourPoints);
	
	////////////////////////////////////////////////////////////////
	// Biot Savart inetgration
	BiotSavart biotSavart;
	double current_Density = 3.0*pow(10, 7);
	vector<double>centerPosition = { 0.0,0.0,0.0 };
	int volIDCurrente = 10000;
	int alongAxis = 2;
	vector<vector<double>> field = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_points, gaussPointsCoord, path);

	
	gaussPoints = { 1,2,3 };
	
	//////////////////////////////////////////////////////////////
	//Volume of all the elements
	Matrix jac(3, 3);
	int mumberOfElements = elem_ID_list.size();
	map<int, double> volume;
	Operations oper;
	vector<double> pFielduv = { 0.0,0.0,0.0 };
	GaussLegendrePoints thisElemGauss(4);
	double weight = thisElemGauss.weights[0];

	int counter = 0;
	for each (int elem in elem_ID_list)
	{
		int thisElemType = mesh_points.elemTypes[elem];
		jac = oper.Jacobian(thisElemType, elem, mesh_points, pFielduv);
		volume[elem] = abs(jac.Det_3x3())*weight;
		counter++;
	}


	//////////////////////////////////////////////////////////////////
	// Energy for the elements
	vector<double> energy;
	vector<pair<double, int>> energyElements;
	Vector1D thisMath;
	int numPoints = thisElemGauss.weights.size();
	double totalEnergy = 0;

	counter = 0;
	for each (int i in elem_ID_list)
	{
		 //Average mahgnetic field inside the element
		double absH = 0;
		for each (int thisPointID in pointsIDPerElement[counter])
		{
			absH = absH + thisMath.Abs(field[thisPointID]);
		}

		absH = numPoints*absH / ((double)pointsIDPerElement[counter].size());

		//Energy list
		double thisEnergy = pow(absH, 2.0)*volume[i] * 4.0* 3.14*pow(10, -7);
		totalEnergy += thisEnergy;
		energyElements.push_back(make_pair(thisEnergy, i));
		counter = counter + 1;
	}


	//////////////////////////////////////////////////////////////////
	//Select the elements with higher energy
	sort(energyElements.begin(), energyElements.end());
	vector<int> selectedElements;
	vector <int> elementsByEnergy;

	double selectedEnergy = 0;
	for (int i = elem_ID_list.size() - 1; i >= 0; i--)
	{
		if (selectedEnergy < totalEnergy*fw)
		{
			selectedElements.push_back(energyElements[i].second);
			selectedEnergy += energyElements[i].first;
		}
		else
			break;

	}
	///////////////////////////////////////////////////////////////////
	//// Adaptive process
	map<string, string> x;
	gaussPoints = { 1,2,3 };
	//int iteration_counter = 0;
	//int numberAddElements = 0.02*elem_ID_list.size();
	////error > error_energy
	double Wc = 0;
	////while (selectedEnergy > Wc)
	////{
	vector<int> newSelectedElements;
	////if (iteration_counter == 0)
	////{
		newSelectedElements = selectedElements;
	////}
	////else
	////{
	//	//Selects more elements
	//	int start = selectedElements.size();
	//	for (int i = elem_ID_list.size() - 1 - start; i >= elem_ID_list.size() - 1.0 - start - numberAddElements; i--)
	//	{
	//		newSelectedElements.push_back(energyElements[i].second);
	//		selectedElements.push_back(energyElements[i].second);
	//	}
	////}

	//Gets those new Gauss Points
	vector<vector<double>>gaussPointsCoordTwo;
	points.getGaussPointsAdaptive(selectedElements, point_counter, gaussPointsCoordTwo, pointsIDPerElement, mesh_points, gaussPoints, forceFourPoints);

	vector<vector<double>> field2 = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_points, gaussPointsCoordTwo, path);
	field.insert(field.end(), field2.begin(), field2.end());
	gaussPointsCoord.insert(gaussPointsCoord.end(), gaussPointsCoordTwo.begin(), gaussPointsCoordTwo.end());
	
	//int start3D = elem_ID_list[0];
	////Recalculates the energy
	//for each (int i in newSelectedElements)
	//{
	//	int elem_ID = i - start3D;
	//	// Average mahgnetic field inside the element
	//	double absH = 0;
	//	for each (int thisPointID in pointsIDPerElement[elem_ID])
	//	{
	//		absH = absH + thisMath.Abs(field[thisPointID]);
	//	}

	//	absH = numPoints*absH / ((double)pointsIDPerElement[elem_ID].size());

	//	//Energy list
	//	double thisEnergy = pow(absH, 2.0)*volume[i] * 4.0* 3.14*pow(10, -7);
	//	Wc += thisEnergy;
	//}

	counter = 0;
	for each (vector<int> elem in pointsIDPerElement)
	{
		if (elem.size() == 1)
		{
			vector<int> newElem = { elem[0],elem[0],elem[0],elem[0] };
			pointsIDPerElement[counter] = newElem;
		}
		counter++;
	}

	//////////////////////////////////////////////////////////////////
	// Post processing
	PostProcessing post;

	//writes resulting mesh
	mesh_points.writeMesh(mesh_points, path + "mesh_selected_" + to_string(fw) + ".txt", selectedElements, volID);

	//writes the field solution to Gmsh
	post.writeVectorField(gaussPointsCoord, field, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");

	//Writes the Gauss points information
	post.writeGaussPointsIDs(elem_ID_list, pointsIDPerElement, gaussPointsCoord, path);

	//Writes the magnetic field
	post.writeDataResults(field, path, "Gauss_Points_H_field");


}



void TestsBiotSavart::AdaptiveProcess(double fw)
{

	////////////////////////////////////////////////////////////////
	//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\37_Atuador_hor\\02_BS_algo\\");
	string meshFile_points("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh_points(path + meshFile_points);

	////////////////////////////////////////////////////////////////
	//Inegration Points
	vector<vector<double>>gaussPointsCoord;
	vector<vector<int>>pointsIDPerElement;
	Operations points;
	vector<int>volID = { 1000,1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008,1009 };
	vector<int> elem_ID_list;
	//std::iota(begin(elem_ID_list), end(elem_ID_list), 0);


	int point_counter = 0;
	vector<int> gaussPoints = { 0 };
	//points.getGaussPointsVol(elem_ID_list, gaussPointsCoord, pointsIDPerElement, mesh_points, volID);
	//points.getGaussPointsAdaptive(elem_ID_list, point_counter, gaussPointsCoord, pointsIDPerElement, mesh_points, gaussPoints);
	//int start3D = elem_ID_list[0];

	//////////////////////////////////////////////////////////////////
	//// Biot Savart inetgration
	//BiotSavart biotSavart;
	//double current_Density = 3.0*pow(10, 7);
	//vector<double>centerPosition = { 0.01,0.00,0.01 };
	//int volIDCurrente = 1006;
	//int alongAxis = 0;
	//vector<vector<double>> field = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_points, gaussPointsCoord, path);

	//vector<vector<double>> field2;
	//int counter_teste = 1;
	//for each (vector<int> points in pointsIDPerElement)
	//{
	//	int ID = (counter_teste * 4) - 1;
	//	field2.push_back(field[ID]);
	//	field2.push_back(field[ID]);
	//	field2.push_back(field[ID]);
	//	field2.push_back(field[ID]);
	//	counter_teste++;
	//}

	////////////////////////////////////////////////////////////////
	// Volume of all the elements
	//Matrix jac(3, 3);
	//int mumberOfElements = elem_ID_list.size();
	//map<int, double> volume;
	//Operations oper;
	//vector<double> pFielduv = { 0,0,0 };
	//GaussLegendrePoints thisElemGauss(4);
	//double weight = thisElemGauss.weights[0];

	//int counter = 0;
	//for each (int elem in elem_ID_list)
	//{
	//	int thisElemType = mesh_points.elemTypes[elem];
	//	jac = oper.Jacobian(thisElemType, elem, mesh_points, pFielduv);
	//	volume[elem] = (abs(jac.Det_3x3())*weight);
	//	counter++;
	//}

	//////////////////////////////////////////////////////////////////
	// Energy for the elements
	//vector<double> energy;
	//vector<pair<double, int>> energyElements;
	//Vector1D thisMath;
	//int numPoints = thisElemGauss.weights.size();
	//double totalEnergy = 0;

	//counter = 0;
	//for each (int i in elem_ID_list)
	//{
	//	 Average mahgnetic field inside the element
	//	double absH = 0;
	//	for each (int thisPointID in pointsIDPerElement[counter])
	//	{
	//		absH = absH + thisMath.Abs(field[thisPointID]);
	//	}

	//	absH = numPoints*absH / ((double)pointsIDPerElement[counter].size());

	//	Energy list
	//	double thisEnergy = pow(absH, 2.0)*volume[i] * 4.0* 3.14*pow(10, -7);
	//	totalEnergy += thisEnergy;

	//	energyElements.push_back(make_pair(thisEnergy, i));

	//	counter = counter + 1;
	//}

	//////////////////////////////////////////////////////////////////
	//Select the elements with higher energy
	Messages messages;
	//sort(energyElements.begin(), energyElements.end());
	//vector<int> selectedElements;
	//vector <int> elementsByEnergy;

	//double selectedEnergy = 0;
	//for (int i = elem_ID_list.size() - 1; i >= 0; i--)
	//{
	//	if (selectedEnergy < totalEnergy*fw)
	//	{
	//		selectedElements.push_back(energyElements[i].second);
	//		selectedEnergy += energyElements[i].first;
	//	}
	//	else
	//		break;

	//}

	///////////////////////////////////////////////////////////////////
	//// Adaptive process
	////map<string, string> x;
	gaussPoints = { 1,2,3 };
	//int iteration_counter = 0;
	//int numberAddElements = 0.02*elem_ID_list.size();
	////error > error_energy
	//double Wc = 0;
	////while (selectedEnergy > Wc)
	////{
	//vector<int> newSelectedElements;
	////if (iteration_counter == 0)
	////{
	//	newSelectedElements = selectedElements;
	////}
	////else
	////{
	//	//Selects more elements
	//	int start = selectedElements.size();
	//	for (int i = elem_ID_list.size() - 1 - start; i >= elem_ID_list.size() - 1.0 - start - numberAddElements; i--)
	//	{
	//		newSelectedElements.push_back(energyElements[i].second);
	//		selectedElements.push_back(energyElements[i].second);
	//	}
	////}



	//Gets those new Gauss Points
	//vector<vector<double>>gaussPointsCoordTwo;
	//points.getGaussPointsAdaptive(selectedElements, point_counter, gaussPointsCoordTwo, pointsIDPerElement, mesh_points, gaussPoints);
	//vector<vector<double>> field2 = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_points, gaussPointsCoordTwo, path);
	//field.insert(field.end(), field2.begin(), field2.end());
	//gaussPointsCoord.insert(gaussPointsCoord.end(), gaussPointsCoordTwo.begin(), gaussPointsCoordTwo.end());
	//pointsIDPerElement.insert(pointsIDPerElement.end(), pointsIDPerElementTwo.begin(), pointsIDPerElementTwo.end());

	////Recalculates the energy
	//for each (int i in newSelectedElements)
	//{
	//	int elem_ID = i - start3D;
	//	// Average mahgnetic field inside the element
	//	double absH = 0;
	//	for each (int thisPointID in pointsIDPerElement[elem_ID])
	//	{
	//		absH = absH + thisMath.Abs(field[thisPointID]);
	//	}

	//	absH = numPoints*absH / ((double)pointsIDPerElement[elem_ID].size());

	//	//Energy list
	//	double thisEnergy = pow(absH, 2.0)*volume[i] * 4.0* 3.14*pow(10, -7);
	//	Wc += thisEnergy;
	//}

	messages.logMessage("----------------------------------");
	//messages.logMessage("Iteraction: " + to_string(iteration_counter));
	//messages.logMessage("Pre calculated energy :" + to_string(selectedEnergy));
	//messages.logMessage("Real selected energy :" + to_string(Wc));
	messages.logMessage("# elements in the mesh:" + to_string(elem_ID_list.size()));

	//	iteration_counter++;
	//}

	//int counter = 0;
	//for each (vector<int> elem in pointsIDPerElement)
	//{
	//	if (elem.size() == 1)
	//	{
	//		vector<int> newElem = { elem[0],elem[0],elem[0],elem[0] };
	//		pointsIDPerElement[counter] = newElem;
	//	}
	//	counter++;
	//}


	////////////////////////////////////////////////////////////////////
	//// Post processing
	//PostProcessing post;

	////writes resulting mesh
	////mesh_points.writeMesh(mesh_points, path + "mesh_selected_" + to_string(fw) + ".txt", selectedElements, volID);

	////writes the field solution to Gmsh
	//post.writeVectorField(gaussPointsCoord, field, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");

	////Writes the Gauss points information
	//post.writeGaussPointsIDs(elem_ID_list, pointsIDPerElement, gaussPointsCoord, path);

	////Writes the magnetic field
	//post.writeDataResults(field, path, "Gauss_Points_H_field");


}


void TestsBiotSavart::AtuadorHor()
{

	////////////////////////////////////////////////////////////////
	//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\37_Atuador_hor\\04_FFEM_coarse\\");
	string meshFile_points("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh_points(path + meshFile_points);

	////////////////////////////////////////////////////////////////
	//Inegration Points
	vector<vector<int>>pointsID;
	Operations points;
	vector<int>volID = { 1000,1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008,1009 };
	vector<int> elem_ID_list;

	vector<int> Elem_IDs;
	vector<vector<double>>pointsCoordinates;
	points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh_points, volID);


	////////////////////////////////////////////////////////////////
	// Biot Savart inetgration
	BiotSavart biotSavart;
	double current_Density = 3.0*pow(10, 7);
	vector<double>centerPosition = { 0.01,0.00,0.01 };
	int volIDCurrente = 1006;
	int alongAxis = 0;
	vector<vector<double>> field = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh_points, pointsCoordinates, path);


	//////////////////////////////////////////////////////////////////
	// Post processing
	PostProcessing post;
	

	//writes the field solution to Gmsh
	post.writeVectorField(pointsCoordinates, field, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");

	//Writes the Gauss points information
	post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);

	//Writes the magnetic field
	post.writeDataResults(field, path, "Gauss_Points_H_field");
}





void TestsBiotSavart::AtuadorVert()
{

	////////////////////////////////////////////////////////////////
	//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\34_Atuador_vert\\teste\\");
	string meshFile_windings("windings.msh");
	string meshFile_model("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh_windings(path + meshFile_windings);
	GetMesh mesh_model(path + meshFile_model);
	Vector1D thisMath;


	////////////////////////////////////////////////////////////////
	//Inegration Points
	vector<vector<int>>pointsID;
	Operations points;
	vector<int>volID = { 5000,5001, 5002,5003,5004,5005,5006,5007,5008,5009 };
	//vector<int>volID = {5004,5005,5006,5007};


	vector<int> elem_ID_list;

	vector<int> Elem_IDs;
	vector<vector<double>>pointsCoordinates;
	points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh_model, volID);


	////////////////////////////////////////////////////////////////
	// Biot Savart inetgration
	BiotSavart biotSavart;
	vector<double>centerPosition_1 = { -0.04,0.00,0.00 };
	vector<double>centerPosition_2 = { 0.04,0.00,0.00 };
	int alongAxis = 2;
	double current_Density =2.727272*pow(10, 7);

	vector<vector<double>> field_1 = biotSavart.integrateSolidWinding(5002, -current_Density, centerPosition_1, alongAxis, mesh_windings, pointsCoordinates, path);
	vector<vector<double>> field_2 = biotSavart.integrateSolidWinding(5003, current_Density, centerPosition_2, alongAxis, mesh_windings, pointsCoordinates, path);
	vector<vector<double>> field_3 = biotSavart.integrateSolidWinding(5000, -current_Density, centerPosition_1, alongAxis, mesh_windings, pointsCoordinates, path);
	vector<vector<double>> field_4 = biotSavart.integrateSolidWinding(5001, current_Density, centerPosition_2, alongAxis, mesh_windings, pointsCoordinates, path);

	vector<vector<double>> fieldTotal = thisMath.sum(field_1, field_2);
	fieldTotal = thisMath.sum(fieldTotal, field_3);
	fieldTotal = thisMath.sum(fieldTotal, field_4);


	//////////////////////////////////////////////////////////////////
	// Post processing
	PostProcessing post;

	//writes the field solution to Gmsh
	post.writeVectorField(pointsCoordinates, fieldTotal, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");

	//Writes the Gauss points information
	post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);

	//Writes the magnetic field
	post.writeDataResults(fieldTotal, path, "Gauss_Points_H_field");


}


void TestsBiotSavart::AtuadorVert_VS()
{

	////////////////////////////////////////////////////////////////
	//Files
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\34_Atuador_vert\\FEMM_3D_VS\\");
	string meshFile_windings("windings.msh");
	string meshFile_model("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh_windings(path + meshFile_windings);
	GetMesh mesh_model(path + meshFile_model);
	Vector1D thisMath;


	//*********************************************************************************************
	//VS
	//*********************************************************************************************
	vector<vector<int>>pointsID;
	Operations points;
	vector<int>volID = {5004,5005,5006,5007};
	//vector<int>volID = { 5005,5006,5007 };


	//Points
	vector<int> elem_ID_list;
	vector<int> Elem_IDs;
	vector<vector<double>>pointsCoordinates;
	points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh_model, volID);

	// Biot Savart integration
	BiotSavart biotSavart;
	vector<double>centerPosition_1 = { -0.04,0.00,0.00 };
	vector<double>centerPosition_2 = { 0.04,0.00,0.00 };
	int alongAxis = 2;
	double current_Density = 2.727272*pow(10, 7);

	vector<vector<double>> field_1 = biotSavart.integrateSolidWinding(5002, -current_Density, centerPosition_1, alongAxis, mesh_windings, pointsCoordinates, path);
	vector<vector<double>> field_2 = biotSavart.integrateSolidWinding(5003, current_Density, centerPosition_2, alongAxis, mesh_windings, pointsCoordinates, path);
	vector<vector<double>> field_3 = biotSavart.integrateSolidWinding(5000, -current_Density, centerPosition_1, alongAxis, mesh_windings, pointsCoordinates, path);
	vector<vector<double>> field_4 = biotSavart.integrateSolidWinding(5001, current_Density, centerPosition_2, alongAxis, mesh_windings, pointsCoordinates, path);

	vector<vector<double>> fieldTotal = thisMath.sum(field_1, field_2);
	fieldTotal = thisMath.sum(fieldTotal, field_3);
	fieldTotal = thisMath.sum(fieldTotal, field_4);

	PostProcessing post;
	//writes the field solution to Gmsh
	post.writeVectorField(pointsCoordinates, fieldTotal, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");

	//Writes the Gauss points information
	post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);

	//Writes the magnetic field
	post.writeDataResults(fieldTotal, path, "Gauss_Points_H_field");


	//*********************************************************************************************
	//BS correction
	//*********************************************************************************************
	//vector<int> surf_ID = { 5015 };
	vector<int> surf_ID = { 5010 };
	vector<int> Elem_IDsSF;

	vector<vector<double>>pointsCoordinatesSurface;
	vector<vector<int>>pointsIDSF;
	points.getGaussPointsVol(Elem_IDsSF, pointsCoordinatesSurface, pointsIDSF, mesh_model, surf_ID);

	// Biot Savart inetgration
	vector<vector<double>> field_1_Surface = biotSavart.integrateSolidWinding(5002, -current_Density, centerPosition_1, alongAxis, mesh_windings, pointsCoordinatesSurface, path);
	vector<vector<double>> field_2_Surface = biotSavart.integrateSolidWinding(5003, current_Density, centerPosition_2, alongAxis, mesh_windings, pointsCoordinatesSurface, path);
	vector<vector<double>> field_3_Surface = biotSavart.integrateSolidWinding(5000, -current_Density, centerPosition_1, alongAxis, mesh_windings, pointsCoordinatesSurface, path);
	vector<vector<double>> field_4_Surface = biotSavart.integrateSolidWinding(5001, current_Density, centerPosition_2, alongAxis, mesh_windings, pointsCoordinatesSurface, path);

	vector<vector<double>> fieldTotal_Surface = thisMath.sum(field_1_Surface, field_2_Surface);
	fieldTotal_Surface = thisMath.sum(fieldTotal_Surface, field_3_Surface);
	fieldTotal_Surface = thisMath.sum(fieldTotal_Surface, field_4_Surface);

	//get component of the field
	vector<vector<double>> tangent;
	vector<vector<double>> normal;
	vector<double>areasSF;

	post.getFieldComponents(normal, tangent, areasSF, fieldTotal_Surface, pointsCoordinatesSurface, Elem_IDsSF, pointsIDSF, "Gauss_Points_H_field_surface", path, mesh_model);
	
	//*********************************************************************************************
	//BS correction
	//*********************************************************************************************
	vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	// Biot Savart inetgration
	vector<vector<double>> field_1_Line = biotSavart.integrateSolidWinding(5002, -current_Density, centerPosition_1, alongAxis, mesh_windings, points_list, path);
	vector<vector<double>> field_2_Line = biotSavart.integrateSolidWinding(5003, current_Density, centerPosition_2, alongAxis, mesh_windings, points_list, path);
	vector<vector<double>> field_3_Line = biotSavart.integrateSolidWinding(5000, -current_Density, centerPosition_1, alongAxis, mesh_windings, points_list, path);
	vector<vector<double>> field_4_Line = biotSavart.integrateSolidWinding(5001, current_Density, centerPosition_2, alongAxis, mesh_windings, points_list, path);

	vector<vector<double>> fieldTotal_Line = thisMath.sum(field_1_Line, field_2_Line);
	fieldTotal_Line = thisMath.sum(fieldTotal_Line, field_3_Line);
	fieldTotal_Line = thisMath.sum(fieldTotal_Line, field_4_Line);


	//Writes the magnetic field
	post.writeDataResults(fieldTotal_Line, path, "line_field_BS");

	//writes the field solution to Gmsh
	post.writeVectorField(points_list, fieldTotal_Line, "H", path + "\\results\\line_field_BS_Gmsh.txt");

}

void TestsBiotSavart::PerfectElectric()
{
	//*********************************************************************************************
	// VS corrected
	//*********************************************************************************************
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\35_Subdomain_Dular_2009_electric\\FFEM_Complete\\");
	string meshFile("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh(path + meshFile);
	PostProcessing post;
	BiotSavart biotSavart;
	Vector1D thisMath;

	//vector<int> vol = { 174 };
	//vector<int> Elem_IDs;

	//vector<vector<double>>pointsCoordinates;
	//vector<vector<int>>pointsID;

	Operations points;
	//points.getGaussPointsVol(Elem_IDs, pointsCoordinates, pointsID, mesh, vol);

	//-------------------------------------------------------------------------------------
	// Biot Savart integration
	double current_Density = pow(10, 7)*1000.0;
	vector<double>centerPosition = { 0.01,0.00,0.01 };
	int volIDCurrente = 173;
	int alongAxis = 0;
	//vector<vector<double>> field1 = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> field2 = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	//vector<vector<double>> fieldTotal = thisMath.sum(field1, field2);
	////-------------------------------------------------------------------------------------
	////Post processing
	//post.writeVectorField(pointsCoordinates, fieldTotal, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	//post.writeGaussPointsIDs(Elem_IDs, pointsID, pointsCoordinates, path);
	//post.writeDataResults(fieldTotal, path, "Gauss_Points_H_field");

	//-------------------------------------------------------------------------------------
	// Line
	vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	// Biot Savart inetgration
	vector<vector<double>> field1Line = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, points_list, path);
	vector<vector<double>> field2Line = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, points_list, path);
	vector<vector<double>> field1LineTotal = thisMath.sum(field1Line, field2Line);


	//Writes the magnetic field
	post.writeDataResults(field1LineTotal, path, "line_field_BS");

	//writes the field solution to Gmsh
	post.writeVectorField(points_list, field1LineTotal, "H", path + "\\results\\line_field_BS_Gmsh.txt");

	//*********************************************************************************************
	//BS correction and Bn
	//*********************************************************************************************
	vector<int> surf_ID = {177,178 };
	vector<int> Elem_IDsSF;

	vector<vector<double>>pointsCoordinatesSF;
	vector<vector<int>>pointsIDSF;
	points.getGaussPointsVol(Elem_IDsSF, pointsCoordinatesSF, pointsIDSF, mesh, surf_ID);

	//-------------------------------------------------------------------------------------
	//External surface

	// Biot Savart inetgration
	vector<vector<double>> field1LineBS = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinatesSF, path);
	vector<vector<double>> field2LineBS = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinatesSF, path);
	vector<vector<double>> field1LineTotalBS = thisMath.sum(field1LineBS, field2LineBS);


	//get component of the field
	vector<vector<double>> tangent;
	vector<vector<double>> normal;
	vector<double> areas;

	post.getFieldComponents(normal, tangent, areas, field1LineTotalBS, pointsCoordinatesSF, Elem_IDsSF, pointsIDSF, "Gauss_Points_H_field_surface", path, mesh);

}


void TestsBiotSavart::PerfectMagnetic()
{
	//*********************************************************************************************
	// Convert the mesh
	//*********************************************************************************************
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\36_Subdomain_Dular_2009_magnetic\\FFEM_Complete_finite\\");
	string meshFileTemp("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh meshTemp(path + meshFileTemp);
	PostProcessing post;
	BiotSavart biotSavart;
	Vector1D thisMath;

	// Layer of elements - Inner
	vector<int> vol = {175};
	int sufaceIDTemp = 177;
	vector<vector<double>>normaVectorsTemp;
	vector<double>areasTemp;
	vector<int> elemIDLayerTemp;
	

	//*********************************************************************************************
	// Fiels calculation
	//*********************************************************************************************
	string meshFile("model.msh");
	GetMesh mesh(path + meshFile);
	vector<int> elemIDLayer;
	vector<vector<int>>pointsID;
	vector<vector<double>>pointsCoordinates;
	Operations points;
	elemIDLayer = mesh.defineVolumeBoundary(mesh, sufaceIDTemp, vol, 2, normaVectorsTemp, areasTemp, path);
	points.getGaussPointsVol(elemIDLayer, pointsCoordinates, pointsID, mesh, vol);
	meshTemp.writeMesh(mesh, path + "\\layer_mesh.msh", elemIDLayer,vol);

	//-------------------------------------------------------------------------------------
	// Biot Savart integration
	double current_Density = pow(10, 7)*1000.0;
	vector<double>centerPosition = { 0.01,0.00,0.01 };
	int volIDCurrente = 173;
	int alongAxis = 0;
	vector<vector<double>> field1 = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	vector<vector<double>> field2 = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	vector<vector<double>> fieldTotal = thisMath.sum(field1, field2);

	//-------------------------------------------------------------------------------------
	//Field components along the element surface
	int counter = 0;
	vector<vector<double>>normalField;
	vector<vector<double>>tangentlField;
	vector<vector<double>> fieldTotalInv;
	
	for each (vector<double> n in normaVectorsTemp)
	{
		double vdotn = thisMath.dot(fieldTotal[counter], n);
		vector<double> this_normalField = thisMath.multiScal(n,vdotn);
		vector<double>this_tangentField = thisMath.subtract(fieldTotal[counter], this_normalField);
		normalField.push_back(this_normalField);
		tangentlField.push_back(thisMath.multiScal(this_tangentField, -1));
		fieldTotalInv.push_back(thisMath.multiScal(fieldTotal[counter], -1));
		counter++;
	}
	////-------------------------------------------------------------------------------------
	////Post processing
	post.writeVectorField(pointsCoordinates, tangentlField, "H", path + "\\results\\Gauss_Points_H_field_Gmsh.txt");
	post.writeGaussPointsIDs(elemIDLayer, pointsID, pointsCoordinates, path);
	post.writeDataResults(tangentlField, path, "Gauss_Points_H_field");
	/*post.writeVectorField(pointsCoordinates, normalField, "H", path + "\\results\\Gauss_Points_H_field_normal_Gmsh.txt");
	post.writeVectorField(pointsCoordinates, tangentlField, "H", path + "\\results\\Gauss_Points_H_field_tangential_Gmsh.txt");*/

	//-------------------------------------------------------------------------------------
	// Line
	//vector<vector<double>> points_list = post.readDataFile(path, "line_coordinates");

	//// Biot Savart inetgration
	//vector<vector<double>> field1Line = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, points_list, path);
	//vector<vector<double>> field2Line = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, points_list, path);
	//vector<vector<double>> field1LineTotal = thisMath.sum(field1Line, field2Line);


	////Writes the magnetic field
	//post.writeDataResults(field1LineTotal, path, "line_field_BS");

	////writes the field solution to Gmsh
	//post.writeVectorField(points_list, field1LineTotal, "H", path + "\\results\\line_field_BS_Gmsh.txt");

	//*********************************************************************************************
	//BS correction
	//*********************************************************************************************
	vector<int> surf_ID = { 178,177 };
	vector<int> Elem_IDsSF;

	vector<vector<double>>pointsCoordinatesSF;
	vector<vector<int>>pointsIDSF;
	points.getGaussPointsVol(Elem_IDsSF, pointsCoordinatesSF, pointsIDSF, mesh, surf_ID);

	//-------------------------------------------------------------------------------------
	//External surface

	// Biot Savart inetgration
	vector<vector<double>> field1LineBS = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinatesSF, path);
	vector<vector<double>> field2LineBS = biotSavart.integrateSolidWinding(176, -current_Density, centerPosition, alongAxis, mesh, pointsCoordinatesSF, path);
	vector<vector<double>> field1LineTotalBS = thisMath.sum(field1LineBS, field2LineBS);


	//get component of the field
	vector<vector<double>> tangent;
	vector<vector<double>> normal;
	vector<double>areasSF;

	post.getFieldComponents(normal, tangent, areasSF, field1LineTotalBS, pointsCoordinatesSF, Elem_IDsSF, pointsIDSF, "Gauss_Points_H_field_surface", path, mesh);

}