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


TestsBiotSavart::TestsBiotSavart()
{

	Messages messages;
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\23_Biot_Savart_3D\\Solid winding\\");
	string meshFile("model.msh");
	string HFieldFile("results\\H_Field.txt");
	GetMesh mesh(path + meshFile);

//Points
	vector<vector<double>>pointsCoordinates;
	vector<vector<int>>pointsID;
	Operations points;
	points.getGaussPoints(pointsCoordinates, pointsID, mesh, 101);


// Analytical
	BiotSavartAnalyt analy;
	vector<double> analyPoints;
	vector<vector<double>> plotAnaly;
	for each (vector<double> thisPoint in pointsCoordinates)
	{
		analyPoints.push_back(thisPoint[2]);
	}
	
	vector<double> Hanaly = analy.getHcoil(50.0/1000.0, analyPoints,0.5,1.0,100.0);
	int counter = 0;
	for each (double thisData in Hanaly)
	{
		vector<double> line;
		line.push_back(analyPoints[counter]);
		line.push_back(thisData);
		plotAnaly.push_back(line);
		counter++;
	}

	
	// Biot Savart inetgration
	BiotSavart biotSavart;
	double current_Density = 20000.0;
	vector<double>centerPosition = { 0.0,0.0,0.0 };
	int volIDCurrente = 100;
	int alongAxis = 0;
	vector<vector<double>> field = biotSavart.integrateSolidWinding(volIDCurrente, current_Density, centerPosition, alongAxis, mesh, pointsCoordinates, path);
	vector<vector<double>> BiotSavart;

	counter = 0;
	for each (vector<double> thisData in field)
	{
		vector<double> line;
		line.push_back(pointsCoordinates[counter][2]);
		line.push_back(abs(thisData[2]));
		BiotSavart.push_back(line);
		counter++;
	}


	// Post processing
	PostProcessing post;
	post.writeVectorField(pointsCoordinates, field, "H", path + "\\results\\Gmsh_H_vector.txt");
	post.writeDataResults(plotAnaly, path, "analyt");
	post.writeDataResults(BiotSavart, path, "Biot");



}
