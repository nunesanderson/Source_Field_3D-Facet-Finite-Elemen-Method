#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include<math.h>
using namespace std;
#include<vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>

/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include "Biot-Savart.h"
#include "Matrix.h"
#include "ShapeFunctions.h"
#include "Messages.h"


const double PI = 3.1415926535897;

vector<double> BiotSavart::biotSavartEquation(vector <double> fieldPoint, vector <double> currentPoint, vector <double> dlUnitary, double current)
{
	Vector1D thisMath;
	vector<double> R = thisMath.subtract(fieldPoint, currentPoint);
	double R_abs = thisMath.Abs(R);
	vector<double> vector_product = thisMath.crossProduct(dlUnitary, R);
	double k;
	if (R_abs == 0.0)
	{
		k = 0;
	}
	else
	{

		k = current / (4.0*PI*pow(R_abs, 3.0));
	}

	vector<double> dH_P = thisMath.multiScal(vector_product, k);

	return dH_P;
}
BiotSavart::~BiotSavart() {}

vector<vector<double>>  BiotSavart::integrateSolidWinding(int volID, double current_Density, vector<double> centerPosition, int alongAxis, GetMesh mesh, vector<vector<double>> gaussPointsData, string filePath) {
	
	Messages messages;
	messages.logMessage("BiotSavart - Solid domain integration");
	//Get the mesh information
	vector<vector<double>>  nodesCoordinates = mesh.nodesCoordinates;
	vector<vector<int>> elemNodes = mesh.elemNodes;
	vector<int> elemTypes = mesh.elemTypes;
	vector<int> physicalTag = mesh.physicalTags;
	vector<int> elementaryTags = mesh.elementaryTags;
	vector<int> numNodesPerElem = mesh.numNodesPerElement;
	int numNodes = mesh.numNodes;
	int numElements = mesh.numElements;

	//classes instantiation
	Operations oper;
	Vector1D thisMath;

	//Initialize data structures
	int numberGaussPoints = gaussPointsData.size();
	vector<vector<double>> Hresults(numberGaussPoints, vector<double>(3, 0));
	vector<vector<double>> currentDensityListCoord;
	vector<vector<double>> currentDensity;

	Matrix jac(3, 3);

	//Types of elements
	vector<int> validElemTypes = { 4,5,6,7,11 };

	//Integration loop
	for (int i = 0; i < numElements; i++)
	{
		int thisElemType = elemTypes[i];

		if ((std::find(validElemTypes.begin(), validElemTypes.end(), thisElemType) != validElemTypes.end()) && (physicalTag[i] == volID))
		{

			// loop for the Gauss points
			GaussLegendrePoints thisElemGauss(thisElemType);
			vector<double> dHElement = { 0, 0, 0 };
			for (int pointCounter = 0; pointCounter < thisElemGauss.pointsCoordinates.rows; pointCounter++)
			{

				//Integration point @UVP
				vector<double>pFielduv;
				pFielduv = thisElemGauss.pointsCoordinates.mat[pointCounter];

				//Integration point @XYZ
				vector<double> pFieldxy = oper.scalLocalToReal(thisElemType, i, mesh, pFielduv);

				//current density vector
				double theta = atan2(pFieldxy[1], pFieldxy[0]);
				vector<double> currentDensityVec = { current_Density*sin(theta) ,-current_Density*cos(theta),0 };

				//save the current density plot
				currentDensity.push_back(currentDensityVec);
				currentDensityListCoord.push_back(pFieldxy);

				// Get dH using the Biot-Savart equation
				jac = oper.Jacobian(thisElemType, i, mesh, pFielduv);
				double detJac = abs(jac.Det_3x3());
				double weight = thisElemGauss.weights[pointCounter];

				for (int pointCounter = 0; pointCounter < numberGaussPoints; pointCounter++)
				{
					//Point to get H
					vector <double> P_coord = { gaussPointsData[pointCounter][0], gaussPointsData[pointCounter][1], gaussPointsData[pointCounter][2] };
					vector<double> dHPoint = biotSavartEquation(P_coord, pFieldxy, currentDensityVec, 1);
					Hresults[pointCounter][0] += dHPoint[0]* weight*detJac;
					Hresults[pointCounter][1] += dHPoint[1]*weight*detJac;
					Hresults[pointCounter][2] += dHPoint[2]*weight*detJac;
				}

			}

		}
	}
	messages.logMessage("BiotSavart - Solid domain integration: Done");

	PostProcessing post;
	post.writeVectorField(currentDensityListCoord, currentDensity, "Current Density", filePath + "\\results\\Gmsh_Current_Density.txt");

	return Hresults;
}

vector<vector<double>>  BiotSavart::integrateLine(GetMesh mesh, vector<vector<double>> gaussPointsData, string path) {

	//Get the mesh information
	vector<vector<double>>  nodesCoordinates = mesh.nodesCoordinates;
	vector<vector<int>> elemNodes = mesh.elemNodes;
	vector<int> elemTypes = mesh.elemTypes;
	vector<int> physicalTag = mesh.physicalTags;
	vector<int> elementaryTags = mesh.elementaryTags;
	vector<int> numNodesPerElem = mesh.numNodesPerElement;
	int numNodes = mesh.numNodes;
	int numElements = mesh.numElements;

	//classes instantiation
	Operations oper;
	Vector1D thisMath;

	//Initialize data structures
	int numberGaussPoints = gaussPointsData.size();
	vector<vector<double>> Hresults(numberGaussPoints, vector<double>(3, 0));
	vector<vector<double>> dlListCoord;
	vector<vector<double>> dlListField;

	bool plotDl = true;
	double current = 1;

	// Loop to get H for all the points at gaussPointsData
	for (int pointCounter = 0; pointCounter < numberGaussPoints; pointCounter++)
	{
		//Point to get H
		vector <double> P_coord = { gaussPointsData[pointCounter][0], gaussPointsData[pointCounter][1], gaussPointsData[pointCounter][2] };

		//H result at point P_coord
		vector <double> Hp = { 0, 0, 0 };

		//Integratio loop
		for (int i = 0; i < numElements; i++)
		{
			vector<vector<double>> path;

			int thisElemType = elemTypes[i];
			if (thisElemType == 1 || thisElemType == 8)
			{

				// Coordinates for first and last element points
				int nodeA = elemNodes[i][0];
				int nodeB = elemNodes[i][1];

				vector <double> dl_start;
				dl_start = nodesCoordinates[nodeA];
				vector <double> dl_end;
				dl_start = nodesCoordinates[nodeB];
				path.push_back(dl_start);

				// loop for the Gauss points
				GaussLegendrePoints thisElemGauss(thisElemType);
				vector<double> dH = { 0, 0, 0 };
				for (int pointCounter = 0; pointCounter < thisElemGauss.pointsCoordinates.rows; pointCounter++)
				{
					//Integration point @UVP
					vector<double>pFielduv;
					pFielduv = thisElemGauss.pointsCoordinates.mat[pointCounter];

					//Integration point @XYZ
					vector<double> pFieldxy = oper.scalLocalToReal(thisElemType, i, mesh, pFielduv);

					//Jacobian evaluated at Gauss point
					Matrix jac = oper.Jacobian(thisElemType, i, mesh, pFielduv);

					//dl vector
					vector<double> dlVect = { jac.mat[0][0], jac.mat[0][1], jac.mat[0][2] };
					vector<double> dlUnitary = thisMath.multiScal(dlVect, 1. / thisMath.Abs(dlVect));

					//save the informqtion to plot the dl
					if (plotDl)
					{
						dlListField.push_back(dlUnitary);
						dlListCoord.push_back(pFieldxy);
					}

					path.push_back(pFieldxy);
					dH = thisMath.sum(dH, biotSavartEquation(P_coord, pFieldxy, dlUnitary, current));
				}

				//Get the lenght of the line element
				path.push_back(dl_end);
				double lenght = 0;
				for (int i = 0; i < path.size() - 1; i++)
				{
					lenght += thisMath.distance(path[i], path[i + 1]);
				}

				//Integration sum
				dH = thisMath.multiScal(dH, 1.0 / (double)thisElemGauss.pointsCoordinates.rows *lenght);
				Hp = thisMath.sum(Hp, dH);
			}
			else
			{
				break;
			}
		}
		Hresults[pointCounter] = Hp;
	}


	PostProcessing teste;
	teste.writeVectorField(dlListCoord, dlListField, "dl's", path + "\\results\\dl.txt");

	return Hresults;
}