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

vector<double> BiotSavart::biotSavartEquationThreeD(vector <double> fieldPoint, vector <double> currentPoint, vector <double> dlUnitary, double current)
{
	Vector1D thisMath;
	vector<double> R = thisMath.subtract(fieldPoint, currentPoint);
	double R_abs = thisMath.Abs(R);
	vector<double> vector_product = thisMath.crossProduct(dlUnitary, R);
	double k;
	
	if (R_abs==0.0)
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

vector<double> BiotSavart::biotSavartEquationTwoD(vector <double> fieldPoint, vector <double> currentPoint, vector <double> dlUnitary, double current)
{
	Vector1D thisMath;
	vector<double> R = thisMath.subtract(fieldPoint, currentPoint);
	double R_abs = thisMath.Abs(R);
	vector<double> vector_product = thisMath.crossProduct(dlUnitary, R);
	double k;
	double lim = 1*pow(10, -4);
	if (R_abs < lim)
	{
		k = 0;
	}
	else
	{

		k = current / (2.0*PI*pow(R_abs, 2.0));
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
				double theta = atan2(pFieldxy[2]-centerPosition[2], pFieldxy[0]-centerPosition[0]);
				vector<double> currentDensityVec = { current_Density*sin(theta) ,0,-current_Density*cos(theta) };

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
					vector<double> dHPoint = biotSavartEquationThreeD(P_coord, pFieldxy, currentDensityVec, 1);
					Hresults[pointCounter][0] += dHPoint[0] * weight*detJac;
					Hresults[pointCounter][1] += dHPoint[1] * weight*detJac;
					Hresults[pointCounter][2] += dHPoint[2] * weight*detJac;
				}

			}

		}
	}
	messages.logMessage("BiotSavart - Solid domain integration: Done");

	PostProcessing post;
	post.writeVectorField(currentDensityListCoord, currentDensity, "Current Density", filePath + "\\results\\Gmsh_Current_Density.txt");

	return Hresults;
}

void BiotSavart::integrateTwoD(vector<vector<double>> &Hresults, int volID, double current_Density, GetMesh mesh, vector<vector<double>> gaussPointsData, string filePath)
{
	Messages messages;
	messages.logMessage("BiotSavart - 2D domain integration");
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
	vector<vector<double>> currentDensityListCoord;
	vector<vector<double>> currentDensity;

	Matrix jac(3, 3);

	//Types of elements
	vector<int> validElemTypes = { 2,3,9,10};

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
				vector<double> currentDensityVec = { 0 ,0,current_Density };

				//save the current density plot
				currentDensity.push_back(currentDensityVec);
				currentDensityListCoord.push_back(pFieldxy);

				// Get dH using the Biot-Savart equation
				jac = oper.Jacobian(thisElemType, i, mesh, pFielduv);
				jac.resize(2, 2);
				double detJac = abs(jac.Det_3x3());
				double weight = thisElemGauss.weights[pointCounter];

				for (int pointCounter = 0; pointCounter < numberGaussPoints; pointCounter++)
				{
					//Point to get H
					vector <double> P_coord = { gaussPointsData[pointCounter][0], gaussPointsData[pointCounter][1], gaussPointsData[pointCounter][2] };
					vector<double> dHPoint = biotSavartEquationTwoD(P_coord, pFieldxy, currentDensityVec, 1);
					Hresults[pointCounter][0] += dHPoint[0] * weight*detJac;
					Hresults[pointCounter][1] += dHPoint[1] * weight*detJac;

				}

			}

		}
	}
	messages.logMessage("BiotSavart - 2D domain integration: Done");

	PostProcessing post;
	post.writeVectorField(currentDensityListCoord, currentDensity, "Current Density", filePath + "\\results\\Gmsh_Current_Density.txt");

}

vector<vector<double>>  BiotSavart::integrateLine(double current,int volID,GetMesh mesh, vector<vector<double>> gaussPointsData, string path) {
	Messages messages;
	messages.logMessage("BiotSavart - Line integration");
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

	//Integration loop
	for (int i = 0; i < numElements; i++)
	{
		vector<vector<double>> path;

		int thisElemType = elemTypes[i];
		if (physicalTag[i]==volID && (thisElemType == 1 || thisElemType == 8))
		{

		// loop for the Gauss points
			GaussLegendrePoints thisElemGauss(thisElemType);
			for (int pointCounter = 0; pointCounter < thisElemGauss.pointsCoordinates.rows; pointCounter++)
			{
				//Integration point @UVP
				vector<double>pFielduv;
				pFielduv = thisElemGauss.pointsCoordinates.mat[pointCounter];
				double weight = thisElemGauss.weights[pointCounter];

				//Integration point @XYZ
				vector<double> pFieldxy = oper.scalLocalToReal(thisElemType, i, mesh, pFielduv);

				//Jacobian evaluated at Gauss point
				Matrix jac = oper.Jacobian(thisElemType, i, mesh, pFielduv);
				double detJac = oper.getDetJac1D(jac);
				//dl vector
				vector<double> dlVect = { jac.mat[0][0], jac.mat[0][1], jac.mat[0][2] };
				vector<double> dlUnitary = thisMath.multiScal(dlVect, 1. / thisMath.Abs(dlVect));

				//save the information to plot the dl
				if (plotDl)
				{
					dlListField.push_back(dlUnitary);
					dlListCoord.push_back(pFieldxy);
				}

				for (int pointCounter = 0; pointCounter < numberGaussPoints; pointCounter++)
				{
					//Point to get H
					vector <double> P_coord = { gaussPointsData[pointCounter][0], gaussPointsData[pointCounter][1], gaussPointsData[pointCounter][2] };
					vector<double> dHPoint = biotSavartEquationThreeD(P_coord, pFieldxy, dlUnitary, current);
					Hresults[pointCounter][0] += dHPoint[0] * weight*detJac;
					Hresults[pointCounter][1] += dHPoint[1] * weight*detJac;
					Hresults[pointCounter][2] += dHPoint[2] * weight*detJac;
				}
			}
		}
		else
		{
			break;
		}
	}

	PostProcessing teste;
	teste.writeVectorField(dlListCoord, dlListField, "dl's", path + "\\results\\dl.txt");

	messages.logMessage("BiotSavart - Line integration: Done");
	return Hresults;
}