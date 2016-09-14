#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include<math.h>
using namespace std;
#include<vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include "Gmsh.h"
#include "Messages.h"
#include "ShapeFunctions.h"

GetMesh::~GetMesh() {}


GetMesh::GetMesh(string filePath)
{
	Messages messages;
	messages.logMessage("Reading Gmsh mesh");

	GetTxtData datafile(filePath);
	vector<string> data = datafile.lines;
	int rows = datafile.numLines;
	vector<vector <string>> allLinesList;

	size_t maxNumNodes = 0;
	bool checkMax = false;
	int startElements;
	int startNodes = 0;

	for each (string str in data)
	{

		istringstream iss(str);
		vector<string> thisLine;
		int counter = 0;
		do
		{
			string sub;
			iss >> sub;
			thisLine.push_back(sub);
		} while (iss);
		allLinesList.push_back(thisLine);

	}

	int counter = 0;
	size_t numberLines = allLinesList.size();
	int jumpLines = 0;
	for (int counter = 0; counter < numberLines; counter++)
	{
		vector <string> thisLine = allLinesList[counter];
		//Nodes coordinates
		if (thisLine[0] == "$Nodes")
		{
			string strNumNodes = allLinesList[counter + 1][0];
			numNodes = atoi(strNumNodes.c_str());
			startNodes = counter;

		}

		//Elements nodes
		if (thisLine[0] == "$Elements")
		{
			string strNustrnumElements = allLinesList[counter + 1][0];
			numElements = atoi(strNustrnumElements.c_str());
			checkMax = true;
			startElements = counter;
		}
		size_t thisLineLen = thisLine.size();
		if (thisLineLen > maxNumNodes && checkMax)
			maxNumNodes = thisLineLen;

	}

	//cout << "Nodes coordinates" << endl;
	nodesCoordinates.resize(numNodes);
	for (int i = 0; i < numNodes; i++)
	{
		nodesCoordinates[i].resize(3);
		nodesCoordinates[i][0] = stod(allLinesList[startNodes + 2 + i][1].c_str());
		nodesCoordinates[i][1] = stod(allLinesList[startNodes + 2 + i][2].c_str());
		nodesCoordinates[i][2] = stod(allLinesList[startNodes + 2 + i][3].c_str());
		//cout << nodesCoordinates[i][0] << " , " << nodesCoordinates[i][1] << " , " << nodesCoordinates[i][2] << endl;
		jumpLines++;
	}
	elemNodes.resize(numElements);
	for (int i = 0; i < numElements; i++)
	{
		//Nodes of the elements
		vector<string> thisElement = allLinesList[i + startElements + 2];
		size_t numStrs = allLinesList[i + startElements + 2].size();
		int colCounter = 0;
		for (int j = 5; j < numStrs - 1; j++)
		{
			string thisNode = thisElement[j];
			elemNodes[i].push_back(stoi(thisNode.c_str()) - 1);
			//cout << elemNodes[i][colCounter] << ",";
			colCounter++;
		}
		numNodesPerElement.push_back(colCounter);

		//Type of the elements
		elemTypes.push_back(stoi(thisElement[1].c_str()));

		//Tags of the elements
		physicalTags.push_back(stoi(thisElement[3].c_str()));
		elementaryTags.push_back(stoi(thisElement[4].c_str()));
	}
	messages.logMessage("Reading Gmsh mesh: Done");
}


void GetMesh::writeMesh(GetMesh mesh, string filePath, vector<int> IDs) {

	Messages messages;
	messages.logMessage("Writing Gmsh mesh");

	vector<vector<double>>  nodesCoordinates = mesh.nodesCoordinates;
	vector<vector<int>> elemNodes = mesh.elemNodes;
	vector<int> elemTypes = mesh.elemTypes;
	vector<int> physicalTasg = mesh.physicalTags;
	vector<int> elementaryTags = mesh.elementaryTags;
	vector<int> numNodesPerElem = mesh.numNodesPerElement;
	int numNodes = mesh.numNodes;
	int numElements = mesh.numElements;

	ofstream myfile;
	myfile.open(filePath);

	myfile << "$MeshFormat\n";
	myfile << "2.2 0 8\n";
	myfile << "$EndMeshFormat\n";
	myfile << "$Nodes\n";
	myfile << numNodes << "\n";
	for (int i = 0; i < numNodes; i++)
	{
		myfile << i + 1 << " " << nodesCoordinates[i][0] << " " << nodesCoordinates[i][1] << " " << nodesCoordinates[i][2] << "\n";
	}

	size_t sizeIDs = IDs.size();
	myfile << "$EndNodes\n";
	myfile << "$Elements\n";
	myfile << sizeIDs << "\n";

	for (int i = 0; i < sizeIDs; i++)
	{
		int thisID = IDs[i];
		myfile << i + 1 << " " << elemTypes[thisID] << " " << "2" << " " << physicalTags[thisID] << " " << elementaryTags[thisID] << " ";

		for (size_t k = 0; k < elemNodes[thisID].size(); k++)
		{
			myfile << elemNodes[thisID][k] + 1;
			if (k < elemNodes[thisID].size() - 1)
			{
				myfile << " ";
			}
		}
		myfile << "\n";
	}
	myfile << "$EndElements\n";

	myfile.close();
	messages.logMessage("Writing Gmsh mesh: Done");
}



vector<int> GetMesh::defineBoundary(GetMesh mesh, int phySurfaceFilter, vector <int> phyVolumeFilter, int atLeastNumNodes, vector<vector<double>> &normalVectors) {

	vector<vector<double>>  nodesCoordinates = mesh.nodesCoordinates;
	vector<vector<int>> elemNodes = mesh.elemNodes;
	vector<int> elemTypes = mesh.elemTypes;
	vector<int> physicalTasg = mesh.physicalTags;
	vector<int> elementaryTags = mesh.elementaryTags;
	vector<int> numNodesPerElem = mesh.numNodesPerElement;
	int numNodes = mesh.numNodes;
	int numElements = mesh.numElements;
	Operations oper;
	vector<int> elementsBothDomains;

	vector < vector<int>>  twoDElementNodes;
	vector<int> globalID;
	vector<vector<double>> printNormalVectorPosition;


	for (int i = 0; i < numElements; i++)
	{

		if (elemTypes[i] < 3 && physicalTags[i] == phySurfaceFilter)
		{

			//Get the nodes of this 2D element 
			vector<int> this2DElemNodes;
			for (size_t k = 0; k < elemNodes[i].size(); k++)
			{
				this2DElemNodes.push_back(elemNodes[i][k]);
			}

			for (int elemCounter = 0; elemCounter < numElements; elemCounter++)
			{

				if (elemTypes[elemCounter] >= 4 && std::find(phyVolumeFilter.begin(), phyVolumeFilter.end(), physicalTags[elemCounter]) != phyVolumeFilter.end())
				{
					vector<int> this3DElemNodes;

					//Get the nodes of this 3D element
					int numNodes3D = numNodesPerElem[elemCounter];
					for (int k = 0; k < numNodes3D; k++)
					{
						this3DElemNodes.push_back(elemNodes[elemCounter][k]);
					}

					//Verify if the 3D element contains all the nodes of the 2D element
					bool addThisElem = false;
					size_t numNodes2D = this2DElemNodes.size();
					int foundCounter = 0;
					for (int DElemCounter = 0; DElemCounter < numNodes2D; DElemCounter++)
					{
						if (std::find(this3DElemNodes.begin(), this3DElemNodes.end(), this2DElemNodes[DElemCounter]) != this3DElemNodes.end())
							foundCounter++;

						if (foundCounter == atLeastNumNodes)
						{
							addThisElem = true;
							break;
						}
					}

					if (addThisElem)
					{
						//Gets this element ID
						elementsBothDomains.push_back(elemCounter);
						twoDElementNodes.push_back(this2DElemNodes);
						globalID.push_back(i);
					}
				}
			}
		}
	}

	//Computes the normal vector
	int counter = 0;
	for each (int elemID in elementsBothDomains)
	{
		vector<int> this2DElemNodes = twoDElementNodes[counter];
		int PointID = this2DElemNodes[0];
		vector<double> P1 = { nodesCoordinates[PointID][0],nodesCoordinates[PointID][1],nodesCoordinates[PointID][2] };

		PointID = this2DElemNodes[1];
		vector<double> P2 = { nodesCoordinates[PointID][0],nodesCoordinates[PointID][1],nodesCoordinates[PointID][2] };

		PointID = this2DElemNodes[2];
		vector<double> P3 = { nodesCoordinates[PointID][0],nodesCoordinates[PointID][1],nodesCoordinates[PointID][2] };

		Vector1D thisMath;
		vector<double> A = thisMath.subtract(P2, P1);
		vector<double> B = thisMath.subtract(P3, P1);
		vector<double> AxB = thisMath.crossProduct(A, B);
		double absAxB = thisMath.Abs(AxB);
		vector<double> normal = thisMath.multiScal(AxB, 1.0 / absAxB);

		int thisElemType = elemTypes[elemID];
		GaussLegendrePoints thisElemGauss(thisElemType);
		for (int pointCounter = 0; pointCounter < thisElemGauss.pointsCoordinates.rows; pointCounter++)
		{
			//UVP
			vector<double>pFielduv;
			pFielduv = thisElemGauss.pointsCoordinates.mat[pointCounter];

			//XYZ
			vector<double> pFieldxy = oper.scalLocalToReal(thisElemType, elemID, mesh, pFielduv);

			normalVectors.push_back(normal);
			printNormalVectorPosition.push_back(pFieldxy);
		}

		counter++;
	}

	PostProcessing post;
	string path("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\31_Subdomain_Dular_2009\\Subdomain");
	post.writeVectorField(printNormalVectorPosition, normalVectors, "NormalVector", path + "\\results\\Gmsh_Normal_Vector.txt");

	return elementsBothDomains;
}

//vector<vector<double>> GetMesh::defineNormaVectors(GetMesh mesh, int phySurfaceFilter,vector<int> boundaryElementsList)
//{
//	vector<vector<double>>  nodesCoordinates = mesh.nodesCoordinates;
//	vector<vector<int>> elemNodes = mesh.elemNodes;
//	vector<int> elemTypes = mesh.elemTypes;
//	vector<int> physicalTasg = mesh.physicalTags;
//	vector<int> elementaryTags = mesh.elementaryTags;
//	vector<int> numNodesPerElem = mesh.numNodesPerElement;
//	int numNodes = mesh.numNodes;
//	int numElements = mesh.numElements;
//
//	vector<int> elementsBothDomains;
//
//	for (int i = 0; i < numElements; i++)
//	{
//
//		if (elemTypes[i] < 3 && physicalTags[i] == phySurfaceFilter)
//		{
//
//			Get the nodes of this 2D element 
//			vector<int> this2DElemNodes;
//			for (size_t k = 0; k < elemNodes[i].size(); k++)
//			{
//				this2DElemNodes.push_back(elemNodes[i][k]);
//
//			}
//
//	return vector<vector<double>>();
//}

GetTxtData::~GetTxtData(void) {}
GetTxtData::GetTxtData(string filePath) {
	string line;
	ifstream myfile(filePath);
	numLines = 0;
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			lines.push_back(line);
			numLines++;

		}
		myfile.close();
	}
	else cout << "Unable to open file";
}

void PostProcessing::writeVectorField(vector<vector<double>> coordinates, vector<vector<double>> fields, string fieldName, string path) {

	Messages messages;
	messages.logMessage("Writing Gmsh vector field " + fieldName);

	ofstream myfile;
	myfile.open(path);
	myfile << "View \"" << fieldName << "\" { \n";
	for (size_t i = 0; i < coordinates.size(); i++)
	{
		myfile << "VP(" << coordinates[i][0] << "," << coordinates[i][1] << "," << coordinates[i][2] << ")";
		myfile << "{" << fields[i][0] << "," << fields[i][1] << "," << fields[i][2] << "};\n";
	}

	myfile << "TIME{ 1 };\n};";
	myfile.close();
	messages.logMessage("Writing Gmsh vector field " + fieldName + ":Done");


}



void PostProcessing::writeGaussPointsIDs(vector<int>elemIDs, vector<vector<int>> pointsIDPerElement, vector<vector<double>> pointsCoordinates, string path)
{
	Messages messages;
	messages.logMessage("Writing Gauss points");

	string fileNameID("results//Gauss_Points_IDs.txt");
	string fileNameCoord("results//Gauss_Points_Coordinates.txt");
	ofstream myfile;
	string sep = " ";

	/* ------------------------------------------------------------------------
	Writes the IDS
	---------------------------------------------------------------------------*/
	myfile.open(path + fileNameID);
	int elemCounter = 0;
	size_t sizeElemIDs = elemIDs.size();
	for each (vector<int> thisElemPoints in pointsIDPerElement)
	{
		if (sizeElemIDs > 0)
			myfile << elemIDs[elemCounter] << sep;
		else
			myfile << elemCounter << sep;

		size_t counter = 0;
		for each (int thisID in thisElemPoints)
		{
			myfile << thisID;
			if (counter < thisElemPoints.size() - 1)
			{
				myfile << sep;
			}
			counter++;
		}
		myfile << "\n";
		elemCounter++;
	}
	myfile.close();

	/* ------------------------------------------------------------------------
	Writes the coordinates
	---------------------------------------------------------------------------*/
	myfile.open(path + fileNameCoord);
	for each (vector<int> thisElemPoints in pointsIDPerElement)
	{
		for each (int i in thisElemPoints)
		{
			myfile << pointsCoordinates[i][0] << " " << pointsCoordinates[i][1] << " " << pointsCoordinates[i][2] << "\n";
		}
	}
	myfile.close();

	messages.logMessage("Writing Gauss points: Done");
}

void PostProcessing::writeDataResults(vector<vector<double>> twoDArrayData, string path, string fileName)
{
	Messages messages;
	messages.logMessage("Writing data file " + fileName);

	string filePath(path + "results//" + fileName + ".txt");
	ofstream myfile;
	myfile.open(filePath);
	for each (vector<double> thisRow in twoDArrayData)
	{
		size_t counter = 0;
		for each (double thisData in thisRow)
		{
			myfile << thisData;
			if (counter < thisRow.size() - 1)
			{
				myfile << " ";
			}
			counter++;
		}
		myfile << "\n";
	}
	myfile.close();

	messages.logMessage("Writing data file " + fileName + ": Done");

}

vector<vector<double>> PostProcessing::readDataFile(string path, string fileName)
{

	Messages messages;
	messages.logMessage("Reading data file");

	string filePath(path + "results//" + fileName + ".txt");
	GetTxtData datafile(filePath);
	vector<string> data = datafile.lines;
	int rows = datafile.numLines;
	vector<vector <double>> allLinesList;


	for each (string str in data)
	{
		istringstream iss(str);
		vector<double> thisLine;
		int counter = 0;
		do
		{
			string sub;
			iss >> sub;
			if (sub != "")
			{
				thisLine.push_back(stod(sub));

			}
		} while (iss);
		allLinesList.push_back(thisLine);
	}


	return allLinesList;
}
