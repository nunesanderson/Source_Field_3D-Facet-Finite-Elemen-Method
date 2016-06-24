#include <string>
#include<vector>
using namespace std;
using std::string;

/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include "Matrix.h"

#ifndef _GET_MESH_INCLUDED_
#define _GET_MESH_INCLUDED_

class GetMesh
{
public:
	GetMesh(string filePath);
	void writeMesh(GetMesh mesh, string filePath, vector<int> IDs);
	vector<int> defineBoundary(GetMesh meshData, int phySurfaceFilter, int phyVolumeFilter, int atLeastNumNodes);
	~GetMesh();
	vector<int> physicalTags;
	vector<vector<double>> nodesCoordinates;
	vector<vector<int>> elemNodes;
	vector<int> elemTypes;
	vector<int> numNodesPerElement;
	vector<int> elementaryTags;
	int numNodes;
	int numElements;
};

class GetTxtData {
public:
	GetTxtData(string filePath);
	~GetTxtData();
	vector<string> lines;
	int numLines;
private:
};

class PostProcessing {
public:
	void writeVectorField(vector<vector<double>> coordinates, vector<vector<double>> fields, string fieldName, string path);
	void writeGaussPointsIDs(vector<vector<int>> pointsIDPerElement, vector<vector<double>> pointsCoordinates, string path);
	void writeDataResults(vector<vector<double>> twoDArrayData, string path,string fileName);

};
#endif