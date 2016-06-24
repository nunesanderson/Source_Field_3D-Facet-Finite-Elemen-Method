#include <string>
#include<vector>
using namespace std;
using std::string;

/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include"Gmsh.h"


#ifndef _BIOT_SAVART_INCLUDED_
#define _BIOT_SAVART_INCLUDED_

class BiotSavart
{
public:
	vector<double> biotSavartEquation(vector <double> fieldPoint, vector <double> currentPoint, vector <double> dlUnitary, double current);
	vector<vector<double>>  integrateSolidWinding(int volID, double current_Density, vector<double> centerPosition, int alongAxis, GetMesh meshdata, vector<vector<double>> gaussPointsData, string filePath);
	vector<vector<double>>  integrateLine(GetMesh meshdata, vector<vector<double>> gaussPointsData, string path);
	~BiotSavart();
};
#endif