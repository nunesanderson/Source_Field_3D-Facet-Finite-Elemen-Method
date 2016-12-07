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
	vector<double> biotSavartEquationThreeD(vector <double> fieldPoint, vector <double> currentPoint, vector <double> dlUnitary, double current);
	vector<double> biotSavartEquationTwoD(vector <double> fieldPoint, vector <double> currentPoint, vector <double> dlUnitary, double current);

	vector<vector<double>>  integrateSolidWinding(int volID, double current_Density, vector<double> centerPosition, int alongAxis, GetMesh meshdata, vector<vector<double>> gaussPointsData, string filePath);
	void integrateTwoD(vector<vector<double>> &Hresults, int volID, double current_Density, GetMesh meshdata, vector<vector<double>> gaussPointsData, string filePath);
	vector<vector<double>>  integrateLine(double current,int volID,GetMesh meshdata, vector<vector<double>> gaussPointsData, string path);
	~BiotSavart();
};
#endif


class CheckPair {
	int i;
public:
	CheckPair(int j) : i(j) { }
	bool operator()(std::pair<int, double> p)
	{
		return i == p.first;
	}

};