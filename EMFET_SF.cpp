#include "stdafx.h"
#include <iostream>
#include<math.h>
using namespace std;
#include<vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include "Biot-Savart.h"
#include "ShapeFunctions.h"
#include "Matrix.h"
#include "Gmsh.h"
#include "Messages.h"
#include "Tests-Biot-Savart.h"


int main()
{
	Messages messages;

	double teste = pow(4, 0.5);
	
	TestsBiotSavart testBS;
	testBS;


	
	//*******************************
	// Mesh layer
	//*******************************
	
	//GetMesh mesh("C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\24_Test boundary\\model.msh");
	//int surf = 100;
	//int vol = 80;
	//vector<int> IDs = mesh.defineBoundary(mesh, surf, vol, 2);
	//mesh.writeMesh(mesh, "C:\\Anderson\\Pessoal\\01_Doutorado\\10_Testes\\24_Test boundary\\model2.msh", IDs);
	
	getchar();  
	return 0;
}

