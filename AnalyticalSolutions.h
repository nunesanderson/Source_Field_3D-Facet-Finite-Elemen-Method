#include <string>
#include<vector>
using namespace std;
using std::string;

/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include "Matrix.h"


#ifndef _ANALYTICAL_SOLUTIONS_INCLUDED_
#define _ANALYTICAL_SOLUTIONS_INCLUDED_

class BiotSavartAnalyt
{
public:
	vector<double> getHcoil(double radius,vector<double> listZPos, double len,double current,double turns);
	vector<double> getHCurrentLoop(double radius, vector<double> listZPos, double current);

private:

};



#endif