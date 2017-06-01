#include <string>
#include<vector>
using namespace std;
using std::string;

/* ------------------------------------------------------------------------
Internal includes
---------------------------------------------------------------------------*/
#include "Matrix.h"
#include "Gmsh.h"

#ifndef _TESTS_BIOT_SAVART_INCLUDED_
#define _TESTS_BIOT_SAVART_INCLUDED_

class TestsBiotSavart
{
public:
	void volumeDomain();
	void currentLoop();
	void TwoD();
	void Atuador3D();
	void AtuadorHor();
	void AtuadorVert();
	void AtuadorVert_VS();
	void Teste3D();
	void H_discon();
	void Teste3D_core();
	void Subdomain();
	void AdaptiveProcessWinding(double fw);
	void AdaptiveProcess(double fw);
	void PerfectElectric();
	void PerfectMagnetic();
private:

};




#endif