#ifndef DRUDE_LORENTZ_H
#define DRUDE_LORENTZ_H

#include "Structure.h"
#include "Fields.h"

struct DL_material
{
	double *JDx, *JDy, *JLx, *JLy;
	double *JLx_1, *JLy_1, Ex_1, Ey_1;
	int ix, jy;
	double temp1x, temp2x, temp1y, temp2y;
};

class Drude_Lorentz
{
public:
	Drude_Lorentz(Structure &sS, char *fileName);
	~Drude_Lorentz();
	double sigma, epsInf, *omega_D, *gamma_D, *delta_eps_L, *omega_L, *delta_L;
	double **Jx, **Jy;
	int num_of_D_poles, num_of_L_poles, num_of_DLM_gridcells;
	DL_material *DLM;
	void EFieldUpdate(Field &fF, Structure &sS);

private:
	double *beta, *kappa, *alpha, *xi, *gamma_L;
	double C1, C2, C3, *C4, *C5, *C6;
	double denominator, temp1, temp2;

};



#endif
