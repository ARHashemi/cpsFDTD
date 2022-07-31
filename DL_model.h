#ifndef DL_MODEL_H
#define DL_MODEL_H

#include "Structure.h"
#include "Fields.h"

struct DL_material
{
	double Jx, Jy, Jz, Px, Py, Pz;
	double Px_1, Py_1, Pz_1, Ex, Ey, Ez, Ex_1, Ey_1, Ez_1;
	int ix, jy, kz;
};

class DL_model
{
public:
	DL_model(Structure &sS, char *fileName);
	~DL_model();
	double sigma, epsInf, omega_D, gamma_D, delta_eps_L, omega_L, gamma_L;
	int num_of_D_poles, num_of_L_poles, num_of_DLM_gridcells;
	DL_material *DLM;
	void EFieldUpdate(Field &fF, Structure &sS);

private:
	double alpha, beta, OMEGA, zeta, tau, rho;
	double C1, C2, C3, C4, C5, C6, C7, C8, C9;
	double denominator;//, temp1, temp2, temp3;
};

#endif
