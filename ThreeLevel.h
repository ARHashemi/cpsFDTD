#ifndef THREELEVELLASER_H
#define THREELEVELLASER_H

#include "Structure.h"
#include "Fields.h"

class ThreeLevel
{
public:
	ThreeLevel(Structure &sS, char *fileName);
	~ThreeLevel();
	void UpdatePolarizationDensities(Field &fF);
	void UpdateElectricField(Field &fF);
	void UpdateStatesPopulations(Field &fF);

	int **location, NumSources;
	int i1, i2, j1, j2, idif, jdif;
	double *N2, *N1, *N0;
	double N_0, N_1;
	double *Pax, *Pay, *Paz, *Pbx, *Pby, *Pbz;
	double gamma_a, gamma_b, lambda_a, lambda_b, omega_a, omega_b;
	double t20, t21, t10;
	char distri[40];

private:
	int i, j, ii;
	double *Pax_1, *Pay_1, *Paz_1, *Pbx_1, *Pby_1, *Pbz_1;
	double *Pax_2, *Pay_2, *Paz_2, *Pbx_2, *Pby_2, *Pbz_2;
	double *Ex_1, *Ey_1, *Ez_1;
	double *N2_1, *N1_1;
	double CPa1, CPa2, CPa3, CPb1, CPb2, CPb3, CEP;
	double CN2N2, CN2E, CN1N1, CN1N2, CN1E, CN0N2, CN0N1, CN0Ea, CN0Eb;
};

#endif