#ifndef FOURLEVELLASER_H
#define FOURLEVELLASER_H

#include "Structure.h"
#include "Fields.h"

class FourLevelLaser
{
public:
	FourLevelLaser(Structure &sS, char *fileName);
	~FourLevelLaser();
	void UpdatePolarizationDensities(Field &fF);
	void UpdateElectricField(Field &fF);
	void UpdateStatesPopulations(Field &fF);

	int **location, NumSources;
	int i1, i2, j1, j2, idif, jdif;
	double *N3, *N2, *N1, *N0;
	double N_0, N_1;
	double *Pax, *Pay, *Paz, *Pbx, *Pby, *Pbz;
	double gamma_a, gamma_b, lambda_a, lambda_b, omega_a, omega_b;
	double t30, t32, t21, t10;
	char distri[40];
	char const *random_ditri = "random";
	
private:
	double CPa1, CPa2, CPa3, CPb1, CPb2, CPb3, CEP;
	double CN31, CN32, CN33;
	double CN21, CN22, CN23;
	double CN11, CN12, CN13;
	double CN01, CN02, CN03;
	double zeta_a, zeta_b;
	int i, j, ii;
	double *Pax_1, *Pay_1, *Paz_1, *Pbx_1, *Pby_1, *Pbz_1;
	double *Pax_2, *Pay_2, *Paz_2, *Pbx_2, *Pby_2, *Pbz_2;
	double *Ex_1, *Ey_1, *Ez_1;

	double t30t32dt, t30t32_2, t30dt, t32dt, dt_2, homega_bdt;
	double t10t21dt, t10t21_2, t10dt, t21dt, t21_2, homega_adt;
	double t32t21dt, t32t21_2, t32_2;
	double t10t30dt, t10t30_2;

	double dtt32, dtt30, homega_b_2;
	double dtt21, dtt10, homega_a_2;
};

#endif

