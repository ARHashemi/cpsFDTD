#ifndef DRUDE_H
#define DRUDE_H

#include "Structure.h"
#include "Fields.h"

class Drude
{
public:
	Drude(Structure &sS, char *fileName);
	~Drude();
	double epsInf, omegaP, sigma, gamma;
	double **Jx, **Jy;
	void EFieldUpdate(Field &fF, Structure &sS);

private:
	double beta, kappa, C1Drude, C2Drude, C3Drude;
	
};

#endif

