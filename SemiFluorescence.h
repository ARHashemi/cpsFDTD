#ifndef SEMIFLUORESCENCE_H
#define SEMIFLUORESCENCE_H

#include "Structure.h"
#include "Fields.h"



class SemiFluorescence
{
public:
	SemiFluorescence(Structure &sS, char *fileName);
	~SemiFluorescence();
	void UpdateEmitters(Field &fF);
	double lambda, nu, T, EmissionThershold, EmissionAmp;
	int nT, nw, nd;
	int NumSources, i, i1, i2, j1, j2, idif, jdif;
	int **location;
	double *Amps;
	
private:
	double dt, pinudt;
};

#endif
