#ifndef CPML_ABC_H
#define CPML_ABC_H

#include "Structure.h"
#include "Fields.h"

struct PML_point
{
	int i, j, k, grade;
	double QH1, QH2, QE1, QE2;
};

class CPML_ABC
{
public:
	CPML_ABC(Structure &sS, char *fileName);
	~CPML_ABC();
	void UpdateCPMLEFields(Field &fF);
	void UpdateCPMLHFields(Field &fF);

	PML_point *PMLpoint;

	double *bxE, *cxE, *byE, *cyE, *bzE, *czE, *bxH, *cxH, *byH, *cyH, *bzH, *czH;
	double cpmlE, cpmlH;
	double m, malp, kapmax, alpmax, sigmax;
	int Npmlx, Npmly, Npmlz;
	int Num_of_PML_points;
	int iE1, iE2, iE3;//, iH1, iH2, iH3;

private:
	int xN, yN, zN;
	double *sigmaEx, *kappaEx, *alphaEx;
	double *sigmaEy, *kappaEy, *alphaEy;
	double *sigmaEz, *kappaEz, *alphaEz;
	double *sigmaHx, *kappaHx, *alphaHx;
	double *sigmaHy, *kappaHy, *alphaHy;
	double *sigmaHz, *kappaHz, *alphaHz;
};

#endif // CPML_ABC_H
