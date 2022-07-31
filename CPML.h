#ifndef CPML_H
#define CPML_H

#include "Structure.h"
#include "Fields.h"

class CPML
{
public:
	CPML(Structure &sS, char *fileName);
	~CPML();
	void InitializeCPML(Structure &sS);
	void UpdateCPMLEFields(Field &fF);
	void UpdateCPMLHFields(Field &fF);


private:
	double ***QHxz1, ***QHxz2;
	double ***QHxy1, ***QHxy2;
	double ***QHyz1, ***QHyz2;
	double ***QHyx1, ***QHyx2;
	double ***QHzx1, ***QHzx2;
	double ***QHzy1, ***QHzy2;

	double ***QExz1, ***QExz2;
	double ***QExy1, ***QExy2;
	double ***QEyz1, ***QEyz2;
	double ***QEyx1, ***QEyx2;
	double ***QEzx1, ***QEzx2;
	double ***QEzy1, ***QEzy2;

	double *bxE1, *cxE1, *byE1, *cyE1, *bzE1, *czE1, *bxH1, *cxH1, *byH1, *cyH1, *bzH1, *czH1;
	double *bxE2, *cxE2, *byE2, *cyE2, *bzE2, *czE2, *bxH2, *cxH2, *byH2, *cyH2, *bzH2, *czH2;
	double cpmlE, cpmlH;

	void E_pml_Coef(double* &kappa, double* &bcoef, double* &ccoef, bool firstend, int npml, double dt, double del, double m, double ma, double sigmax, double kapmax, double alpmax);
	void H_pml_Coef(double* &kappa, double* &bcoef, double* &ccoef, bool firstend, int npml, double dt, double del, double m, double ma, double sigmax, double kapmax, double alpmax);
	int Nx, Ny, Nz;
	int Npmlx1, Npmly1, Npmlz1;
	int Npmlx2, Npmly2, Npmlz2;
	double m, malp, kapmax, alpmax, sigmax, sig_ratio, *kappa;
//	double *sigmaEx, *kappaEx, *alphaEx;
//	double *sigmaEy, *kappaEy, *alphaEy;
//	double *sigmaEz, *kappaEz, *alphaEz;
//	double *sigmaHx, *kappaHx, *alphaHx;
//	double *sigmaHy, *kappaHy, *alphaHy;
//	double *sigmaHz, *kappaHz, *alphaHz;
};



#endif
