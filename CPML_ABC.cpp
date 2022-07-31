#include "CPML_ABC.h"
#include "stdafx.h"
#include "Constants.h"
#include "Structure.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <omp.h>
using namespace std;

CPML_ABC::CPML_ABC(Structure &sS, char *fileName)
{
	cout << "Constructing CPML Boundaries...";
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> Npmlx;
	iss >> Npmly;
	iss >> Npmlz;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> m;
	iss >> malp;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> kapmax;
	iss >> alpmax;
	sigmax = (0.8 * (m + 1) / (sS.dx * sqrt(mu0 / eps0)));
	//cout << m << "\t" << kapmax << endl;
	xN = sS.Nx;
	yN = sS.Ny;
	zN = sS.Nz;
	cout << " Done." << endl;

	cout << "Initializing CPML Coefficients...";
	Num_of_PML_points = 2*(Npmlx+1)*yN*zN + 2*xN*(Npmly+1)*zN + 2*xN*yN*(Npmlz+1) + 1;
	PMLpoint = new PML_point[Num_of_PML_points];
	int i, j, k, ii = 0;
	for (i = 1; i <= (Npmlx); ++i){
		for (j = 1; j<yN; ++j){
			for (k = 1; k<zN; ++k){
                PMLpoint[ii].i = i;
                PMLpoint[ii].j = j;
                PMLpoint[ii].k = k;
                PMLpoint[ii].QE1 = 0.0;
                PMLpoint[ii].QE2 = 0.0;
                PMLpoint[ii].QH1 = 0.0;
                PMLpoint[ii].QH2 = 0.0;
                PMLpoint[ii].grade = i;

                PMLpoint[ii+1].i = xN - i;
                PMLpoint[ii+1].j = j;
                PMLpoint[ii+1].k = k;
                PMLpoint[ii+1].QE1 = 0.0;
                PMLpoint[ii+1].QE2 = 0.0;
                PMLpoint[ii+1].QH1 = 0.0;
                PMLpoint[ii+1].QH2 = 0.0;
                PMLpoint[ii+1].grade = i;

                ii = ii + 2;
			}
		}
	}
	iE1 = ii;
	for (i = 1; i < xN; ++i){
		for (j = 1; j <= (Npmly); ++j){
			for (k = 1; k<zN; ++k){
                PMLpoint[ii].i = i;
                PMLpoint[ii].j = j;
                PMLpoint[ii].k = k;
                PMLpoint[ii].QE1 = 0.0;
                PMLpoint[ii].QE2 = 0.0;
                PMLpoint[ii].QH1 = 0.0;
                PMLpoint[ii].QH2 = 0.0;
                PMLpoint[ii].grade = j;

                PMLpoint[ii+1].i = i;
                PMLpoint[ii+1].j = yN - j;
                PMLpoint[ii+1].k = k;
                PMLpoint[ii+1].QE1 = 0.0;
                PMLpoint[ii+1].QE2 = 0.0;
                PMLpoint[ii+1].QH1 = 0.0;
                PMLpoint[ii+1].QH2 = 0.0;
                PMLpoint[ii+1].grade = j;

                ii = ii + 2;
			}
		}
	}
	iE2 = ii;
	for (i = 1; i < xN; ++i){
		for (j = 1; j < yN ; ++j){
			for (k = 1; k<=(Npmlz); ++k){
                PMLpoint[ii].i = i;
                PMLpoint[ii].j = j;
                PMLpoint[ii].k = k;
                PMLpoint[ii].QE1 = 0.0;
                PMLpoint[ii].QE2 = 0.0;
                PMLpoint[ii].QH1 = 0.0;
                PMLpoint[ii].QH2 = 0.0;
                PMLpoint[ii].grade = k;

                PMLpoint[ii+1].i = i;
                PMLpoint[ii+1].j = j;
                PMLpoint[ii+1].k = zN - k;
                PMLpoint[ii+1].QE1 = 0.0;
                PMLpoint[ii+1].QE2 = 0.0;
                PMLpoint[ii+1].QH1 = 0.0;
                PMLpoint[ii+1].QH2 = 0.0;
                PMLpoint[ii+1].grade = k;

                ii = ii + 2;
			}
		}
	}
	iE3 = ii;

	sigmaEx = new double[Npmlx + 1];
	kappaEx = new double[Npmlx + 1];
	alphaEx = new double[Npmlx + 1];
	bxE = new double[Npmlx + 1];
	cxE = new double[Npmlx + 1];
	sigmaEy = new double[Npmly + 1];
	kappaEy = new double[Npmly + 1];
	alphaEy = new double[Npmly + 1];
	byE = new double[Npmly + 1];
	cyE = new double[Npmly + 1];
	sigmaHx = new double[Npmlx + 1];
	kappaHx = new double[Npmlx + 1];
	alphaHx = new double[Npmlx + 1];
	bxH = new double[Npmlx + 1];
	cxH = new double[Npmlx + 1];
	sigmaHy = new double[Npmly + 1];
	kappaHy = new double[Npmly + 1];
	alphaHy = new double[Npmly + 1];
	byH = new double[Npmly + 1];
	cyH = new double[Npmly + 1];
	sigmaEz = new double[Npmlz + 1];
	kappaEz = new double[Npmlz + 1];
	alphaEz = new double[Npmlz + 1];
	bzE = new double[Npmlz + 1];
	czE = new double[Npmlz + 1];
	sigmaHz = new double[Npmlz + 1];
	kappaHz = new double[Npmlz + 1];
	alphaHz = new double[Npmlz + 1];
	bzH = new double[Npmlz + 1];
	czH = new double[Npmlz + 1];

	double polyScaleSK, polyScaleA, temp, dtotau;
	//x directed boundary for E field
	for (i = 1; i <= (Npmlx); ++i){
		polyScaleSK = pow(((double)i / (double)Npmlx), m);
		polyScaleA = pow(((double)(Npmlx - i) / (double)Npmlx), malp);
		sigmaEx[i] = sigmax * polyScaleSK;
		kappaEx[i] = 1.0 + (kapmax - 1.0) * polyScaleSK;
		alphaEx[i] = alpmax * polyScaleA;

		temp = kappaEx[i] * alphaEx[i] + sigmaEx[i];
		dtotau = temp * sS.dt / (kappaEx[i] * eps0);
		bxE[i] = exp(-dtotau);
		if (temp == 0.0){
			cxE[i] = 0.0;
		}
		else{
			cxE[i] = sigmaEx[i] * (1.0 - bxE[i]) / (kappaEx[i] * temp);
		}
		sS.cEyz[i] = 1.0 / (kappaEx[i] * sS.dx);
		sS.cEyz[sS.Nx - i] = 1.0 / (kappaEx[i] * sS.dx);
		sS.cEzy[i] = 1.0 / (kappaEx[i] * sS.dx);
		sS.cEzy[sS.Nx - i] = 1.0 / (kappaEx[i] * sS.dx);
	}
	//y directed boundary for E field
	for (j = 1; j <= (Npmly); ++j){
		polyScaleSK = pow(((double)j / (double)Npmly), m);
		polyScaleA = pow(((double)(Npmly - j) / (double)Npmly), malp);
		sigmaEy[j] = sigmax * polyScaleSK;
		kappaEy[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
		alphaEy[j] = alpmax * polyScaleA;

		temp = kappaEy[j] * alphaEy[j] + sigmaEy[j];
		dtotau = temp * sS.dt / (kappaEy[j] * eps0);
		byE[j] = exp(-dtotau);
		if (temp == 0.0){
			cyE[j] = 0.0;
		}
		else{
			cyE[j] = sigmaEy[j] * (1.0 - byE[j]) / (kappaEy[j] * temp);
		}
		sS.cExz[j] = 1.0 / (kappaEy[j] * sS.dy);
		sS.cExz[sS.Ny - j] = 1.0 / (kappaEy[j] * sS.dy);
		sS.cEzx[j] = 1.0 / (kappaEy[j] * sS.dy);
		sS.cEzx[sS.Ny - j] = 1.0 / (kappaEy[j] * sS.dy);
	}
	//z directed boundary for E field
	for (j = 1; j <= (Npmlz); ++j){
		polyScaleSK = pow(((double)j / (double)Npmlz), m);
		polyScaleA = pow(((double)(Npmlz - j) / (double)Npmlz), malp);
		sigmaEz[j] = sigmax * polyScaleSK;
		kappaEz[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
		alphaEz[j] = alpmax * polyScaleA;

		temp = kappaEz[j] * alphaEz[j] + sigmaEz[j];
		dtotau = temp * sS.dt / (kappaEz[j] * eps0);
		bzE[j] = exp(-dtotau);
		if (temp == 0.0){
			czE[j] = 0.0;
		}
		else{
			czE[j] = sigmaEz[j] * (1.0 - bzE[j]) / (kappaEz[j] * temp);
		}
		sS.cExy[j] = 1.0 / (kappaEz[j] * sS.dz);
		sS.cExy[sS.Nz - j] = 1.0 / (kappaEz[j] * sS.dz);
		sS.cEyx[j] = 1.0 / (kappaEz[j] * sS.dz);
		sS.cEyx[sS.Nz - j] = 1.0 / (kappaEz[j] * sS.dz);
	}
	//x directed boundary for H field
	for (i = 1; i < (Npmlx); ++i){
		polyScaleSK = pow(((double)(i + 0.5) / (double)Npmlx), m);
		polyScaleA = pow(((double)(Npmlx - i - 0.5) / (double)Npmlx), malp);
		sigmaHx[i] = sigmax * polyScaleSK;
		kappaHx[i] = 1.0 + (kapmax - 1.0) * polyScaleSK;
		alphaHx[i] = alpmax * polyScaleA;

		temp = kappaHx[i] * alphaHx[i] + sigmaHx[i];
		dtotau = temp * sS.dt / (kappaHx[i] * eps0);
		bxH[i] = exp(-dtotau);
		if (temp == 0.0){
			cxH[i] = 0.0;
		}
		else{
			cxH[i] = sigmaHx[i] * (1.0 - bxH[i]) / (kappaHx[i] * temp);
		}
		sS.cHzy[i] = 1.0 / (kappaHx[i] * sS.dx);
		sS.cHzy[sS.Nx - i] = 1.0 / (kappaHx[i] * sS.dx);
		sS.cHyz[i] = 1.0 / (kappaHx[i] * sS.dx);
		sS.cHyz[sS.Nx - i] = 1.0 / (kappaHx[i] * sS.dx);
	}
//	i = Npmlx-1;
//	polyScaleSK = pow(((Npmlx + 0.5 - i) / Npmlx), m);
//	polyScaleA = pow(((i - 0.5) / Npmlx), malp);
//	sigmaHx[i] = sigmax * polyScaleSK;
//	kappaHx[i] = 1.0 + (kapmax - 1.0) * polyScaleSK;
//	alphaHx[i] = alpmax * polyScaleA;
	//y directed boundary for H field
	for (j = 1; j < (Npmly); ++j){
		polyScaleSK = pow(((double)(j + 0.5) / (double)Npmly), m);
		polyScaleA = pow(((double)(Npmly - j - 0.5) / (double)Npmly), malp);
		sigmaHy[j] = sigmax * polyScaleSK;
		kappaHy[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
		alphaHy[j] = alpmax * polyScaleA;

		temp = kappaHy[j] * alphaHy[j] + sigmaHy[j];
		dtotau = temp * sS.dt / (kappaHy[j] * eps0);
		byH[j] = exp(-dtotau);
		if (temp == 0.0){
			cyH[j] = 0.0;
		}
		else{
			cyH[j] = sigmaHy[j] * (1.0 - byH[j]) / (kappaHy[j] * temp);
		}
		sS.cHzx[j] = 1.0 / (kappaHy[j] * sS.dy);
		sS.cHzx[sS.Ny - j] = 1.0 / (kappaHy[j] * sS.dy);
		sS.cHxz[j] = 1.0 / (kappaHy[j] * sS.dy);
		sS.cHxz[sS.Ny - j] = 1.0 / (kappaHy[j] * sS.dy);
	}
//	j = Npmly-1;
//	polyScaleSK = pow(((Npmly + 0.5 - j) / Npmly), m);
//	polyScaleA = pow(((j - 0.5) / Npmly), malp);
//	sigmaHy[j] = sigmax * polyScaleSK;
//	kappaHy[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
//	alphaHy[j] = alpmax * polyScaleA;
	//z directed boundary for E field
	for (j = 1; j < (Npmlz); ++j){
		polyScaleSK = pow(((double)(j + 0.5) / (double)Npmlz), m);
		polyScaleA = pow(((double)(Npmlz - j - 0.5) / (double)Npmlz), malp);
		sigmaHz[j] = sigmax * polyScaleSK;
		kappaHz[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
		alphaHz[j] = alpmax * polyScaleA;

		temp = kappaHz[j] * alphaHz[j] + sigmaHz[j];
		dtotau = temp * sS.dt / (kappaHz[j] * eps0);
		byH[j] = exp(-dtotau);
		if (temp == 0.0){
			czH[j] = 0.0;
		}
		else{
			czH[j] = sigmaHz[j] * (1.0 - bzH[j]) / (kappaHz[j] * temp);
		}
		sS.cHxy[j] = 1.0 / (kappaHz[j] * sS.dz);
		sS.cHxy[sS.Nz - j] = 1.0 / (kappaHz[j] * sS.dz);
		sS.cHyx[j] = 1.0 / (kappaHz[j] * sS.dz);
		sS.cHyx[sS.Nz - j] = 1.0 / (kappaHz[j] * sS.dz);
	}
//	j = Npmlz-1;
//	polyScaleSK = pow(((Npmlz + 0.5 - j) / Npmlz), m);
//	polyScaleA = pow(((j - 0.5) / Npmlz), malp);
//	sigmaHz[j] = sigmax * polyScaleSK;
//	kappaHz[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
//	alphaHz[j] = alpmax * polyScaleA;

	cpmlE = sS.bE[0];
	cpmlH = sS.bH[0];
	cout << " Done." << endl;
}

CPML_ABC::~CPML_ABC()
{
	//dtor
}


void CPML_ABC::UpdateCPMLEFields(Field &fF)
{
	int ii;
	//#pragma omp parallel for default(shared) private(ii)
	for (ii = 0; ii <iE1; ++ii){
		PMLpoint[ii].QH1 = bxE[PMLpoint[ii].grade] * PMLpoint[ii].QH1 + cxE[PMLpoint[ii].grade] * (fF.Hz[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - fF.Hz[PMLpoint[ii].i - 1][PMLpoint[ii].j][PMLpoint[ii].k]);
		fF.Ey[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Ey[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - PMLpoint[ii].QH1 * cpmlE;
		PMLpoint[ii].QH2 = bxE[PMLpoint[ii].grade] * PMLpoint[ii].QH2 + cxE[PMLpoint[ii].grade] * (fF.Hy[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - fF.Hy[PMLpoint[ii].i - 1][PMLpoint[ii].j][PMLpoint[ii].k]);
		fF.Ez[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Ez[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] + PMLpoint[ii].QH2 * cpmlE;
	}

	//#pragma omp parallel for default(shared) private(ii)
	for (ii = iE1; ii <iE2; ++ii){
		PMLpoint[ii].QH1 = byE[PMLpoint[ii].grade] * PMLpoint[ii].QH1 + cyE[PMLpoint[ii].grade] * (fF.Hx[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - fF.Hx[PMLpoint[ii].i][PMLpoint[ii].j - 1][PMLpoint[ii].k]);
		fF.Ez[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Ez[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - PMLpoint[ii].QH1 * cpmlE;
		PMLpoint[ii].QH2 = byE[PMLpoint[ii].grade] * PMLpoint[ii].QH2 + cyE[PMLpoint[ii].grade] * (fF.Hz[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - fF.Hz[PMLpoint[ii].i][PMLpoint[ii].j - 1][PMLpoint[ii].k]);
		fF.Ex[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Ex[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] + PMLpoint[ii].QH2 * cpmlE;
	}

	//#pragma omp parallel for default(shared) private(ii)
	for (ii = iE2; ii <iE3; ++ii){
		PMLpoint[ii].QH1 = bzE[PMLpoint[ii].grade] * PMLpoint[ii].QH1 + czE[PMLpoint[ii].grade] * (fF.Hy[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - fF.Hy[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k - 1]);
		fF.Ex[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Ex[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - PMLpoint[ii].QH1 * cpmlE;
		PMLpoint[ii].QH2 = bzE[PMLpoint[ii].grade] * PMLpoint[ii].QH2 + czE[PMLpoint[ii].grade] * (fF.Hx[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - fF.Hx[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k - 1]);
		fF.Ey[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Ey[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] + PMLpoint[ii].QH2 * cpmlE;
	}

}

void CPML_ABC::UpdateCPMLHFields(Field &fF)
{
	int ii;
	//#pragma omp parallel for default(shared) private(ii)
	for (ii = 0; ii <iE1; ++ii){
		PMLpoint[ii].QE1 = bxH[PMLpoint[ii].grade-1] * PMLpoint[ii].QE1 + cxH[PMLpoint[ii].grade-1] * (fF.Ez[PMLpoint[ii].i + 1][PMLpoint[ii].j][PMLpoint[ii].k] - fF.Ez[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k]);
		fF.Hy[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Hy[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] + PMLpoint[ii].QE1 * cpmlH;
		PMLpoint[ii].QE2 = bxH[PMLpoint[ii].grade-1] * PMLpoint[ii].QE2 + cxH[PMLpoint[ii].grade-1] * (fF.Ey[PMLpoint[ii].i + 1][PMLpoint[ii].j][PMLpoint[ii].k] - fF.Ey[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k]);
		fF.Hz[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Hz[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - PMLpoint[ii].QE2 * cpmlH;
	}

	//#pragma omp parallel for default(shared) private(ii)
	for (ii = iE1; ii <iE2; ++ii){
		PMLpoint[ii].QE1 = byH[PMLpoint[ii].grade-1] * PMLpoint[ii].QE1 + cyH[PMLpoint[ii].grade-1] * (fF.Ex[PMLpoint[ii].i][PMLpoint[ii].j + 1][PMLpoint[ii].k] - fF.Ex[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k]);
		fF.Hz[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Hz[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] + PMLpoint[ii].QE1 * cpmlH;
		PMLpoint[ii].QE2 = byH[PMLpoint[ii].grade-1] * PMLpoint[ii].QE2 + cyH[PMLpoint[ii].grade-1] * (fF.Ez[PMLpoint[ii].i][PMLpoint[ii].j + 1][PMLpoint[ii].k] - fF.Ez[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k]);
		fF.Hx[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Hx[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - PMLpoint[ii].QE2 * cpmlH;
	}

	//#pragma omp parallel for default(shared) private(ii)
	for (ii = iE2; ii <iE3; ++ii){
		PMLpoint[ii].QE1 = bzH[PMLpoint[ii].grade-1] * PMLpoint[ii].QE1 + czH[PMLpoint[ii].grade-1] * (fF.Ey[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k + 1] - fF.Ey[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k]);
		fF.Hx[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Hx[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] + PMLpoint[ii].QE1 * cpmlH;
		PMLpoint[ii].QE2 = bzH[PMLpoint[ii].grade-1] * PMLpoint[ii].QE2 + czH[PMLpoint[ii].grade-1] * (fF.Ex[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k + 1] - fF.Ex[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k]);
		fF.Hy[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] = fF.Hy[PMLpoint[ii].i][PMLpoint[ii].j][PMLpoint[ii].k] - PMLpoint[ii].QE2 * cpmlH;
	}

}


//for (i = 1; i <= (Npmlx); ++i){
//		polyScaleSK = pow(((Npmlx - i) / (Npmlx-1)), m);
//		polyScaleA = pow(((i - 1.0) / (Npmlx-1)), malp);
//		sigmaEx[i] = sigmax * polyScaleSK;
//		kappaEx[i] = 1.0 + (kapmax - 1.0) * polyScaleSK;
//		alphaEx[i] = alpmax * polyScaleA;
//
//		temp = kappaEx[i] * alphaEx[i] + sigmaEx[i];
//		dtotau = temp * sS.dt / (kappaEx[i] * eps0);
//		bxE[i] = exp(-dtotau);
//		if (temp == 0.0){
//			cxE[i] = 0.0;
//		}
//		else{
//			cxE[i] = sigmaEx[i] * (bxE[i] - 1.0) / (kappaEx[i] * temp * sS.dx);
//		}
//		sS.cEyz[i] = 1.0 / (kappaEx[i] * sS.dx);
//		sS.cEyz[sS.Nx - i] = 1.0 / (kappaEx[i] * sS.dx);
//		sS.cEzy[i] = 1.0 / (kappaEx[i] * sS.dx);
//		sS.cEzy[sS.Nx - i] = 1.0 / (kappaEx[i] * sS.dx);
//	}
//	//y directed boundary for E field
//	for (j = 1; j <= (Npmly); ++j){
//		polyScaleSK = pow(((Npmly - j) / (Npmly-1)), m);
//		polyScaleA = pow(((j - 1.0) / (Npmly-1)), malp);
//		sigmaEy[j] = sigmax * polyScaleSK;
//		kappaEy[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
//		alphaEy[j] = alpmax * polyScaleA;
//
//		temp = kappaEy[j] * alphaEy[j] + sigmaEy[j];
//		dtotau = temp * sS.dt / (kappaEy[j] * eps0);
//		byE[j] = exp(-dtotau);
//		if (temp == 0.0){
//			cyE[j] = 0.0;
//		}
//		else{
//			cyE[j] = sigmaEy[j] * (byE[j] - 1.0) / (kappaEy[j] * temp * sS.dy);
//		}
//		sS.cExz[j] = 1.0 / (kappaEy[j] * sS.dy);
//		sS.cExz[sS.Ny - j] = 1.0 / (kappaEy[j] * sS.dy);
//		sS.cEzx[j] = 1.0 / (kappaEy[j] * sS.dy);
//		sS.cEzx[sS.Ny - j] = 1.0 / (kappaEy[j] * sS.dy);
//	}
//	//z directed boundary for E field
//	for (j = 1; j <= (Npmlz); ++j){
//		polyScaleSK = pow(((Npmlz - j) / (Npmlz-1)), m);
//		polyScaleA = pow(((j - 1.0) / (Npmlz-1)), malp);
//		sigmaEz[j] = sigmax * polyScaleSK;
//		kappaEz[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
//		alphaEz[j] = alpmax * polyScaleA;
//
//		temp = kappaEz[j] * alphaEz[j] + sigmaEz[j];
//		dtotau = temp * sS.dt / (kappaEz[j] * eps0);
//		bzE[j] = exp(-dtotau);
//		if (temp == 0.0){
//			czE[j] = 0.0;
//		}
//		else{
//			czE[j] = sigmaEz[j] * (bzE[j] - 1.0) / (kappaEz[j] * temp * sS.dz);
//		}
//		sS.cExy[j] = 1.0 / (kappaEz[j] * sS.dz);
//		sS.cExy[sS.Nz - j] = 1.0 / (kappaEz[j] * sS.dz);
//		sS.cEyx[j] = 1.0 / (kappaEz[j] * sS.dz);
//		sS.cEyx[sS.Nz - j] = 1.0 / (kappaEz[j] * sS.dz);
//	}
//	//x directed boundary for H field
//	for (i = 1; i < (Npmlx); ++i){
//		polyScaleSK = pow(((Npmlx + 0.5 - i) / (Npmlx-1)), m);
//		polyScaleA = pow(((i - 0.5) / (Npmlx-1)), malp);
//		sigmaHx[i] = sigmax * polyScaleSK;
//		kappaHx[i] = 1.0 + (kapmax - 1.0) * polyScaleSK;
//		alphaHx[i] = alpmax * polyScaleA;
//
//		temp = kappaHx[i] * alphaHx[i] + sigmaHx[i];
//		dtotau = temp * sS.dt / (kappaHx[i] * eps0);
//		bxH[i] = exp(-dtotau);
//		if (temp == 0.0){
//			cxH[i] = 0.0;
//		}
//		else{
//			cxH[i] = sigmaHx[i] * (bxH[i] - 1.0) / (kappaHx[i] * temp * sS.dx);
//		}
//		sS.cHzy[i] = 1.0 / (kappaHx[i] * sS.dx);
//		sS.cHzy[sS.Nx - i] = 1.0 / (kappaHx[i] * sS.dx);
//		sS.cHyz[i] = 1.0 / (kappaHx[i] * sS.dx);
//		sS.cHyz[sS.Nx - i] = 1.0 / (kappaHx[i] * sS.dx);
//	}
////	i = Npmlx-1;
////	polyScaleSK = pow(((Npmlx + 0.5 - i) / Npmlx), m);
////	polyScaleA = pow(((i - 0.5) / Npmlx), malp);
////	sigmaHx[i] = sigmax * polyScaleSK;
////	kappaHx[i] = 1.0 + (kapmax - 1.0) * polyScaleSK;
////	alphaHx[i] = alpmax * polyScaleA;
//	//y directed boundary for H field
//	for (j = 1; j < (Npmly); ++j){
//		polyScaleSK = pow(((Npmly + 0.5 - j) / (Npmly-1)), m);
//		polyScaleA = pow(((j - 0.5) / (Npmly-1)), malp);
//		sigmaHy[j] = sigmax * polyScaleSK;
//		kappaHy[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
//		alphaHy[j] = alpmax * polyScaleA;
//
//		temp = kappaHy[j] * alphaHy[j] + sigmaHy[j];
//		dtotau = temp * sS.dt / (kappaHy[j] * eps0);
//		byH[j] = exp(-dtotau);
//		if (temp == 0.0){
//			cyH[j] = 0.0;
//		}
//		else{
//			cyH[j] = sigmaHy[j] * (byH[j] - 1.0) / (kappaHy[j] * temp * sS.dy);
//		}
//		sS.cHzx[j] = 1.0 / (kappaHy[j] * sS.dy);
//		sS.cHzx[sS.Ny - j] = 1.0 / (kappaHy[j] * sS.dy);
//		sS.cHxz[j] = 1.0 / (kappaHy[j] * sS.dy);
//		sS.cHxz[sS.Ny - j] = 1.0 / (kappaHy[j] * sS.dy);
//	}
////	j = Npmly-1;
////	polyScaleSK = pow(((Npmly + 0.5 - j) / Npmly), m);
////	polyScaleA = pow(((j - 0.5) / Npmly), malp);
////	sigmaHy[j] = sigmax * polyScaleSK;
////	kappaHy[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
////	alphaHy[j] = alpmax * polyScaleA;
//	//z directed boundary for E field
//	for (j = 1; j < (Npmlz); ++j){
//		polyScaleSK = pow(((Npmlz + 0.5 - j) / (Npmlz-1)), m);
//		polyScaleA = pow(((j - 0.5) / (Npmlz-1)), malp);
//		sigmaHz[j] = sigmax * polyScaleSK;
//		kappaHz[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
//		alphaHz[j] = alpmax * polyScaleA;
//
//		temp = kappaHz[j] * alphaHz[j] + sigmaHz[j];
//		dtotau = temp * sS.dt / (kappaHz[j] * eps0);
//		byH[j] = exp(-dtotau);
//		if (temp == 0.0){
//			czH[j] = 0.0;
//		}
//		else{
//			czH[j] = sigmaHz[j] * (bzH[j] - 1.0) / (kappaHz[j] * temp * sS.dz);
//		}
//		sS.cHxy[j] = 1.0 / (kappaHz[j] * sS.dz);
//		sS.cHxy[sS.Nz - j] = 1.0 / (kappaHz[j] * sS.dz);
//		sS.cHyx[j] = 1.0 / (kappaHz[j] * sS.dz);
//		sS.cHyx[sS.Nz - j] = 1.0 / (kappaHz[j] * sS.dz);
//	}
////	j = Npmlz-1;
////	polyScaleSK = pow(((Npmlz + 0.5 - j) / Npmlz), m);
////	polyScaleA = pow(((j - 0.5) / Npmlz), malp);
////	sigmaHz[j] = sigmax * polyScaleSK;
////	kappaHz[j] = 1.0 + (kapmax - 1.0) * polyScaleSK;
////	alphaHz[j] = alpmax * polyScaleA;
