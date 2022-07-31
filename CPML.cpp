#include "stdafx.h"
#include "Constants.h"
#include "Structure.h"
#include "CPML.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <omp.h>
using namespace std;

CPML::CPML(Structure &sS, char *fileName)
{
	cout << BOLDRED << ">" << RESET << "Constructing " << BOLDBLUE << "CPML " << RESET << "boundaries...";
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> Npmlx1;
	iss >> Npmlx2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> Npmly1;
	iss >> Npmly2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> Npmlz1;
	iss >> Npmlz2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> m;
	iss >> malp;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> kapmax;
	iss >> sig_ratio;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> alpmax;

	Nx = sS.Nx;
	Ny = sS.Ny;
	Nz = sS.Nz;
	cout << " Done." << endl;
}

CPML::~CPML()
{
	cout << "Deconstructing CPML...";
//	delete[] sigmaEx;
//	delete[] kappaEx;
//	delete[] alphaEx;
//	delete[] bxE;
//	delete[] cxE;
//	delete[] sigmaEy;
//	delete[] kappaEy;
//	delete[] alphaEy;
//	delete[] byE;
//	delete[] cyE;
//	delete[] sigmaHx;
//	delete[] kappaHx;
//	delete[] alphaHx;
//	delete[] bxH;
//	delete[] cxH;
//	delete[] sigmaHy;
//	delete[] kappaHy;
//	delete[] alphaHy;
//	delete[] byH;
//	delete[] cyH;
	int i, j;
	for (i = 0; i <= Npmlx1; i++){
		for(j = 0; j<=Npmly1; j++){
			delete[] QHxz1[i][j];
			delete[] QHxz2[i][j];
			delete[] QExy1[i][j];
			delete[] QExy2[i][j];
		}
		delete[] QHxz1[i];
		delete[] QHxz2[i];
		delete[] QExy1[i];
		delete[] QExy2[i];
	}
	delete[] QHxz1[Npmlx1 + 1];
	delete[] QHxz2[Npmlx1 + 1];
	delete[] QHxz1;
	delete[] QHxz2;
	delete[] QExy1;
	delete[] QExy2;
	QHxz1 = 0;
	QHxz2 = 0;
	QExy1 = 0;
	QExy2 = 0;
	for (i = 0; i <= Nx; i++){
		for(j = 0; j<=Npmly1; j++){
			delete[] QHyz1[i][j];
			delete[] QHyz2[i][j];
			delete[] QEyx1[i][j];
			delete[] QEyx2[i][j];
		}
		delete[] QHyz1[i];
		delete[] QHyz2[i];
		delete[] QEyx1[i];
		delete[] QEyx2[i];
	}
	delete[] QHyz1;
	delete[] QHyz2;
	delete[] QEyx1;
	delete[] QEyx2;
	QHyz1 = 0;
	QHyz2 = 0;
	QEyx1 = 0;
	QEyx2 = 0;
	cout << " Done." << endl;
}

void CPML::E_pml_Coef(double* &kappa, double* &bcoef, double* &ccoef, bool firstend, int npml, double dt, double del, double m, double ma, double sigmax, double kapmax, double alpmax)
{
	double polyscalesk, polyscalea, *sigma, *alpha, dtotau, temp;
	int i;
	sigma = new double[npml+1];
	alpha = new double[npml+1];
	if(firstend){
		for (i = 1; i <= (npml+1); ++i){
			polyscalesk = pow((((float)npml+1.0-i)/((float)npml)),m);
			polyscalea = pow((((float)i -1.0)/((float)npml)),ma);
			sigma[i-1] = sigmax * polyscalesk;
			kappa[i-1] = 1.0 + (kapmax-1.0) * polyscalesk;
			alpha[i-1] = alpmax * polyscalea;
		}
	}else{
		for (i = 1; i <= (npml+1); ++i){
			polyscalesk = pow((((float)i -1.0)/npml),m);
			polyscalea = pow((((float)npml+1.0-i)/npml),ma);
			sigma[i-1] = sigmax* polyscalesk;
			kappa[i-1] = 1.0 + (kapmax-1.0)* polyscalesk;
			alpha[i-1] = alpmax * polyscalea;
		}
	}
	for (i = 0; i <= npml; ++i){
//		temp = kappa[i]*alpha[i]+sigma[i];
//		dtotau = temp*dt/(kappa[i]*eps0);
//		bcoef[i] = exp(-dtotau);
//		if(temp == 0.0)
//			ccoef[i] = 0.0;
//		else
//			ccoef[i] = sigma[i]*(bcoef[i]-1.0)/(kappa[i]*temp*del);
		bcoef[i] = exp(-(sigma[i] / kappa[i] + alpha[i]) * dt / eps0);
		if ((sigma[i] == 0.0) && (alpha[i] == 0.0))
            ccoef[i] = 0.0;
		else
			ccoef[i] = sigma[i] * (bcoef[i] - 1.0) / (sigma[i] + kappa[i] * alpha[i]) / kappa[i] / del;//;


		//cout << ccoef[i] << "\t";
	}
	delete[] sigma;
	delete[] alpha;
	sigma = 0;
	alpha = 0;
}

void CPML::H_pml_Coef(double* &kappa, double* &bcoef, double* &ccoef, bool firstend, int npml, double dt, double del, double m, double ma, double sigmax, double kapmax, double alpmax)
{
	double polyscalesk, polyscalea, *sigma, *alpha, dtotau, temp;
	int i;
	sigma = new double[npml];
	alpha = new double[npml];
	if(firstend){
		for (i = 1; i <= npml; ++i){
			polyscalesk = pow(((npml-i+0.5)/npml),m);
			polyscalea = pow(((i - 0.5)/npml),ma);
			sigma[i-1] = sigmax * polyscalesk;
			kappa[i-1] = 1.0 + (kapmax-1.0) * polyscalesk;
			alpha[i-1] = alpmax * polyscalea;
		}

	}else{
		for(i = 1; i <= npml; ++i){
			polyscalesk = pow(((i - 0.5)/npml),m);
			polyscalea = pow(((npml-i+0.5)/npml),ma);
			sigma[i-1] = sigmax * polyscalesk;
			kappa[i-1] = 1.0 + (kapmax-1.0)* polyscalesk;
			alpha[i-1] = alpmax * polyscalea;
		}
	}

	for (i = 0; i < (npml-1); ++i){
//		temp = kappa[i]*alpha[i]+sigma[i];
//		dtotau = temp*dt/(kappa[i]*eps0);
//		bcoef[i] = exp(-dtotau);
//		if(temp == 0.0)
//			ccoef[i] = 0.0;
//		else
//			ccoef[i] = sigma[i]*(bcoef[i]-1.0)/(kappa[i]*temp*del);
		bcoef[i] = exp(-(sigma[i] / kappa[i] + alpha[i]) * dt / eps0);
		if ((sigma[i] == 0.0) && (alpha[i] == 0.0))
            ccoef[i] = 0.0;
		else
			ccoef[i] = sigma[i] * (bcoef[i] - 1.0) / (sigma[i] + kappa[i] * alpha[i]) / kappa[i] / del;//

	}
	delete[] sigma;
	delete[] alpha;
	sigma = 0;
	alpha = 0;
}

void CPML::InitializeCPML(Structure &sS)
{
	cout << BOLDYELLOW << " \u2514" << RESET << "Initializing " << BOLDBLUE << "CPML " << RESET << "coefficients...";
	int i, j, k;
	QHxz1 = new double**[Npmlx1 + 1];
	QHxz2 = new double**[Npmlx2 + 1];
	QExy1 = new double**[Npmlx1 + 1];
	QExy2 = new double**[Npmlx2 + 1];
	QHxy1 = new double**[Npmlx1 + 1];
	QHxy2 = new double**[Npmlx2 + 1];
	QExz1 = new double**[Npmlx1 + 1];
	QExz2 = new double**[Npmlx2 + 1];
	for (i = 0; i <= Npmlx1; i++){
		QHxz1[i] = new double*[Ny + 1];
		QHxz2[i] = new double*[Ny + 1];
		QExy1[i] = new double*[Ny + 1];
		QExy2[i] = new double*[Ny + 1];
		QHxy1[i] = new double*[Ny + 1];
		QHxy2[i] = new double*[Ny + 1];
		QExz1[i] = new double*[Ny + 1];
		QExz2[i] = new double*[Ny + 1];
		for (j = 0; j <= Ny; j++){
			QHxz1[i][j] = new double[Nz + 1];
			QHxz2[i][j] = new double[Nz + 1];
			QExy1[i][j] = new double[Nz + 1];
			QExy2[i][j] = new double[Nz + 1];
			QHxy1[i][j] = new double[Nz + 1];
			QHxy2[i][j] = new double[Nz + 1];
			QExz1[i][j] = new double[Nz + 1];
			QExz2[i][j] = new double[Nz + 1];
			for (k = 0; k <= Nz; k++){
				QHxz1[i][j][k] = 0.0;
				QHxz2[i][j][k] = 0.0;
				QExy1[i][j][k] = 0.0;
				QExy2[i][j][k] = 0.0;
				QHxy1[i][j][k] = 0.0;
				QHxy2[i][j][k] = 0.0;
				QExz1[i][j][k] = 0.0;
				QExz2[i][j][k] = 0.0;
			}
		}
	}

	QHyz1 = new double**[Nx + 1];
	QHyz2 = new double**[Nx + 1];
	QEyz1 = new double**[Nx + 1];
	QEyz2 = new double**[Nx + 1];
	QHyx1 = new double**[Nx + 1];
	QHyx2 = new double**[Nx + 1];
	QEyx1 = new double**[Nx + 1];
	QEyx2 = new double**[Nx + 1];
	for (i = 0; i <= Nx; i++){
		QHyz1[i] = new double*[Npmly1 + 1];
		QHyz2[i] = new double*[Npmly2 + 1];
		QEyz1[i] = new double*[Npmly1 + 1];
		QEyz2[i] = new double*[Npmly2 + 1];
		QHyx1[i] = new double*[Npmly1 + 1];
		QHyx2[i] = new double*[Npmly2 + 1];
		QEyx1[i] = new double*[Npmly1 + 1];
		QEyx2[i] = new double*[Npmly2 + 1];
		for (j = 0; j <= Npmly1; j++){
			QHyz1[i][j] = new double[Nz + 1];
			QHyz2[i][j] = new double[Nz + 1];
			QEyz1[i][j] = new double[Nz + 1];
			QEyz2[i][j] = new double[Nz + 1];
			QHyx1[i][j] = new double[Nz + 1];
			QHyx2[i][j] = new double[Nz + 1];
			QEyx1[i][j] = new double[Nz + 1];
			QEyx2[i][j] = new double[Nz + 1];
			for (k = 0; k <= Nz; k++){
				QHyz1[i][j][k] = 0.0;
				QHyz2[i][j][k] = 0.0;
				QEyz1[i][j][k] = 0.0;
				QEyz2[i][j][k] = 0.0;
				QHyx1[i][j][k] = 0.0;
				QHyx2[i][j][k] = 0.0;
				QEyx1[i][j][k] = 0.0;
				QEyx2[i][j][k] = 0.0;
			}
		}
	}

	QHzy1 = new double**[Nx + 1];
	QHzy2 = new double**[Nx + 1];
	QEzy1 = new double**[Nx + 1];
	QEzy2 = new double**[Nx + 1];
	QHzx1 = new double**[Nx + 1];
	QHzx2 = new double**[Nx + 1];
	QEzx1 = new double**[Nx + 1];
	QEzx2 = new double**[Nx + 1];
	for (i = 0; i <= Nx; i++){
		QHzy1[i] = new double*[Ny + 1];
		QHzy2[i] = new double*[Ny + 1];
		QEzy1[i] = new double*[Ny + 1];
		QEzy2[i] = new double*[Ny + 1];
		QHzx1[i] = new double*[Ny + 1];
		QHzx2[i] = new double*[Ny + 1];
		QEzx1[i] = new double*[Ny + 1];
		QEzx2[i] = new double*[Ny + 1];
		for (j = 0; j <= Ny; j++){
			QHzy1[i][j] = new double[Npmlz1 + 1];
			QHzy2[i][j] = new double[Npmlz2 + 1];
			QEzy1[i][j] = new double[Npmlz1 + 1];
			QEzy2[i][j] = new double[Npmlz2 + 1];
			QHzx1[i][j] = new double[Npmlz1 + 1];
			QHzx2[i][j] = new double[Npmlz2 + 1];
			QEzx1[i][j] = new double[Npmlz1 + 1];
			QEzx2[i][j] = new double[Npmlz2 + 1];
			for (k = 0; k <= Npmlz1; k++){
				QHzy1[i][j][k] = 0.0;
				QHzy2[i][j][k] = 0.0;
				QEzy1[i][j][k] = 0.0;
				QEzy2[i][j][k] = 0.0;
				QHzx1[i][j][k] = 0.0;
				QHzx2[i][j][k] = 0.0;
				QEzx1[i][j][k] = 0.0;
				QEzx2[i][j][k] = 0.0;
			}
		}
	}


	double eta0 = sqrt(mu0/eps0);
//E
//x-directed
	sigmax = sig_ratio*(0.8 * (m + 1 ) / (sS.dx * eta0));
	kappa = new double[Npmlx1+1];
	bxE1 = new double[Npmlx1+1];
	cxE1 = new double[Npmlx1+1];
	E_pml_Coef(kappa, bxE1, cxE1, true, Npmlx1, sS.dt, sS.dx, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i <= Npmlx1; ++i){
		sS.cEzy[i] = sS.cEzy[i]/kappa[i];
		sS.cEyz[i] = sS.cEyz[i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
	//cout << endl;
	delete[] kappa;
	kappa = 0;

	kappa = new double[Npmlx2+1];
	bxE2 = new double[Npmlx2+1];
	cxE2 = new double[Npmlx2+1];
	E_pml_Coef(kappa, bxE2, cxE2, false, Npmlx2, sS.dt, sS.dx, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmlx1; ++i){
		sS.cEzy[sS.Nx-Npmlx2+i] = sS.cEzy[sS.Nx-Npmlx2+i]/kappa[i];
		sS.cEyz[sS.Nx-Npmlx2+i] = sS.cEyz[sS.Nx-Npmlx2+i]/kappa[i];
	//cout << cxE2[i] << "\t";
	}
	//cout << endl;
	delete[] kappa;
	kappa = 0;

//y-directed
	sigmax = sig_ratio*(0.8 * (m + 1 ) / (sS.dy * eta0));
	kappa = new double[Npmly1+1];
	byE1 = new double[Npmly1+1];
	cyE1 = new double[Npmly1+1];
	E_pml_Coef(kappa, byE1, cyE1, true, Npmly1, sS.dt, sS.dy, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i <= Npmly1; ++i){
		sS.cEzx[i] = sS.cEzx[i]/kappa[i];
		sS.cExz[i] = sS.cExz[i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

	kappa = new double[Npmly2+1];
	byE2 = new double[Npmly2+1];
	cyE2 = new double[Npmly2+1];
	E_pml_Coef(kappa, byE2, cyE2, false, Npmly2, sS.dt, sS.dy, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmly1; ++i){
		sS.cEzx[sS.Ny-Npmly2+i] = sS.cEzx[sS.Ny-Npmly2+i]/kappa[i];
		sS.cExz[sS.Ny-Npmly2+i] = sS.cExz[sS.Ny-Npmly2+i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

//z-directed
	sigmax = sig_ratio*(0.8 * (m + 1 ) / (sS.dz * eta0));
	kappa = new double[Npmlz1+1];
	bzE1 = new double[Npmlz1+1];
	czE1 = new double[Npmlz1+1];
	E_pml_Coef(kappa, bzE1, czE1, true, Npmlz1, sS.dt, sS.dz, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i <= Npmlz1; ++i){
		sS.cExy[i] = sS.cExy[i]/kappa[i];
		sS.cEyx[i] = sS.cEyx[i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

	kappa = new double[Npmlz2+1];
	bzE2 = new double[Npmlz2+1];
	czE2 = new double[Npmlz2+1];
	E_pml_Coef(kappa, bzE2, czE2, false, Npmlz2, sS.dt, sS.dz, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmlz1; ++i){
		sS.cExy[sS.Nz-Npmlz2+i] = sS.cExy[sS.Nz-Npmlz2+i]/kappa[i];
		sS.cEyx[sS.Nz-Npmlz2+i] = sS.cEyx[sS.Nz-Npmlz2+i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

//H
//x-directed
	sigmax = sig_ratio*(0.8 * (m + 1 ) / (sS.dx * eta0));
	kappa = new double[Npmlx1];
	bxH1 = new double[Npmlx1];
	cxH1 = new double[Npmlx1];
	H_pml_Coef(kappa, bxH1, cxH1, true, Npmlx1, sS.dt, sS.dx, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmlx1; ++i){
		sS.cHzy[i] = sS.cHzy[i]/kappa[i];
		sS.cHyz[i] = sS.cHyz[i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

	kappa = new double[Npmlx2];
	bxH2 = new double[Npmlx2];
	cxH2 = new double[Npmlx2];
	H_pml_Coef(kappa, bxH2, cxH2, false, Npmlx2, sS.dt, sS.dx, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmlx1; ++i){
		sS.cHzy[sS.Nx-Npmlx2+i] = sS.cHzy[sS.Nx-Npmlx2+i]/kappa[i];
		sS.cHyz[sS.Nx-Npmlx2+i] = sS.cHyz[sS.Nx-Npmlx2+i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

//y-directed
	sigmax = sig_ratio*(0.8 * (m + 1 ) / (sS.dy * eta0));
	kappa = new double[Npmly1];
	byH1 = new double[Npmly1];
	cyH1 = new double[Npmly1];
	H_pml_Coef(kappa, byH1, cyH1, true, Npmly1, sS.dt, sS.dy, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmly1; ++i){
		sS.cHzx[i] = sS.cHzx[i]/kappa[i];
		sS.cHxz[i] = sS.cHxz[i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

	kappa = new double[Npmly2];
	byH2 = new double[Npmly2];
	cyH2 = new double[Npmly2];
	H_pml_Coef(kappa, byH2, cyH2, false, Npmly2, sS.dt, sS.dy, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmly1; ++i){
		sS.cHzx[sS.Ny-Npmly2+i] = sS.cHzx[sS.Ny-Npmly2+i]/kappa[i];
		sS.cHxz[sS.Ny-Npmly2+i] = sS.cHxz[sS.Ny-Npmly2+i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

//z-directed
	sigmax = sig_ratio*(0.8 * (m + 1 ) / (sS.dz * eta0));
	kappa = new double[Npmlz1];
	bzH1 = new double[Npmlz1];
	czH1 = new double[Npmlz1];
	H_pml_Coef(kappa, bzH1, czH1, true, Npmlz1, sS.dt, sS.dz, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmlz1; ++i){
		sS.cHxy[i] = sS.cHxy[i]/kappa[i];
		sS.cHyx[i] = sS.cHyx[i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

	kappa = new double[Npmlz2];
	bzH2 = new double[Npmlz2];
	czH2 = new double[Npmlz2];
	H_pml_Coef(kappa, bzH2, czH2, false, Npmlz2, sS.dt, sS.dz, m, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmlz1; ++i){
		sS.cHxy[sS.Nz-Npmlz2+i] = sS.cHxy[sS.Nz-Npmlz2+i]/kappa[i];
		sS.cHyx[sS.Nz-Npmlz2+i] = sS.cHyx[sS.Nz-Npmlz2+i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

	cpmlE = sS.dt/eps0;//sS.bE[0];
	cpmlH = sS.bH[0];//sS.dt/mu0;
	cout << " Done." << endl;
}

void CPML::UpdateCPMLEFields(Field &fF)
{
	int i, j, k, ii, jj, kk;
// x-pml
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 1; i <= Npmlx1; ++i){
		for (j = 1; j < Ny; ++j){
			for (k = 1; k < Nz; ++k){
				QHxz1[i][j][k] = bxE1[i] * QHxz1[i][j][k] + cxE1[i] * (fF.Hz[i][j][k] - fF.Hz[i - 1][j][k]);
				fF.Ey[i][j][k] = fF.Ey[i][j][k] - QHxz1[i][j][k] * cpmlE;
				QHxy1[i][j][k] = bxE1[i] * QHxy1[i][j][k] + cxE1[i] * (fF.Hy[i][j][k] - fF.Hy[i - 1][j][k]);
				fF.Ez[i][j][k] = fF.Ez[i][j][k] + QHxy1[i][j][k] * cpmlE;
			}
		}
	}
	i = Nx-Npmlx2;//1;
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (ii = 0; ii < Npmlx2; ++ii){
		for (j = 1; j < Ny; ++j){
			for (k = 1; k < Nz; ++k){
				QHxz2[ii][j][k] = bxE2[ii] * QHxz2[ii][j][k] + cxE2[ii] * (fF.Hz[i][j][k] - fF.Hz[i - 1][j][k]);
				fF.Ey[i][j][k] = fF.Ey[i][j][k] - QHxz2[ii][j][k] * cpmlE;
				QHxy2[ii][j][k] = bxE2[ii] * QHxy2[ii][j][k] + cxE2[ii] * (fF.Hy[i][j][k] - fF.Hy[i - 1][j][k]);
				fF.Ez[i][j][k] = fF.Ez[i][j][k] + QHxy2[ii][j][k] * cpmlE;
			}
		}
		i = i + 1;
	}
// y-pml
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 1; i < Nx; ++i){
		for (j = 1; j <= (Npmly1); ++j){
			for (k = 1; k < Nz; ++k){
				QHyx1[i][j][k] = byE1[j] * QHyx1[i][j][k] + cyE1[j] * (fF.Hx[i][j][k] - fF.Hx[i][j - 1][k]);
				fF.Ez[i][j][k] = fF.Ez[i][j][k] - QHyx1[i][j][k] * cpmlE;
				QHyz1[i][j][k] = byE1[j] * QHyz1[i][j][k] + cyE1[j] * (fF.Hz[i][j][k] - fF.Hz[i][j - 1][k]);
				fF.Ex[i][j][k] = fF.Ex[i][j][k] + QHyz1[i][j][k] * cpmlE;
			}
		}
	}
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 1; i< Nx; ++i){
		j = Ny-Npmly2;//1;
		for (jj = 0; jj < Npmly2; ++jj){
			for (k = 1; k < Nz; ++k){
				QHyx2[i][jj][k] = byE2[jj] * QHyx2[i][jj][k] + cyE2[jj] * (fF.Hx[i][j][k] - fF.Hx[i][j - 1][k]);
				fF.Ez[i][j][k] = fF.Ez[i][j][k] - QHyx2[i][jj][k] * cpmlE;
				QHyz2[i][jj][k] = byE2[jj] * QHyz2[i][jj][k] + cyE2[jj] * (fF.Hz[i][j][k] - fF.Hz[i][j - 1][k]);
				fF.Ex[i][j][k] = fF.Ex[i][j][k] + QHyz2[i][jj][k] * cpmlE;
			}
			j = j + 1;
		}
	}
// z-pml
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 1; i < Nx; ++i){
		for (j = 1; j < Ny ; ++j){
			for (k = 1; k <= (Npmlz1); ++k){
				QHzx1[i][j][k] = bzE1[k] * QHzx1[i][j][k] + czE1[k] * (fF.Hy[i][j][k] - fF.Hy[i][j][k - 1]);
				fF.Ex[i][j][k] = fF.Ex[i][j][k] - QHzx1[i][j][k] * cpmlE;
				QHzy1[i][j][k] = bzE1[k] * QHzy1[i][j][k] + czE1[k] * (fF.Hx[i][j][k] - fF.Hx[i][j][k - 1]);
				fF.Ey[i][j][k] = fF.Ey[i][j][k] + QHzy1[i][j][k] * cpmlE;
			}
		}
	}
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 1; i < Nx; ++i){
		for (j = 1; j < Ny ; ++j){
			k = Nz-Npmlz2;//1;
			for (kk = 0; kk < Npmlz2; ++kk){
				QHzx2[i][j][kk] = bzE2[kk] * QHzx2[i][j][kk] + czE2[kk] * (fF.Hy[i][j][k] - fF.Hy[i][j][k - 1]);
				fF.Ex[i][j][k] = fF.Ex[i][j][k] - QHzx2[i][j][kk] * cpmlE;
				QHzy2[i][j][kk] = bzE2[kk] * QHzy2[i][j][kk] + czE2[kk] * (fF.Hx[i][j][k] - fF.Hx[i][j][k - 1]);
				fF.Ey[i][j][k] = fF.Ey[i][j][k] + QHzy2[i][j][kk] * cpmlE;
				k = k + 1;
			}
		}
	}
}

void CPML::UpdateCPMLHFields(Field &fF)
{
	int i, j, k, ii, jj, kk;
//x-pml
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 0; i < Npmlx1 ; ++i){
		for (j = 0; j < Ny; ++j){
			for (k = 0; k < Nz; ++k){
				QExz1[i][j][k] = bxH1[i] * QExz1[i][j][k] + cxH1[i] * (fF.Ez[i + 1][j][k] - fF.Ez[i][j][k]);
				fF.Hy[i][j][k] = fF.Hy[i][j][k] + QExz1[i][j][k] * cpmlH;
				QExy1[i][j][k] = bxH1[i] * QExy1[i][j][k] + cxH1[i] * (fF.Ey[i + 1][j][k] - fF.Ey[i][j][k]);
				fF.Hz[i][j][k] = fF.Hz[i][j][k] - QExy1[i][j][k] * cpmlH;
			}
		}
	}
	i = Nx-Npmlx2;//1;
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (ii = 0; ii < Npmlx2; ++ii){
		for (j = 0; j < Ny; ++j){
			for (k = 0; k < Nz; ++k){//
				QExz2[ii][j][k] = bxH2[ii] * QExz2[ii][j][k] + cxH2[ii] * (fF.Ez[i + 1][j][k] - fF.Ez[i][j][k]);
				fF.Hy[i][j][k] = fF.Hy[i][j][k] + QExz2[ii][j][k] * cpmlH;
				QExy2[ii][j][k] = bxH2[ii] * QExy2[ii][j][k] + cxH2[ii] * (fF.Ey[i + 1][j][k] - fF.Ey[i][j][k]);
				fF.Hz[i][j][k] = fF.Hz[i][j][k] - QExy2[ii][j][k] * cpmlH;
			}
		}
		i = i + 1;
	}
//y-pml
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 0; i < Nx ; ++i){
		for (j = 0; j < Npmly1; ++j){
			for (k = 0; k < Nz; ++k){
				QEyx1[i][j][k] = byH1[j] * QEyx1[i][j][k] + cyH1[j] * (fF.Ex[i][j + 1][k] - fF.Ex[i][j][k]);
				fF.Hz[i][j][k] = fF.Hz[i][j][k] - QEyx1[i][j][k] * cpmlH;
				QEyz1[i][j][k] = byH1[j] * QEyz1[i][j][k] + cyH1[j] * (fF.Ez[i][j + 1][k] - fF.Ez[i][j][k]);
				fF.Hx[i][j][k] = fF.Hx[i][j][k] + QEyz1[i][j][k] * cpmlH;
			}
		}
	}
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 0; i < Nx; ++i){
		j = Ny-Npmly2;//1;
		for (jj = 0; jj < Npmly2; ++jj){
			for (k = 0; k < Nz; ++k){
				QEyx2[i][jj][k] = byH2[jj] * QEyx2[i][jj][k] + cyH2[jj] * (fF.Ex[i][j + 1][k] - fF.Ex[i][j][k]);
				fF.Hz[i][j][k] = fF.Hz[i][j][k] - QEyx2[i][jj][k] * cpmlH;
				QEyz2[i][jj][k] = byH2[jj] * QEyz2[i][jj][k] + cyH2[jj] * (fF.Ez[i][j + 1][k] - fF.Ez[i][j][k]);
				fF.Hx[i][j][k] = fF.Hx[i][j][k] + QEyz2[i][jj][k] * cpmlH;
			}
			j = j + 1;
		}
	}
//z-pml
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 0; i < Nx ; ++i){
		for (j = 0; j < Ny; ++j){
			for (k = 0; k < Npmlz1; ++k){
				QEzy1[i][j][k] = bzH1[k] * QEzy1[i][j][k] + czH1[k] * (fF.Ey[i][j][k + 1] - fF.Ey[i][j][k]);
				fF.Hx[i][j][k] = fF.Hx[i][j][k] - QEzy1[i][j][k] * cpmlH;
				QEzx1[i][j][k] = bzH1[k] * QEzx1[i][j][k] + czH1[k] * (fF.Ex[i][j][k + 1] - fF.Ex[i][j][k]);
				fF.Hy[i][j][k] = fF.Hy[i][j][k] + QEzx1[i][j][k] * cpmlH;
			}
		}
	}
//#pragma omp parallel for default(shared) private(i,j,k,ii,jj,kk)
	for (i = 0; i < Nx ; ++i){
		for (j = 0; j < Ny; ++j){
			k = Nz-Npmlz2;//1;
			for (kk = 0; kk < Npmlz2; ++kk){
				QEzy2[i][j][kk] = bzH2[kk] * QEzy2[i][j][kk] + czH2[kk] * (fF.Ey[i][j][k + 1] - fF.Ey[i][j][k]);
				fF.Hx[i][j][k] = fF.Hx[i][j][k] - QEzy2[i][j][kk] * cpmlH;
				QEzx2[i][j][kk] = bzH2[kk] * QEzx2[i][j][kk] + czH2[kk] * (fF.Ex[i][j][k + 1] - fF.Ex[i][j][k]);
				fF.Hy[i][j][k] = fF.Hy[i][j][k] + QEzx2[i][j][kk] * cpmlH;
				k = k + 1;
			}
		}
	}


}
