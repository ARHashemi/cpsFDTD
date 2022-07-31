#include "CFS_PML.h"
#include "stdafx.h"
#include "Constants.h"
#include "Structure.h"
#include "Fields.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <omp.h>

using namespace std;

CFS_PML::CFS_PML(Structure &sS, char *fileName)
{
	cout << BOLDRED << ">" << RESET << "Constructing " << BOLDBLUE << "CPML " << RESET << "boundaries..." << endl;
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
	iss >> p_order;
	iss >> malp;
	//cout << p_order << "\t" << malp << "\t";
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> kapmax;
	iss >> sig_ratio;
	//cout << kapmax << "\t" << sig_ratio << "\t";
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> alpmax;
	//cout << alpmax << "\t";
	//cout << m << "\t" << kapmax << endl;
	Nx = sS.Nx;
	Ny = sS.Ny;
	Nz = sS.Nz;
	dx = sS.dx;
	dy = sS.dy;
	dz = sS.dz;

	CpmlE = sS.bE[0];
	CpmlH = sS.bH[0];

	//Num_of_PML_points = Npmlx1*(Ny*Nz) + Npmlx2*(Ny*Nz) + Npmly1*(Nx*Nz) + Npmly2*(Nx*Nz) + Npmlz1*(Nx*Ny) + Npmlz2*(Nx*Ny);

	//cout << " Done." << endl;
}

CFS_PML::~CFS_PML()
{
	//dtor
}

void CFS_PML::E_pml_Coef(double* &kappa, double* &bcoef, double* &ccoef, bool firstend, int npml, double dt, double del, double m, double ma, double sigmax, double kapmax, double alpmax)
{
	double polyscalesk, polyscalea, *sigma, *alpha, dtotau, temp;
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
			ccoef[i] = sigma[i] * (bcoef[i] - 1.0) / (sigma[i] + kappa[i] * alpha[i]) / kappa[i];


		//cout << ccoef[i] << "\t";
	}
	delete[] sigma;
	delete[] alpha;
	sigma = 0;
	alpha = 0;
}

void CFS_PML::H_pml_Coef(double* &kappa, double* &bcoef, double* &ccoef, bool firstend, int npml, double dt, double del, double m, double ma, double sigmax, double kapmax, double alpmax)
{
	double polyscalesk, polyscalea, *sigma, *alpha, dtotau, temp;
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
			ccoef[i] = sigma[i] * (bcoef[i] - 1.0) / (sigma[i] + kappa[i] * alpha[i]) / kappa[i];

	}
	delete[] sigma;
	delete[] alpha;
	sigma = 0;
	alpha = 0;
}

void CFS_PML::Initialize_PML(Structure &sS)
{
	cout << BOLDYELLOW << " \u2514" << RESET << "Initializing " << BOLDBLUE << "CPML " << RESET << "coefficients...";
	double eta0 = sqrt(mu0/eps0);
//E
//x-directed
	sigmax = sig_ratio*(0.8 * (p_order + 1.0 ) / (sS.dx * eta0));
	cout << sigmax << "\t";
	kappa = new double[Npmlx1+1];
	bxE1 = new double[Npmlx1+1];
	cxE1 = new double[Npmlx1+1];
	E_pml_Coef(kappa, bxE1, cxE1, true, Npmlx1, sS.dt, sS.dx, p_order, malp, sigmax, kapmax, alpmax);
	for (i = 0; i <= Npmlx1; ++i){
		sS.cEzy[i] = sS.cEzy[i]/kappa[i];
		sS.cEyz[i] = sS.cEyz[i]/kappa[i];
		//cout << cxE1[i] << "\t";
	}
	//cout << endl;
	delete[] kappa;
	kappa = 0;

	kappa = new double[Npmlx2+1];
	bxE2 = new double[Npmlx2+1];
	cxE2 = new double[Npmlx2+1];
	E_pml_Coef(kappa, bxE2, cxE2, false, Npmlx2, sS.dt, sS.dx, p_order, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmlx1; ++i){
		sS.cEzy[sS.Nx-Npmlx2+i] = sS.cEzy[sS.Nx-Npmlx2+i]/kappa[i];
		sS.cEyz[sS.Nx-Npmlx2+i] = sS.cEyz[sS.Nx-Npmlx2+i]/kappa[i];
	//cout << cxE2[i] << "\t";
	}
	//cout << endl;
	delete[] kappa;
	kappa = 0;

//y-directed
	sigmax = sig_ratio*(0.8 * (p_order + 1 ) / (sS.dy * eta0));
	kappa = new double[Npmly1+1];
	byE1 = new double[Npmly1+1];
	cyE1 = new double[Npmly1+1];
	E_pml_Coef(kappa, byE1, cyE1, true, Npmly1, sS.dt, sS.dy, p_order, malp, sigmax, kapmax, alpmax);
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
	E_pml_Coef(kappa, byE2, cyE2, false, Npmly2, sS.dt, sS.dy, p_order, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmly1; ++i){
		sS.cEzx[sS.Ny-Npmly2+i] = sS.cEzx[sS.Ny-Npmly2+i]/kappa[i];
		sS.cExz[sS.Ny-Npmly2+i] = sS.cExz[sS.Ny-Npmly2+i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

//z-directed
	sigmax = sig_ratio*(0.8 * (p_order + 1 ) / (sS.dz * eta0));
	kappa = new double[Npmlz1+1];
	bzE1 = new double[Npmlz1+1];
	czE1 = new double[Npmlz1+1];
	E_pml_Coef(kappa, bzE1, czE1, true, Npmlz1, sS.dt, sS.dz, p_order, malp, sigmax, kapmax, alpmax);
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
	E_pml_Coef(kappa, bzE2, czE2, false, Npmlz2, sS.dt, sS.dz, p_order, malp, sigmax, kapmax, alpmax);
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
	sigmax = sig_ratio*(0.8 * (p_order + 1 ) / (sS.dx * eta0));
	kappa = new double[Npmlx1];
	bxH1 = new double[Npmlx1];
	cxH1 = new double[Npmlx1];
	H_pml_Coef(kappa, bxH1, cxH1, true, Npmlx1, sS.dt, sS.dx, p_order, malp, sigmax, kapmax, alpmax);
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
	H_pml_Coef(kappa, bxH2, cxH2, false, Npmlx2, sS.dt, sS.dx, p_order, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmlx1; ++i){
		sS.cHzy[sS.Nx-Npmlx2+i] = sS.cHzy[sS.Nx-Npmlx2+i]/kappa[i];
		sS.cHyz[sS.Nx-Npmlx2+i] = sS.cHyz[sS.Nx-Npmlx2+i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

//y-directed
	sigmax = sig_ratio*(0.8 * (p_order + 1 ) / (sS.dy * eta0));
	kappa = new double[Npmly1];
	byH1 = new double[Npmly1];
	cyH1 = new double[Npmly1];
	H_pml_Coef(kappa, byH1, cyH1, true, Npmly1, sS.dt, sS.dy, p_order, malp, sigmax, kapmax, alpmax);
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
	H_pml_Coef(kappa, byH2, cyH2, false, Npmly2, sS.dt, sS.dy, p_order, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmly1; ++i){
		sS.cHzx[sS.Ny-Npmly2+i] = sS.cHzx[sS.Ny-Npmly2+i]/kappa[i];
		sS.cHxz[sS.Ny-Npmly2+i] = sS.cHxz[sS.Ny-Npmly2+i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;

//z-directed
	sigmax = sig_ratio*(0.8 * (p_order + 1 ) / (sS.dz * eta0));
	kappa = new double[Npmlz1];
	bzH1 = new double[Npmlz1];
	czH1 = new double[Npmlz1];
	H_pml_Coef(kappa, bzH1, czH1, true, Npmlz1, sS.dt, sS.dz, p_order, malp, sigmax, kapmax, alpmax);
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
	H_pml_Coef(kappa, bzH2, czH2, false, Npmlz2, sS.dt, sS.dz, p_order, malp, sigmax, kapmax, alpmax);
	for (i = 0; i < Npmlz1; ++i){
		sS.cHxy[sS.Nz-Npmlz2+i] = sS.cHxy[sS.Nz-Npmlz2+i]/kappa[i];
		sS.cHyx[sS.Nz-Npmlz2+i] = sS.cHyx[sS.Nz-Npmlz2+i]/kappa[i];
	//cout << kappa[i] << "\t";
	}
//	cout << endl;
	delete[] kappa;
	kappa = 0;


    Num_of_PML_points = (Nz-1)*(Ny-1)*(Npmlx1) + (Nz-1)*(Ny-1)*(Npmlx2) + (Nz-1)*(Nx-1)*(Npmly1) + (Nz-1)*(Nx-1)*(Npmly2) + (Nx-1)*(Ny-1)*(Npmlz1) + (Nx-1)*(Ny-1)*(Npmlz2);
    PMLpointE = new PML_point[Num_of_PML_points];
    int ii = 0, jj;
    for (k = 1; k < sS.Nz; ++k){
		for (j = 1; j < sS.Ny; ++j){
			for (i = 1; i <= Npmlx1; ++i){
                PMLpointE[ii].i = i;
                PMLpointE[ii].j = j;
                PMLpointE[ii].k = k;
                PMLpointE[ii].Q1 = 0.0;
                PMLpointE[ii].Q2 = 0.0;
                PMLpointE[ii].grade = i;
                ++ii;
			}
		}
    }
    iixE1 = ii;
    for (k = 1; k < sS.Nz; ++k){
		for (j = 1; j < sS.Ny; ++j){
			jj = Nx - Npmlx2;
			for (i = 0; i < Npmlx2; ++i){
                PMLpointE[ii].i = jj;
                PMLpointE[ii].j = j;
                PMLpointE[ii].k = k;
                PMLpointE[ii].Q1 = 0.0;
                PMLpointE[ii].Q2 = 0.0;
                PMLpointE[ii].grade = i;
                ++ii;
                ++jj;
			}
		}
    }
    iixE2 = ii;

    for (k = 1; k < sS.Nz; ++k){
		for (j = 1; j <= Npmly1; ++j){
			for (i = 1; i < sS.Nx; ++i){
                PMLpointE[ii].i = i;
                PMLpointE[ii].j = j;
                PMLpointE[ii].k = k;
                PMLpointE[ii].Q1 = 0.0;
                PMLpointE[ii].Q2 = 0.0;
                PMLpointE[ii].grade = j;
                ++ii;
			}
		}
    }
    iiyE1 = ii;
    for (k = 1; k < sS.Nz; ++k){
		jj = Ny - Npmly2;
		for (j = 0; j < Npmly2; ++j){
			for (i = 1; i < sS.Nx; ++i){
                PMLpointE[ii].i = i;
                PMLpointE[ii].j = jj;
                PMLpointE[ii].k = k;
                PMLpointE[ii].Q1 = 0.0;
                PMLpointE[ii].Q2 = 0.0;
                PMLpointE[ii].grade = j;
                ++ii;
			}
			++jj;
		}
    }
    iiyE2 = ii;

    for (k = 1; k <= Npmly1; ++k){
		for (j = 1; j < sS.Ny; ++j){
			for (i = 1; i < sS.Nx; ++i){
                PMLpointE[ii].i = i;
                PMLpointE[ii].j = j;
                PMLpointE[ii].k = k;
                PMLpointE[ii].Q1 = 0.0;
                PMLpointE[ii].Q2 = 0.0;
                PMLpointE[ii].grade = k;
                ++ii;
			}
		}
    }
    iizE1 = ii;
    jj = Nz - Npmlz2;
    for (k = 0; k < Npmlz2; ++k){
		for (j = 1; j < sS.Ny; ++j){
			for (i = 1; i < sS.Nx; ++i){
                PMLpointE[ii].i = i;
                PMLpointE[ii].j = j;
                PMLpointE[ii].k = jj;
                PMLpointE[ii].Q1 = 0.0;
                PMLpointE[ii].Q2 = 0.0;
                PMLpointE[ii].grade = k;
                ++ii;
			}
		}
		++jj;
    }
    iizE2 = ii;

    //H
    PMLpointH = new PML_point[Num_of_PML_points];
    ii = 0;
    for (k = 1; k < sS.Nz; ++k){
		for (j = 1; j < sS.Ny; ++j){
			for (i = 0; i < Npmlx1; ++i){
                PMLpointH[ii].i = i;
                PMLpointH[ii].j = j;
                PMLpointH[ii].k = k;
                PMLpointH[ii].Q1 = 0.0;
                PMLpointH[ii].Q2 = 0.0;
                PMLpointH[ii].grade = i;
                ++ii;
			}
		}
    }
    iixH1 = ii;
    for (k = 1; k < sS.Nz; ++k){
		for (j = 1; j < sS.Ny; ++j){
			jj = Nx - Npmlx2;
			for (i = 0; i < Npmlx2; ++i){
                PMLpointH[ii].i = jj;
                PMLpointH[ii].j = j;
                PMLpointH[ii].k = k;
                PMLpointH[ii].Q1 = 0.0;
                PMLpointH[ii].Q2 = 0.0;
                PMLpointH[ii].grade = i;
                ++ii;
                ++jj;
			}
		}
    }
    iixH2 = ii;

     for (k = 1; k < sS.Nz; ++k){
		for (j = 0; j < Npmly1; ++j){
			for (i = 1; i < sS.Nx; ++i){
                PMLpointH[ii].i = i;
                PMLpointH[ii].j = j;
                PMLpointH[ii].k = k;
                PMLpointH[ii].Q1 = 0.0;
                PMLpointH[ii].Q2 = 0.0;
                PMLpointH[ii].grade = j;
                ++ii;
			}
		}
    }
    iiyH1 = ii;
    for (k = 1; k < sS.Nz; ++k){
		jj = Ny - Npmly2;
		for (j = 0; j < Npmly2; ++j){
			for (i = 1; i < sS.Nx; ++i){
                PMLpointH[ii].i = i;
                PMLpointH[ii].j = jj;
                PMLpointH[ii].k = k;
                PMLpointH[ii].Q1 = 0.0;
                PMLpointH[ii].Q2 = 0.0;
                PMLpointH[ii].grade = j;
                ++ii;
			}
			++jj;
		}
    }
    iiyH2 = ii;

    for (k = 0; k < Npmly1; ++k){
		for (j = 1; j < sS.Ny; ++j){
			for (i = 1; i < sS.Nx; ++i){
                PMLpointH[ii].i = i;
                PMLpointH[ii].j = j;
                PMLpointH[ii].k = k;
                PMLpointH[ii].Q1 = 0.0;
                PMLpointH[ii].Q2 = 0.0;
                PMLpointH[ii].grade = k;
                ++ii;
			}
		}
    }
    iizH1 = ii;
    jj = Nz - Npmlz2;
    for (k = 0; k < Npmlz2; ++k){
		for (j = 1; j < sS.Ny; ++j){
			for (i = 1; i < sS.Nx; ++i){
                PMLpointH[ii].i = i;
                PMLpointH[ii].j = j;
                PMLpointH[ii].k = jj;
                PMLpointH[ii].Q1 = 0.0;
                PMLpointH[ii].Q2 = 0.0;
                PMLpointH[ii].grade = k;
                ++ii;
			}
		}
		++jj;
    }
    iizH2 = ii;

	cout << " Done." << endl;
}

void CFS_PML::UpdatePMLEFields(Field &fF)
{
	int ii;
//#pragma omp parallel for default(shared) private(ii)
//	for (ii = 0; ii < iixE1; ++ii){
//		PMLpointE[ii].Q1 = bxE1[PMLpointE[ii].grade] * PMLpointE[ii].Q1 + cxE1[PMLpointE[ii].grade] * (fF.Hz[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hz[PMLpointE[ii].i-1][PMLpointE[ii].j][PMLpointE[ii].k]);
//		fF.Ey[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ey[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - CpmlE * PMLpointE[ii].Q1;
//
//		PMLpointE[ii].Q2 = bxE1[PMLpointE[ii].grade] * PMLpointE[ii].Q2 + cxE1[PMLpointE[ii].grade] * (fF.Hy[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hy[PMLpointE[ii].i-1][PMLpointE[ii].j][PMLpointE[ii].k]);
//		fF.Ez[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ez[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] + CpmlE * PMLpointE[ii].Q2;
//	}
#pragma omp parallel for default(shared) private(ii)
	for (ii = iixE1; ii < iixE2; ++ii){
		PMLpointE[ii].Q1 = bxE2[PMLpointE[ii].grade] * PMLpointE[ii].Q1 + cxE2[PMLpointE[ii].grade] * (fF.Hz[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hz[PMLpointE[ii].i-1][PMLpointE[ii].j][PMLpointE[ii].k]);
		fF.Ey[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ey[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - CpmlE * PMLpointE[ii].Q1;

		PMLpointE[ii].Q2 = bxE2[PMLpointE[ii].grade] * PMLpointE[ii].Q2 + cxE2[PMLpointE[ii].grade] * (fF.Hy[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hy[PMLpointE[ii].i-1][PMLpointE[ii].j][PMLpointE[ii].k]);
		fF.Ez[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ez[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] + CpmlE * PMLpointE[ii].Q2;
	}
	cout << PMLpointE[ii-1].Q1 << "\t";

#pragma omp parallel for default(shared) private(ii)
	for (ii = iixE2; ii < iiyE1; ++ii){
		PMLpointE[ii].Q1 = byE1[PMLpointE[ii].grade] * PMLpointE[ii].Q1 + cyE1[PMLpointE[ii].grade] * (fF.Hx[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hx[PMLpointE[ii].i][PMLpointE[ii].j-1][PMLpointE[ii].k]);
		fF.Ez[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ez[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - CpmlE * PMLpointE[ii].Q1;

		PMLpointE[ii].Q2 = byE1[PMLpointE[ii].grade] * PMLpointE[ii].Q2 + cyE1[PMLpointE[ii].grade] * (fF.Hz[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hz[PMLpointE[ii].i][PMLpointE[ii].j-1][PMLpointE[ii].k]);
		fF.Ex[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ex[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] + CpmlE * PMLpointE[ii].Q2;
	}
#pragma omp parallel for default(shared) private(ii)
	for (ii = iiyE1; ii < iiyE2; ++ii){
		PMLpointE[ii].Q1 = byE2[PMLpointE[ii].grade] * PMLpointE[ii].Q1 + cyE2[PMLpointE[ii].grade] * (fF.Hx[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hx[PMLpointE[ii].i][PMLpointE[ii].j-1][PMLpointE[ii].k]);
		fF.Ez[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ez[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - CpmlE * PMLpointE[ii].Q1;

		PMLpointE[ii].Q2 = byE2[PMLpointE[ii].grade] * PMLpointE[ii].Q2 + cyE2[PMLpointE[ii].grade] * (fF.Hz[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hz[PMLpointE[ii].i][PMLpointE[ii].j-1][PMLpointE[ii].k]);
		fF.Ex[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ex[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] + CpmlE * PMLpointE[ii].Q2;
	}

#pragma omp parallel for default(shared) private(ii)
	for (ii = iiyE2; ii < iizE1; ++ii){
		PMLpointE[ii].Q1 = bzE1[PMLpointE[ii].grade] * PMLpointE[ii].Q1 + czE1[PMLpointE[ii].grade] * (fF.Hy[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hy[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k-1]);
		fF.Ex[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ex[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - CpmlE * PMLpointE[ii].Q1;

		PMLpointE[ii].Q2 = bzE1[PMLpointE[ii].grade] * PMLpointE[ii].Q2 + czE1[PMLpointE[ii].grade] * (fF.Hx[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hx[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k-1]);
		fF.Ey[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ey[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] + CpmlE * PMLpointE[ii].Q2;
	}
#pragma omp parallel for default(shared) private(ii)
	for (ii = iizE1; ii < iizE2; ++ii){
		PMLpointE[ii].Q1 = bzE2[PMLpointE[ii].grade] * PMLpointE[ii].Q1 + czE2[PMLpointE[ii].grade] * (fF.Hy[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hy[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k-1]);
		fF.Ex[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ex[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - CpmlE * PMLpointE[ii].Q1;

		PMLpointE[ii].Q2 = bzE2[PMLpointE[ii].grade] * PMLpointE[ii].Q2 + czE2[PMLpointE[ii].grade] * (fF.Hx[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] - fF.Hx[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k-1]);
		fF.Ey[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] = fF.Ey[PMLpointE[ii].i][PMLpointE[ii].j][PMLpointE[ii].k] + CpmlE * PMLpointE[ii].Q2;
	}
}

void CFS_PML::UpdatePMLHFields(Field &fF)
{
	int ii;
#pragma omp parallel for default(shared) private(ii)
	for (ii = 0; ii < iixH1; ++ii){
		PMLpointH[ii].Q1 = bxH1[PMLpointH[ii].grade] * PMLpointH[ii].Q1 + cxH1[PMLpointH[ii].grade] * (fF.Ez[PMLpointH[ii].i+1][PMLpointH[ii].j][PMLpointH[ii].k] - fF.Ez[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hy[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hy[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] + CpmlE * PMLpointH[ii].Q1;

		PMLpointH[ii].Q2 = bxH1[PMLpointH[ii].grade] * PMLpointH[ii].Q2 + cxH1[PMLpointH[ii].grade] * (fF.Ey[PMLpointH[ii].i+1][PMLpointH[ii].j][PMLpointH[ii].k] - fF.Ey[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hz[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hz[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] - CpmlE * PMLpointH[ii].Q2;
	}
#pragma omp parallel for default(shared) private(ii)
	for (ii = iixH1; ii < iixH2; ++ii){
		PMLpointH[ii].Q1 = bxH2[PMLpointH[ii].grade] * PMLpointH[ii].Q1 + cxH2[PMLpointH[ii].grade] * (fF.Ez[PMLpointH[ii].i+1][PMLpointH[ii].j][PMLpointH[ii].k] - fF.Ez[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hy[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hy[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] + CpmlE * PMLpointH[ii].Q1;

		PMLpointH[ii].Q2 = bxH2[PMLpointH[ii].grade] * PMLpointH[ii].Q2 + cxH2[PMLpointH[ii].grade] * (fF.Ey[PMLpointH[ii].i+1][PMLpointH[ii].j][PMLpointH[ii].k] - fF.Ey[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hz[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hz[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] - CpmlE * PMLpointH[ii].Q2;
	}

#pragma omp parallel for default(shared) private(ii)
	for (ii = iixH2; ii < iiyH1; ++ii){
		PMLpointH[ii].Q1 = byH1[PMLpointH[ii].grade] * PMLpointH[ii].Q1 + cyH1[PMLpointH[ii].grade] * (fF.Ez[PMLpointH[ii].i][PMLpointH[ii].j+1][PMLpointH[ii].k] - fF.Ez[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hx[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hx[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] - CpmlE * PMLpointH[ii].Q1;

		PMLpointH[ii].Q2 = byH1[PMLpointH[ii].grade] * PMLpointH[ii].Q2 + cyH1[PMLpointH[ii].grade] * (fF.Ex[PMLpointH[ii].i][PMLpointH[ii].j+1][PMLpointH[ii].k] - fF.Ex[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hz[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hz[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] + CpmlE * PMLpointH[ii].Q2;
	}
#pragma omp parallel for default(shared) private(ii)
	for (ii = iiyH1; ii < iiyH2; ++ii){
		PMLpointH[ii].Q1 = byH2[PMLpointH[ii].grade] * PMLpointH[ii].Q1 + cyH2[PMLpointH[ii].grade] * (fF.Ez[PMLpointH[ii].i][PMLpointH[ii].j+1][PMLpointH[ii].k] - fF.Ez[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hx[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hx[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] - CpmlE * PMLpointH[ii].Q1;

		PMLpointH[ii].Q2 = byH2[PMLpointH[ii].grade] * PMLpointH[ii].Q2 + cyH2[PMLpointH[ii].grade] * (fF.Ex[PMLpointH[ii].i][PMLpointH[ii].j+1][PMLpointH[ii].k] - fF.Ex[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hz[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hz[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] + CpmlE * PMLpointH[ii].Q2;
	}

#pragma omp parallel for default(shared) private(ii)
	for (ii = iiyH2; ii < iizH1; ++ii){
		PMLpointH[ii].Q1 = bzH1[PMLpointH[ii].grade] * PMLpointH[ii].Q1 + czH1[PMLpointH[ii].grade] * (fF.Ex[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k+1] - fF.Ex[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hy[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hy[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] - CpmlE * PMLpointH[ii].Q1;

		PMLpointH[ii].Q2 = bzH1[PMLpointH[ii].grade] * PMLpointH[ii].Q2 + czH1[PMLpointH[ii].grade] * (fF.Ey[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k+1] - fF.Ey[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hx[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hx[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] + CpmlE * PMLpointH[ii].Q2;
	}
#pragma omp parallel for default(shared) private(ii)
	for (ii = iizH1; ii < iizH2; ++ii){
		PMLpointH[ii].Q1 = bzH2[PMLpointH[ii].grade] * PMLpointH[ii].Q1 + czH2[PMLpointH[ii].grade] * (fF.Ex[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k+1] - fF.Ex[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hy[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hy[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] - CpmlE * PMLpointH[ii].Q1;

		PMLpointH[ii].Q2 = bzH2[PMLpointH[ii].grade] * PMLpointH[ii].Q2 + czH2[PMLpointH[ii].grade] * (fF.Ey[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k+1] - fF.Ey[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k]);
		fF.Hx[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] = fF.Hx[PMLpointH[ii].i][PMLpointH[ii].j][PMLpointH[ii].k] + CpmlE * PMLpointH[ii].Q2;
	}
}

