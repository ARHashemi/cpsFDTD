#include "TF_SF.h"
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

TF_SF::TF_SF(Structure &sS, char *fileName)
{
	cout << BOLDRED << ">" << RESET << "Constructing " << BOLDBLUE << "TF/SF " << RESET << "domains...";
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> lambda;
	k_num = 2.0 * pi / lambda;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> dlambda;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> w0;
	iss >> z0;
	zR = pi * w0*w0 / lambda;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
    iss >> xs;
    iss >> ys;
    iss >> zs;
    i0 = xs/sS.dx;
    j0 = ys/sS.dy;
    z0 = zs/sS.dz;
	//cout << i0 << "\t" << j0 << "\t" << k0 << "\t";
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> theta_inc;
	//cout << theta_inc << "\t";
	theta_inc = theta_inc * pi / 180.0;
	iss >> phi_inc;
	//cout << phi_inc << "\t";
	phi_inc = phi_inc * pi / 180.0;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> E0;
	iss >> psi_E;
	//cout << psi_E << "\t";
	psi_E = psi_E * pi / 180.0;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> nTFSFx1;
	iss >> nTFSFx2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> nTFSFy1;
	iss >> nTFSFy2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> nTFSFz1;
	iss >> nTFSFz2;

	i1 = nTFSFx1;
	i2 = sS.Nx - nTFSFx2;
	j1 = nTFSFy1;
	j2 = sS.Ny - nTFSFy2;
	k1 = nTFSFz1;
	k2 = sS.Nz - nTFSFz2;

	E_phi = E0 * cos(psi_E);
	E_theta = E0 * sin(psi_E);

	eta0 = sqrt(mu0/eps0);

//Elsherbeni
//	E0x = E_theta * cos(theta_inc) * cos(phi_inc) - E_phi * sin(phi_inc);
//	E0y = E_theta * cos(theta_inc) * sin(phi_inc) + E_phi * cos(phi_inc);
//	E0z = -E_theta * sin(theta_inc);
//	H0x = (-1/eta0)*(E_phi * cos(theta_inc) * cos(phi_inc) + E_theta * sin(phi_inc));
//	H0y = (-1/eta0)*(E_phi * cos(theta_inc) * sin(phi_inc) - E_theta * cos(phi_inc));
//	H0z = (1/eta0)*(E_phi * sin(theta_inc));

//Taflove
	E0x = E0 * (cos(psi_E)*sin(phi_inc) - sin(psi_E)*cos(theta_inc)*cos(phi_inc));
	E0y = E0 * (-cos(psi_E)*cos(phi_inc) - sin(psi_E)*cos(theta_inc)*sin(phi_inc));
	E0z = E0 * sin(psi_E)*sin(theta_inc);
	H0x = E0 / eta0 * (sin(psi_E)*sin(phi_inc) + cos(psi_E)*cos(theta_inc)*cos(phi_inc));
	H0y = E0 / eta0 * (-sin(psi_E)*cos(phi_inc) + cos(psi_E)*cos(theta_inc)*sin(phi_inc));
	H0z = -E0 / eta0 * cos(psi_E)*sin(theta_inc);

	kx = sin(theta_inc)*cos(phi_inc);
	ky = sin(theta_inc)*sin(phi_inc);
	kz = cos(theta_inc);
	kabs = sqrt(kx*kx + ky*ky + kz*kz);

	k_dot_r0 = new double[8];
	k_dot_r0[0] = (kx * (i1-i0) * sS.dx + ky * (j1-j0) * sS.dy + kz * (k1-k0) * sS.dz)/C;
	k_dot_r0[1] = (kx * (i1-i0) * sS.dx + ky * (j1-j0) * sS.dy + kz * (k2-k0) * sS.dz)/C;
	k_dot_r0[2] = (kx * (i1-i0) * sS.dx + ky * (j2-j0) * sS.dy + kz * (k1-k0) * sS.dz)/C;
	k_dot_r0[3] = (kx * (i1-i0) * sS.dx + ky * (j2-j0) * sS.dy + kz * (k2-k0) * sS.dz)/C;
	k_dot_r0[4] = (kx * (i2-i0) * sS.dx + ky * (j1-j0) * sS.dy + kz * (k1-k0) * sS.dz)/C;
	k_dot_r0[5] = (kx * (i2-i0) * sS.dx + ky * (j1-j0) * sS.dy + kz * (k2-k0) * sS.dz)/C;
	k_dot_r0[6] = (kx * (i2-i0) * sS.dx + ky * (j2-j0) * sS.dy + kz * (k1-k0) * sS.dz)/C;
	k_dot_r0[7] = (kx * (i2-i0) * sS.dx + ky * (j2-j0) * sS.dy + kz * (k2-k0) * sS.dz)/C;

	l0 = k_dot_r0[0];
	for (int i = 1; i < 8; ++i){
		if (l0 > k_dot_r0[i])
			l0 = k_dot_r0[i];
	}
	//l0 = 0.0;


    num_of_TFSF_points = 2 * ((i2 - i1 + 1)*(j2 - j1 + 1) + (i2 - i1 + 1)*(k2 - k1 + 1) + (j2 - j1 + 1)*(k2 - k1 + 1));
    TFSFp = new TFSF_point[num_of_TFSF_points];

    int i, j, k;
    double x, y, z;
    ii = 0,
    //y-directed
    j = j1;
    for (i=i1; i<i2; ++i){
		for (k=k1; k<k2; ++k){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);
			y = sS.dy * (j-j0-0.5);//-1
			z = sS.dz * (k-k0);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dy;
			y = sS.dy * (j-j0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dy;
			x = sS.dx * (i-i0);//-1
			y = sS.dy * (j-j0-0.5);//-1
			z = sS.dz * (k-k0);
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dy;
			y = sS.dy * (j-j0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dy;
			ii = ii + 1;
		}
    }
    mxz1 = ii;
	k = k2;
    for (i=i1; i<i2; ++i){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);
			y = sS.dy * (j-j0-0.5);//-1
			z = sS.dz * (k-k0);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dy;
			y = sS.dy * (j-j0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dy;
			ii = ii + 1;
    }
    mxz1_1 = ii;
    i = i2;
	for (k=k1; k<k2; ++k){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);//-1
			y = sS.dy * (j-j0-0.5);//-1
			z = sS.dz * (k-k0);
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dy;
			y = sS.dy * (j-j0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dy;
			ii = ii + 1;
    }
    mxz1_2 = ii;

    j = j2;
    for (i=i1; i<i2; ++i){
		for (k=k1; k<k2; ++k){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);
			y = sS.dy * (j-j0+0.5);//+1
			z = sS.dz * (k-k0);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dy;
			y = sS.dy * (j-j0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dy;
			x = sS.dx * (i-i0);//-1
			y = sS.dy * (j-j0+0.5);//+1
			z = sS.dz * (k-k0);
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dy;
			y = sS.dy * (j-j0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dy;
			ii = ii + 1;
		}
    }
    mxz2 = ii;
	k = k2;
    for (i=i1; i<i2; ++i){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);
			y = sS.dy * (j-j0+0.5);//+1
			z = sS.dz * (k-k0);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dy;
			y = sS.dy * (j-j0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dy;
			ii = ii + 1;
    }
    mxz2_1 = ii;
    i = i2;
	for (k=k1; k<k2; ++k){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);//-1
			y = sS.dy * (j-j0+0.5);//+1
			z = sS.dz * (k-k0);
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dy;
			y = sS.dy * (j-j0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dy;
			ii = ii + 1;
    }
    mxz2_2 = ii;

    //x-directed
    i = i1;
    for (j=j1; j<j2; ++j){
		for (k=k1; k<k2; ++k){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0-0.5);//-1
			y = sS.dy * (j-j0);
			z = sS.dz * (k-k0);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dx;
			x = sS.dx * (i-i0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dx;
			x = sS.dx * (i-i0-0.5);//-1
			y = sS.dy * (j-j0);//-1
			z = sS.dz * (k-k0);
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dx;
			x = sS.dx * (i-i0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dx;
			ii = ii + 1;
		}
    }
    mzy1 = ii;
	k = k2;
    for (j=j1; j<j2; ++j){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0-0.5);//-1
			y = sS.dy * (j-j0);
			z = sS.dz * (k-k0);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dx;
			x = sS.dx * (i-i0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dx;
			ii = ii + 1;
    }
    mzy1_1 = ii;
    j = j2;
	for (k=k1; k<k2; ++k){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0-0.5);//-1
			y = sS.dy * (j-j0);//-1
			z = sS.dz * (k-k0);
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dx;
			x = sS.dx * (i-i0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dx;
			ii = ii + 1;
    }
    mzy1_2 = ii;

    i = i2;
    for (j=j1; j<j2; ++j){
		for (k=k1; k<k2; ++k){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0+0.5);//+1
			y = sS.dy * (j-j0);
			z = sS.dz * (k-k0);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dx;
			x = sS.dx * (i-i0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dx;
			x = sS.dx * (i-i0+0.5);//+1
			y = sS.dy * (j-j0);//-1
			z = sS.dz * (k-k0);
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dx;
			x = sS.dx * (i-i0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dx;
			ii = ii + 1;
		}
    }
    mzy2 = ii;
	k = k2;
    for (j=j1; j<j2; ++j){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0+0.5);//+1
			y = sS.dy * (j-j0);
			z = sS.dz * (k-k0);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dx;
			x = sS.dx * (i-i0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dx;
			ii = ii + 1;
    }
    mzy2_1 = ii;
    j = j2;
	for (k=k1; k<k2; ++k){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0+0.5);//+1
			y = sS.dy * (j-j0);//-1
			z = sS.dz * (k-k0);
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dx;
			x = sS.dx * (i-i0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dx;
			ii = ii + 1;
    }
    mzy2_2 = ii;

    //z-directed
    k = k1;
    for (i=i1; i<i2; ++i){
		for (j=j1; j<j2; ++j){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);
			y = sS.dy * (j-j0);//-1
			z = sS.dz * (k-k0-0.5);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dz;
			z = sS.dz * (k-k0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dz;
			x = sS.dx * (i-i0);//-1
			y = sS.dy * (j-j0);
			z = sS.dz * (k-k0-0.5);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dz;
			z = sS.dz * (k-k0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dz;
			ii = ii + 1;
		}
    }
    mxy1 = ii;
	j = j2;
    for (i=i1; i<i2; ++i){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);
			y = sS.dy * (j-j0);//-1
			z = sS.dz * (k-k0-0.5);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dz;
			z = sS.dz * (k-k0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dz;
			ii = ii + 1;
    }
    mxy1_1 = ii;
    i = i2;
	for (j=j1; j<j2; ++j){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);//-1
			y = sS.dy * (j-j0);
			z = sS.dz * (k-k0-0.5);//-1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dz;
			z = sS.dz * (k-k0);//-1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dz;
			ii = ii + 1;
    }
    mxy1_2 = ii;

    k = k2;
    for (i=i1; i<i2; ++i){
		for (j=j1; j<j2; ++j){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);
			y = sS.dy * (j-j0);//-1
			z = sS.dz * (k-k0+0.5);//+1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dz;
			z = sS.dz * (k-k0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dz;
			x = sS.dx * (i-i0);//-1
			y = sS.dy * (j-j0);
			z = sS.dz * (k-k0+0.5);//+1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dz;
			z = sS.dz * (k-k0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dz;
			ii = ii + 1;
		}
    }
    mxy2 = ii;
	j = j2;
    for (i=i1; i<i2; ++i){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);
			y = sS.dy * (j-j0);//-1
			z = sS.dz * (k-k0+0.5);//+1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE1 = TFSFp[ii].k_dot_r_E1 * C / kabs / sS.dz;
			z = sS.dz * (k-k0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H1 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H1 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH1 = TFSFp[ii].k_dot_r_H1 * C / kabs / sS.dz;
			ii = ii + 1;
    }
    mxy2_1 = ii;
    i = i2;
	for (j=j1; j<j2; ++j){
			TFSFp[ii].i = i;
			TFSFp[ii].j = j;
			TFSFp[ii].k = k;
			x = sS.dx * (i-i0);//-1
			y = sS.dy * (j-j0);
			z = sS.dz * (k-k0+0.5);//+1
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_E2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_E2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irE2 = TFSFp[ii].k_dot_r_E2 * C / kabs / sS.dz;
			z = sS.dz * (k-k0);//+1.5
			r2 = x*x + y*y + z*z;
			r_dot_k = (x*kx + y*ky + z*kz) / kabs;
			rho2 = r2 - r_dot_k*r_dot_k;
			r_dot_k = r_dot_k - z0;
            wz = w0 * sqrt(1 + r_dot_k*r_dot_k/(zR*zR));
            Rz = r_dot_k * (1 + (zR*zR)/r_dot_k*r_dot_k);
            phi_z = atan(r_dot_k/zR);
            TFSFp[ii].profile_H2 = exp(-rho2/wz/wz);// * cos(-kabs*rho2/Rz/2.0 + phi_z);// w0/wz *
			TFSFp[ii].k_dot_r_H2 = ((x*kx + y*ky + z*kz)/C - l0) / sS.dt;
//			TFSFp[ii].irH2 = TFSFp[ii].k_dot_r_H2 * C / kabs / sS.dz;
			ii = ii + 1;
    }
    mxy2_2 = ii;

    CEx = sS.dt/(eps0*sS.dx);
    CEy = sS.dt/(eps0*sS.dy);
    CEz = sS.dt/(eps0*sS.dz);
    CHx = sS.dt/(mu0*sS.dx);
    CHy = sS.dt/(mu0*sS.dy);
    CHz = sS.dt/(mu0*sS.dz);


    double eta0 = sqrt(mu0/eps0);
    getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> Nt;
	Einc = new double[Nt+1];
	Hinc = new double[Nt+1];
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> cwmode;
	nu = C / lambda;
	T = 1.0 / nu;
	nT = T / sS.dt;
	double t0 = 1.5*T;//150*sS.dt;//16*T;
	double dnu = C / (lambda - dlambda/2.0) - C / (lambda + dlambda/2.0);
	double domega = 2.0*pi*dnu;
	double sig = domega*domega/(16*log(2));//dnu/(2.0*sqrt(2.0*log(2.0)));
	nd = 5*nT;
	Np = 2.0*nd;
	if (cwmode){
		nd = 5*nT;
		for (int n=0; n<nd; ++n){
			Einc[n] = exp(-((n*sS.dt)-5.0*T)*((n*sS.dt)-5*T)/(4.0*T*T)) * sin(2.0*pi*nu*(n*sS.dt));
			Hinc[n] = exp(-((n*sS.dt+sS.dt/2.0)-5.0*T)*((n*sS.dt+sS.dt/2.0)-5*T)/(4.0*T*T)) * sin(2.0*pi*nu*(n*sS.dt+sS.dt/2.0));
		}
		for (int n=nd; n<=Nt; ++n){
			Einc[n] = sin(2.0*pi*nu*(n*sS.dt));
			Hinc[n] = sin(2.0*pi*nu*(n*sS.dt+sS.dt/2.0));
		}
		ofstream Pout;
		Pout.open("Pout.txt");
		for (int n=0; n<=Nt; ++n){
			Pout << Einc[n] << "\t" << Hinc[n] << endl;
		}
		Pout.close();
	}else{

//		nd = 4.0/sqrt(2.0*sig*sig*pi*pi*dt*dt);
			t0 = 5*T;
		for (int n=0; n<=Nt; ++n){
			Einc[n] = exp(-sig*((double)n*sS.dt-t0)*((double)n*sS.dt-t0))*sin(2.0*pi*nu*((double)n*sS.dt-t0));
			Hinc[n] = exp(-sig*((double)n*sS.dt+sS.dt/2.0-t0)*((double)n*sS.dt+sS.dt/2.0-t0))*sin(2.0*pi*nu*((double)n*sS.dt+sS.dt/2.0-t0));
			//Einc[n] = exp(-((n*sS.dt)-5.0*T)*((n*sS.dt)-5.0*T)/(4.0*T*T)) * sin(2.0*pi*nu*(n*sS.dt));
			//Hinc[n] = exp(-((n*sS.dt+sS.dt/2.0)-5.0*T)*((n*sS.dt+sS.dt/2.0)-5*T)/(4.0*T*T)) * sin(2.0*pi*nu*(n*sS.dt+sS.dt/2.0));

//			Einc[n] = exp(-2.0*sig*sig*pi*pi*sS.dt*sS.dt*(n-nd)*(n-nd)) * sin(2.0*pi*nu*(n*sS.dt));//-nd
//			Hinc[n] = exp(-2.0*sig*sig*pi*pi*sS.dt*sS.dt*(n+0.5-nd)*(n+0.5-nd)) * sin(2.0*pi*nu*(n*sS.dt+sS.dt/2.0));//-nd
			//Ricker
//			Einc[n] = (1.0 - 2.0*pi*pi*nu*nu*(sS.dt*n-t0)*(sS.dt*n-t0)) * exp(-pi*pi*nu*nu*(sS.dt*n-t0)*(sS.dt*n-t0));
//			Hinc[n] = (1.0 - 2.0*pi*pi*nu*nu*(sS.dt*n+sS.dt/2.0-t0)*(sS.dt*n+sS.dt/2.0-t0)) * exp(-pi*pi*nu*nu*(sS.dt*n+sS.dt/2.0-t0)*(sS.dt*n+sS.dt/2.0-t0));
		}
		ofstream Pout;
		Pout.open("Pout.txt");
		for (int n=0; n<=Nt; ++n){
			Pout << Einc[n] << "\t" << Hinc[n] << endl;
		}
		Pout.close();
	}


	cout << " Done." << endl;
}

TF_SF::~TF_SF()
{
	//dtor
}

double TF_SF::Einc_onPoint(double retarded)
{
	if (retarded <= 0)
		return 0.0;
	else{
        ntemp = floor(retarded);
        resid = retarded - (double)ntemp;
       	return ((1.0-resid)*Einc[ntemp]+resid*Einc[ntemp+1]);
	}
}

double TF_SF::Hinc_onPoint(double retarded)
{
	if (retarded <= 0)
		return 0.0;
	else{
        ntemp = floor(retarded);
        resid = retarded - (double)ntemp;
		return ((1.0-resid)*Hinc[ntemp]+resid*Hinc[ntemp+1]);
	}
}

void TF_SF::UpdateEField(int timestep, Field &fF, Structure &sS)
{
	// y-directed
	for (ii = 0; ii<mxz1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHy * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHy * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mxz1; ii<mxz1_1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHy * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mxz1_1; ii<mxz1_2; ++ii){
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHy * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}

	for (ii = mxz1_2; ii<mxz2; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHy * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHy * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mxz2; ii<mxz2_1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHy * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mxz2_1; ii<mxz2_2; ++ii){
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHy * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}

	// x-directed
	for (ii = mxz2_2; ii<mzy1; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mzy1; ii<mzy1_1; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mzy1_1; ii<mzy1_2; ++ii){
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}

	for (ii = mzy1_2; ii<mzy2; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mzy2; ii<mzy2_1; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mzy2_1; ii<mzy2_2; ++ii){
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}

	// z-directed
	for (ii = mzy2_2; ii<mxy1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHz * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHz * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mxy1; ii<mxy1_1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHz * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mxy1_1; ii<mxy1_2; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHz * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}

	for (ii = mxy1_2; ii<mxy2; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHz * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHz * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mxy2; ii<mxy2_1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHz * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mxy2_1; ii<mxy2_2; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHz * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}

}
void TF_SF::x1UpdateEField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mxz2_2; ii<mzy1; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mzy1; ii<mzy1_1; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mzy1_1; ii<mzy1_2; ++ii){
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
}
void TF_SF::x2UpdateEField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mzy1_2; ii<mzy2; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mzy2; ii<mzy2_1; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mzy2_1; ii<mzy2_2; ++ii){
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
}
void TF_SF::z1UpdateEField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mzy2_2; ii<mxy1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1);// * TFSFp[ii].profile_E1
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHz * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2);// * TFSFp[ii].profile_E2
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHz * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mxy1; ii<mxy1_1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1);// * TFSFp[ii].profile_E1
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHz * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mxy1_1; ii<mxy1_2; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2);// * TFSFp[ii].profile_E2
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHz * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
}
void TF_SF::z2UpdateEField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mxy1_2; ii<mxy2; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHz * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHz * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mxy2; ii<mxy2_1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0y * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHz * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mxy2_1; ii<mxy2_2; ++ii){
		fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHz * E0y * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}

}
void TF_SF::y1UpdateEField(int timestep, Field &fF, Structure &sS)
{
	// y-directed
	for (ii = 0; ii<mxz1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHy * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHy * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mxz1; ii<mxz1_1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHy * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mxz1_1; ii<mxz1_2; ++ii){
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHy * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
}
void TF_SF::y2UpdateEField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mxz1_2; ii<mxz2; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHy * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);

		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHy * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
	for (ii = mxz2; ii<mxz2_1; ++ii){
		fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0z * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E1) * TFSFp[ii].profile_E1;
		//fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CHy * E0x * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E1);
	}
	for (ii = mxz2_1; ii<mxz2_2; ++ii){
		fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0x * Hinc_onPoint((double)timestep-TFSFp[ii].k_dot_r_E2) * TFSFp[ii].profile_E2;
		//fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CHy * E0z * Einc_onPoint(timestep+1+0.5-TFSFp[ii].k_dot_r_E2);
	}
}

void TF_SF::UpdateHField(int timestep, Field &fF, Structure &sS)
{
	// y-directed
	for (ii = 0; ii<mxz1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] - CHy * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;

		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] + CHy * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mxz1; ii<mxz1_1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] - CHy * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mxz1_1; ii<mxz1_2; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] + CHy * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}

	for (ii = mxz1_2; ii<mxz2; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hz[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] + CHy * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;

		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] - CHy * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mxz2; ii<mxz2_1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hz[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] + CHy * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mxz2_1; ii<mxz2_2; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] - CHy * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}

	// x-directed
	for (ii = mxz2_2; ii<mzy1; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;

		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mzy1; ii<mzy1_1; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mzy1_1; ii<mzy1_2; ++ii){
		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}

	for (ii = mzy1_2; ii<mzy2; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;

		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hz[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mzy2; ii<mzy2_1; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mzy2_1; ii<mzy2_2; ++ii){
		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hz[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}

	// z-directed
	for (ii = mzy2_2; ii<mxy1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] + CHz * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;

		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] - CHz * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mxy1; ii<mxy1_1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] + CHz * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mxy1_1; ii<mxy1_2; ++ii){
		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] - CHz * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}

	for (ii = mxy1_2; ii<mxy2; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] - CHz * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;

		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] + CHz * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mxy2; ii<mxy2_1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] - CHz * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mxy2_1; ii<mxy2_2; ++ii){
		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] + CHz * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
}
void TF_SF::y1UpdateHField(int timestep, Field &fF, Structure &sS)
{
	// y-directed
	for (ii = 0; ii<mxz1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] - CHy * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;

		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] + CHy * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mxz1; ii<mxz1_1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] - CHy * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mxz1_1; ii<mxz1_2; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j-1][TFSFp[ii].k] + CHy * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
}
void TF_SF::y2UpdateHField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mxz1_2; ii<mxz2; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hz[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] + CHy * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;

		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] - CHy * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mxz2; ii<mxz2_1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEy * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hz[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] + CHy * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mxz2_1; ii<mxz2_2; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEy * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j+1][TFSFp[ii].k] - CHy * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
}
void TF_SF::x1UpdateHField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mxz2_2; ii<mzy1; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;

		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mzy1; ii<mzy1_1; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mzy1_1; ii<mzy1_2; ++ii){
		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i-1][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
}
void TF_SF::x2UpdateHField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mzy1_2; ii<mzy2; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;

		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hz[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mzy2; ii<mzy2_1; ++ii){
		//fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ez[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEx * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hy[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] + CHx * E0z * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mzy2_1; ii<mzy2_2; ++ii){
		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEx * H0z * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hz[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] = fF.Hz[TFSFp[ii].i+1][TFSFp[ii].j][TFSFp[ii].k] - CHx * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
}
void TF_SF::z1UpdateHField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mzy2_2; ii<mxy1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] + CHz * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1);// * TFSFp[ii].profile_H1

		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] - CHz * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2);// * TFSFp[ii].profile_H2
	}
	for (ii = mxy1; ii<mxy1_1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] + CHz * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1);// * TFSFp[ii].profile_H1
	}
	for (ii = mxy1_1; ii<mxy1_2; ++ii){
		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k-1] - CHz * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2);// * TFSFp[ii].profile_H2
	}
}
void TF_SF::z2UpdateHField(int timestep, Field &fF, Structure &sS)
{
	for (ii = mxy1_2; ii<mxy2; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] - CHz * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;

		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] + CHz * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
	for (ii = mxy2; ii<mxy2_1; ++ii){
		//fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ex[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] - CEz * H0y * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H1);
		fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] = fF.Hy[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] - CHz * E0x * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H1) * TFSFp[ii].profile_H1;
	}
	for (ii = mxy2_1; ii<mxy2_2; ++ii){
		//fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] = fF.Ey[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k] + CEz * H0x * Hinc_onPoint(timestep+0.5+1-TFSFp[ii].k_dot_r_H2);
		fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] = fF.Hx[TFSFp[ii].i][TFSFp[ii].j][TFSFp[ii].k+1] + CHz * E0y * Einc_onPoint((double)timestep-TFSFp[ii].k_dot_r_H2) * TFSFp[ii].profile_H2;
	}
}
