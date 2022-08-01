//    This file is part of the arhFDTD package
//    Copyright (C) <2022>  <AliReza Hashemi>

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as
//    published by the Free Software Foundation, either version 3 of the
//    License, or (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.

//    You should have received a copy of the GNU Affero General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#ifndef TF_SF_H
#define TF_SF_H

#include "Structure.h"
#include "Fields.h"

struct TFSF_point
{
    int i, j, k;
    double k_dot_r_E1, k_dot_r_H1;
    double k_dot_r_E2, k_dot_r_H2;
    double profile_E1, profile_H1;
    double profile_E2, profile_H2;
    //double irE1, irE2, irH1, irH2;
    double E0_E1, E0_E2, E0_H1, E0_H2;
};

class TF_SF
{
	public:
		TF_SF(Structure &sS, char *fileName);
		~TF_SF();
		void UpdateEField(int timestep, Field &fF, Structure &sS);
		void UpdateHField(int timestep, Field &fF, Structure &sS);
		void z1UpdateEField(int timestep, Field &fF, Structure &sS);
		void z1UpdateHField(int timestep, Field &fF, Structure &sS);
		void y1UpdateEField(int timestep, Field &fF, Structure &sS);
		void y1UpdateHField(int timestep, Field &fF, Structure &sS);
		void x1UpdateEField(int timestep, Field &fF, Structure &sS);
		void x1UpdateHField(int timestep, Field &fF, Structure &sS);
		void z2UpdateEField(int timestep, Field &fF, Structure &sS);
		void z2UpdateHField(int timestep, Field &fF, Structure &sS);
		void y2UpdateEField(int timestep, Field &fF, Structure &sS);
		void y2UpdateHField(int timestep, Field &fF, Structure &sS);
		void x2UpdateEField(int timestep, Field &fF, Structure &sS);
		void x2UpdateHField(int timestep, Field &fF, Structure &sS);
		double Einc_onPoint(double retarded);
		double Hinc_onPoint(double retarded);
		double theta_inc, phi_inc;
		double psi_E, E0;
		int nTFSFx1, nTFSFx2, nTFSFy1, nTFSFy2, nTFSFz1, nTFSFz2;
		int num_of_TFSF_points;
		TFSF_point *TFSFp;
		double kx, ky, kz, kabs;
		double *Einc, *Hinc, *Einc_tran, *Hinc_tran;
		double lambda, dlambda, nu, T, amp;
		int nT, nw, nd, Np, Nt;
		double k_num, w0, z0, zR, wz, Rz, phi_z, rho2;
        double xs, ys, zs;
		int i0, j0, k0;
		bool cwmode;

	protected:

	private:
		double E_theta, E_phi;
		double dt, resid;
		int ntemp;
		double *k_dot_r0, l0, r2, r_dot_k;
		double E0x, E0y, E0z, H0x, H0y, H0z;
		double eta0, CEx, CEy, CEz, CHx, CHy, CHz;
		int i1, i2, j1, j2, k1, k2, ii;
		int mxz1, mxz1_1, mxz1_2, mxz2, mxz2_1, mxz2_2, mxy1, mxy1_1, mxy1_2, mxy2, mxy2_1, mxy2_2, mzy1, mzy1_1, mzy1_2, mzy2, mzy2_1, mzy2_2;
};

#endif // TF_SF_H
