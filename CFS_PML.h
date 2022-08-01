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
#ifndef CFS_PML_H
#define CFS_PML_H

#include "Structure.h"
#include "Fields.h"

struct PML_point
{
	int i, j, k, grade;
	double Q1, Q2;
};

class CFS_PML
{
public:
	CFS_PML(Structure &sS, char *fileName);
	~CFS_PML();
	void Initialize_PML(Structure &sS);
	void UpdatePMLEFields(Field &fF);
	void UpdatePMLHFields(Field &fF);

	PML_point *PMLpointH, *PMLpointE;
	double p_order, kapmax, alpmin, alpmax, malp, sig_ratio, sigmax;
	int Npmlx1, Npmly1, Npmlz1, Npmlx2, Npmly2, Npmlz2;
	int Num_of_PML_points;

protected:

private:
	void E_pml_Coef(double* &kappa, double* &bcoef, double* &ccoef, bool firstend, int npml, double dt, double del, double m, double ma, double sigmax, double kapmax, double alpmax);
	void H_pml_Coef(double* &kappa, double* &bcoef, double* &ccoef, bool firstend, int npml, double dt, double del, double m, double ma, double sigmax, double kapmax, double alpmax);
	double *bxE1, *cxE1, *byE1, *cyE1, *bzE1, *czE1, *bxH1, *cxH1, *byH1, *cyH1, *bzH1, *czH1;
	double *bxE2, *cxE2, *byE2, *cyE2, *bzE2, *czE2, *bxH2, *cxH2, *byH2, *cyH2, *bzH2, *czH2;
	double *kappa;
	double CpmlE, CpmlH;
	double dx, dy, dz;
	int Nx, Ny, Nz;
	int i, j, k;
	int iixE1, iiyE1, iizE1;
	int iixH1, iiyH1, iizH1;
	int iixE2, iiyE2, iizE2;
	int iixH2, iiyH2, iizH2;
//	int i1Ex, i2Ex, i3Ex, i4Ex, i1Ey, i2Ey, i3Ey, i4Ey, i1Ez, i2Ez, i3Ez, i4Ez;
//	int i1Hx, i2Hx, i3Hx, i4Hx, i1Hy, i2Hy, i3Hy, i4Hy, i1Hz, i2Hz, i3Hz, i4Hz;
};

#endif // CFS_PML_H
