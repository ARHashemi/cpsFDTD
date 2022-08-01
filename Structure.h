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
#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <string>
struct shape
{
    int index;
    int material_index;
    std::string type;
    double theta, phi, psi;
    double x0, y0, z0;
    double a, b, c;
    int ix0, jy0, kz0;
    int ia, jb, kc;
};

class Structure
{
public:
    Structure(const char *fileName);
    ~Structure();
    void generate_mesh(const char *fileName);
    void InitializeBasicFDTDCoefficients();
    double Lx, Ly, Lz;
    double dx, dy, dz;
    double dt;
    int Nx, Ny, Nz;
    short ***ID;
    int NoShapes;
    shape *shapes;
    bool freespace;

    short NoMedia;
	double *epsR;
	double *muR;
	double *sig;
	double *sigM;
	double *aE, *bE, *aH, *bH;
	double *cExy, *cEyx, *cEzx, *cExz, *cEzy, *cEyz;
	double *cHxy, *cHyx, *cHzx, *cHxz, *cHzy, *cHyz;

private:

};

void cuboid(short ***&ID, int Nx, int Ny, int Nz, double dx, double dy, double dz, shape cuboid_shape);
void ellipsoid(short ***&ID, shape ellipsoid_shape, int Nx, int Ny, int Nz, double dx, double dy, double dz);
void cylinder(short ***&ID, shape cylinder_shape, int Nx, int Ny, int Nz, double dx, double dy, double dz);
void multishape(short ***&ID, int Nx, int Ny, int Nz, double dx, double dy, double dz, shape multi_shape, shape single_shape);

#endif // STRUCTURE_H
