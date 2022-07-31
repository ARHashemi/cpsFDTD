#include "Ports.h"
#include "stdafx.h"
#include "Constants.h"
#include "Structure.h"
#include "Fields.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

SphericalPort::SphericalPort(Structure &sS, const char *fileName)
{
	cout << BOLDRED << ">" << RESET << "Constructing " << BOLDBLUE << "Port " << RESET << "from " << fileName << " ..";
	dx = sS.dx;
	dy = sS.dy;
	dz = sS.dz;
	Nx = sS.Nx;
	Ny = sS.Ny;
	Nz = sS.Nz;
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> r;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> xc;
	iss >> yc;
	iss >> zc;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> theta1;
	iss >> theta2;
	theta1 = theta1 * pi / 180.0;
	theta2 = theta2 * pi / 180.0;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> phi1;
	iss >> phi2;
	phi1 = phi1 * pi / 180.0;
	phi2 = phi2 * pi / 180.0;
	i1 = (xc + r*cos(phi1)*sin(theta1))/dx;
	j1 = (yc + r*sin(phi1)*sin(theta1))/dy;
	k1 = (zc + r*cos(theta1))/dz;
	dphi = dx/r;// atan(j1*dy/(i1-1)/dx) - phi1;
	dtheta = dz/r;// acos((k1-1)*dz/sqrt(i1*i1*dx*dx+j1*j1*dy*dy+(k1-1)*(k1-1)*dz*dz)) - theta1;
	int ii = 0;
	for (double phi=phi1; phi<=phi2; phi=phi+dphi){
		for (double theta=theta1; theta<=theta2; theta=theta+dtheta){
			ii = ii + 1;
		}
	}
	Num_of_Port_Points = ii + 1;
	port_point = new Port_Points[Num_of_Port_Points];
	ii = 0;
	for (double phi=phi1; phi<=phi2; phi=phi+dphi){
		for (double theta=theta1; theta<=theta2; theta=theta+dtheta){
			port_point[ii].i = (xc + r*sin(theta)*cos(phi))/dx;
			port_point[ii].j = (yc + r*sin(theta)*sin(phi))/dy;
			port_point[ii].k = (zc + r*cos(theta))/dz;
			ii = ii + 1;
		}
	}
	cout << " Done." << endl;
}

void SphericalPort::Update_Fields_on_Port(Field &fF)
{
	for (int ii = 0; ii < Num_of_Port_Points; ++ii){
		port_point[ii].Ex = fF.Ex[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Ey = fF.Ey[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Ez = fF.Ez[port_point[ii].i][port_point[ii].j][port_point[ii].k];

		port_point[ii].Hx = fF.Hx[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Hy = fF.Hy[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Hz = fF.Hz[port_point[ii].i][port_point[ii].j][port_point[ii].k];
	}
}

void SphericalPort::Export_Port(const char *fileName)
{
    ofstream Out(fileName);
    Out << xc/dx << "\t" << yc/dy << "\t" << zc/dz << endl;
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].i << "\t" << port_point[ii].j << "\t" << port_point[ii].k << endl;
    }
    Out.close();
}

void SphericalPort::Export_E_onPort(const char *fileName)
{
	ofstream Out(fileName);
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].Ex << "\t" << port_point[ii].Ey << "\t" << port_point[ii].Ez << endl;
    }
    Out.close();
}
void SphericalPort::Export_E_onPort(ofstream &OutFile)
{
	for (int ii=0; ii<Num_of_Port_Points; ++ii){
		OutFile << port_point[ii].Ex << "\t" << port_point[ii].Ey << "\t" << port_point[ii].Ez << endl;
    }
    //OutFile << endl;
}

void SphericalPort::Export_Ex_onPort(const char *fileName)
{
	ofstream Out(fileName);
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].Ex << "\t";
    }
    Out.close();
}
void SphericalPort::Export_Ex_onPort(ofstream &OutFile)
{
	for (int ii=0; ii<Num_of_Port_Points; ++ii){
		OutFile << port_point[ii].Ex << "\t";
    }
    OutFile << endl;
}

void SphericalPort::Export_Ey_onPort(const char *fileName)
{
	ofstream Out(fileName);
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].Ey << "\t";
    }
    Out.close();
}
void SphericalPort::Export_Ey_onPort(ofstream &OutFile)
{
	for (int ii=0; ii<Num_of_Port_Points; ++ii){
		OutFile << port_point[ii].Ey << "\t";
    }
    OutFile << endl;
}

void SphericalPort::Export_Ez_onPort(const char *fileName)
{
	ofstream Out(fileName);
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].Ez << "\t";
    }
    Out.close();
}
void SphericalPort::Export_Ez_onPort(ofstream &OutFile)
{
	for (int ii=0; ii<Num_of_Port_Points; ++ii){
		OutFile << port_point[ii].Ez << "\t";
    }
    OutFile << endl;
}

void SphericalPort::Export_EIntensity_onPort(const char *fileName)
{
	ofstream Out(fileName);
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].Ex*port_point[ii].Ex + port_point[ii].Ey*port_point[ii].Ey + port_point[ii].Ez*port_point[ii].Ez << "\t";
    }
    Out.close();
}
void SphericalPort::Export_EIntensity_onPort(ofstream &OutFile)
{
	for (int ii=0; ii<Num_of_Port_Points; ++ii){
		OutFile << port_point[ii].Ex*port_point[ii].Ex + port_point[ii].Ey*port_point[ii].Ey + port_point[ii].Ez*port_point[ii].Ez << "\t";
    }
    OutFile << endl;
}

SphericalPort::~SphericalPort()
{

}






PlainPort::PlainPort(Structure &sS, const char *fileName)
{
	cout << BOLDRED << ">" << RESET << "Constructing " << BOLDBLUE << "Port " << RESET << "from " << fileName << " ..";
	dx = sS.dx;
	dy = sS.dy;
	dz = sS.dz;
	Nx = sS.Nx;
	Ny = sS.Ny;
	Nz = sS.Nz;
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> normal_direction;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> position;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> w1;
	iss >> w2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> h1;
	iss >> h2;


	if (normal_direction == 'x'){
		i1 = position/dx; i2 = i1;
		j1 = w1/dy; j2 = w2/dy;
		k1 = h1/dz; k2 = h2/dz;
		Nw = j2 - j1 + 1;
		Nh = k2 - k1 + 1;
		Num_of_Port_Points = (j2-j1+1)*(k2-k1+1) + 1;
		port_point = new Port_Points[Num_of_Port_Points];
		//port_out = new double*[j2-j1+1];
		ii = 0;
		for (int k=k1; k<=k2; k++){
			//port_out[j-j1] = new double[k2-k1+1];
			for (int j=j1; j<=j2; j++){
				//port_out[j-j1][k-k1] = 0.0;
				port_point[ii].i = i1;
				port_point[ii].j = j;
				port_point[ii].k = k;
				port_point[ii].id = sS.ID[i1][j][k];
				ii = ii + 1;
			}
		}
	} else if (normal_direction == 'y'){
		j1 = position/dy; j2 = j1;
		i1 = w1/dx; i2 = w2/dx;
		k1 = h1/dz; k2 = h2/dz;
		Nw = i2 - i1 + 1;
		Nh = k2 - k1 + 1;
		Num_of_Port_Points = (i2-i1+1)*(k2-k1+1) + 1;
		port_point = new Port_Points[Num_of_Port_Points];
		//port_out = new double*[i2-i1+1];
		ii = 0;
		for (int k=k1; k<=k2; k++){
			//port_out[i-i1] = new double[k2-k1+1];
			for (int i=i1; i<=i2; i++){
				//port_out[i-i1][k-k1] = 0.0;
				port_point[ii].i = i;
				port_point[ii].j = j1;
				port_point[ii].k = k;
				port_point[ii].id = sS.ID[i][j1][k];
				ii = ii + 1;
			}
		}
	} else if (normal_direction == 'z'){
		k1 = position/dz; k2 = k1;
		i1 = w1/dx; i2 = w2/dx;
		j1 = h1/dy; j2 = h2/dy;
		Nw = i2 - i1 + 1;
		Nh = j2 - j1 + 1;
		Num_of_Port_Points = (i2-i1+1)*(j2-j1+1) + 1;
		port_point = new Port_Points[Num_of_Port_Points];
		//port_out = new double*[i2-i1+1];
		ii = 0;
		for (int j=j1; j<=j2; j++){
			//port_out[i-i1] = new double[j2-j1+1];
			for (int i=i1; i<=i2; i++){
				//port_out[i-i1][j-j1] = 0.0;
				port_point[ii].i = i;
				port_point[ii].j = j;
				port_point[ii].k = k1;
				port_point[ii].id = sS.ID[i][j][k1];
				ii = ii + 1;
			}
		}
	}
	port_out = new double*[Nw];
	for (int i=0; i<Nw; ++i){
		port_out[i] = new double[Nh];
		for (int j=0; j<Nh; ++j){
			port_out[i][j] = 0.0;
		}
	}


	cout << " Done." << endl;// Nw = " << Nw << " & Nh = " << Nh
}

void PlainPort::Export_ID_onPort(const char *fileName)
{
	ofstream Out(fileName);
	ii = 0;
	for (int j=0; ii<Num_of_Port_Points; ++j){
		for (int i=0; i<Nw; ++i){
			Out << port_point[ii].id << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
    Out.close();
}

short **PlainPort::Export_ID()
{
	short **Out;
	ii = 0;
	Out = new short*[Nw];
	for (int i=0; i<Nw; ++i){
		Out[i] = new short[Nh];
		for (int j=0; j<Nh; ++j){
			Out[i][j] = 0.0;
		}
	}
	for (int j=0; j<Nh; ++j){
		for (int i=0; i<Nw; ++i){
			Out[i][j] = port_point[ii].id;
			ii = ii+1;
		}
    }
    return Out;
}

void PlainPort::Update_Fields_on_Port(Field &fF)
{
	for (int ii = 0; ii < Num_of_Port_Points; ++ii){
		port_point[ii].Ex = fF.Ex[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Ey = fF.Ey[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Ez = fF.Ez[port_point[ii].i][port_point[ii].j][port_point[ii].k];

		port_point[ii].Hx = fF.Hx[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Hy = fF.Hy[port_point[ii].i][port_point[ii].j][port_point[ii].k];
		port_point[ii].Hz = fF.Hz[port_point[ii].i][port_point[ii].j][port_point[ii].k];
	}
}

void PlainPort::Export_Port(const char *fileName)
{
    ofstream Out(fileName);
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].i << "\t" << port_point[ii].j << "\t" << port_point[ii].k << endl;
    }
    Out.close();
}

void PlainPort::Export_E_onPort(const char *fileName)
{
	ofstream Out(fileName);
    for (int ii=0; ii<Num_of_Port_Points; ++ii){
		Out << port_point[ii].Ex << "\t" << port_point[ii].Ey << "\t" << port_point[ii].Ez << endl;
    }
    Out.close();
}
void PlainPort::Export_E_onPort(ofstream &OutFile)
{
	for (int ii=0; ii<Num_of_Port_Points; ++ii){
		OutFile << port_point[ii].Ex << "\t" << port_point[ii].Ey << "\t" << port_point[ii].Ez << endl;
    }
    //OutFile << endl;
}

void PlainPort::Export_Ex_onPort(const char *fileName)
{
	ofstream Out(fileName);
	ii = 0;
	for (int j=0; ii<Num_of_Port_Points; ++j){
		for (int i=0; i<Nw; ++i){
			Out << port_point[ii].Ex << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
    Out.close();
}
void PlainPort::Export_Ex_onPort(ofstream &OutFile)
{
	ii = 0;
	for (int j=0; ii<Num_of_Port_Points; ++j){
		for (int i=0; i<Nw; ++i){
			OutFile << port_point[ii].Ex << "\t";
			ii = ii+1;
		}
		//OutFile << endl;
    }
    OutFile << endl;
}
double **PlainPort::Export_Ex()
{
	ii = 0;
	for (int j=0; j<Nh; ++j){
		for (int i=0; i<Nw; ++i){
			port_out[i][j] = port_point[ii].Ex;
			ii = ii+1;
		}
    }
    return port_out;
}

void PlainPort::Export_Ey_onPort(const char *fileName)
{
	ofstream Out(fileName);
	ii = 0;
	for (int j=0; ii<Num_of_Port_Points; ++j){
		for (int i=0; i<Nw; ++i){
			Out << port_point[ii].Ey << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
    Out.close();
}
void PlainPort::Export_Ey_onPort(ofstream &OutFile)
{
	ii = 0;
    for (int j=0; ii<Num_of_Port_Points; ++j){
		for (int i=0; i<Nw; ++i){
			OutFile << port_point[ii].Ey << "\t";
			ii = ii+1;
		}
		//OutFile << endl;
    }
    OutFile << endl;
}
double **PlainPort::Export_Ey()
{
	ii = 0;
	for (int j=0; j<Nh; ++j){
		for (int i=0; i<Nw; ++i){
			port_out[i][j] = port_point[ii].Ey;
			ii = ii+1;
		}
    }
    return port_out;
}

void PlainPort::Export_Ez_onPort(const char *fileName)
{
	ofstream Out(fileName);
    ii = 0;
    for (int j=0; ii<Num_of_Port_Points; ++j){
		for (int i=0; i<Nw; ++i){
			Out << port_point[ii].Ez << "\t";
			ii = ii+1;
		}
		Out << endl;
    }
    Out.close();
}
void PlainPort::Export_Ez_onPort(ofstream &OutFile)
{
	ii = 0;
    for (int j=0; ii<Num_of_Port_Points; ++j){
		for (int i=0; i<Nw; ++i){
			OutFile << port_point[ii].Ez << "\t";
			ii = ii+1;
		}
//		OutFile << endl;
    }
    OutFile << endl;
}
double **PlainPort::Export_Ez()
{
	ii = 0;
	for (int j=0; j<Nh; ++j){
		for (int i=0; i<Nw; ++i){
			port_out[i][j] = port_point[ii].Ez;
			ii = ii+1;
		}
    }
    return port_out;
}

void PlainPort::Export_EIntensity_onPort(const char *fileName)
{
	ofstream Out(fileName);
    ii = 0;
    for (int j=0; ii<Num_of_Port_Points; ++j){
		for (int i=0; i<Nw; ++i){
			Out << port_point[ii].Ex*port_point[ii].Ex + port_point[ii].Ey*port_point[ii].Ey + port_point[ii].Ez*port_point[ii].Ez << "\t";
			ii = ii+1;
		}
    }Out << endl;
}
void PlainPort::Export_EIntensity_onPort(ofstream &OutFile)
{
	ii = 0;
    for (int j=0; ii<Num_of_Port_Points; ++j){
		for (int i=0; i<Nw; ++i){
			OutFile << port_point[ii].Ex*port_point[ii].Ex + port_point[ii].Ey*port_point[ii].Ey + port_point[ii].Ez*port_point[ii].Ez << "\t";
			ii = ii+1;
		}
    }OutFile << endl;
}
double **PlainPort::Export_EIntensity()
{
	ii = 0;
	for (int j=0; j<Nh; ++j){
		for (int i=0; i<Nw; ++i){
			port_out[i][j] = port_point[ii].Ex*port_point[ii].Ex + port_point[ii].Ey*port_point[ii].Ey + port_point[ii].Ez*port_point[ii].Ez;
			ii = ii+1;
		}
    }
    return port_out;
}

PlainPort::~PlainPort()
{

}
