#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "Structure.h"
#include "Fields.h"



double ContinuousWave(double Amplitude, double retardedtime, int nd, int nw, double pinudt);
double Pulse(double Amplitude, double retardedtime, int nd, int nw, double pinudt);
double GaussianWaveSpatialDistribution(double r, double C1, double C2, double C3);

class Gaussian
{
public:
	Gaussian(Structure &sS, char *fileName);
	~Gaussian();


	void GaussianPulseUpdate_z_x(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_z_x(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_z_y(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_z_y(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_y_x(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_y_x(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_y_z(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_y_z(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_x_y(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_x_y(int timestep, Field &fF, Structure &sS);
	void GaussianPulseUpdate_x_z(int timestep, Field &fF, Structure &sS);
	void GaussianCWUpdate_x_z(int timestep, Field &fF, Structure &sS);
	void InitializeGaussianSource();

	double lambda, nu, T, K;
	double amp, ampH;
	double l0, w0, z0;
	int i1TFSF, i2TFSF, j1TFSF, j2TFSF, k1TFSF, k2TFSF;
	char direction, polarization;

private:
	double dt, dx, dy, dz, r, r01, r02;
	double retard, pinudt, Cl1, Cl2, Cl3, P;
	int nd, nw, nT, nTrans, nstep;
	double **Er, Einc, Hinc;
};

#endif
