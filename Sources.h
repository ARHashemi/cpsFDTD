#ifndef SOURCES_H
#define SOURCES_H

#include "Constants.h"




class Source
{
public:
	Source(double deltat, char *fileName);
	~Source();

	double Pulse(double Amplitude, double ntimestep, double location);
	double lambda, nu, T, delay, FWHM;

	double ContinuousWave(double Amplitude, double ntimestep, double location);

private:
	double dt, retard, pinudt;
	int nd, nw, nT, nTrans;
	double *OF;

};




#endif