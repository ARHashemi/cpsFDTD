#include "stdafx.h"
#include "Constants.h"
#include "Sources.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;



Source::Source(double deltat, char *fileName)
{
	dt = deltat;
	cout << "Constructing Source";
	ifstream InputFile;
	InputFile.open(fileName);
	string sline, id;
	getline(InputFile, sline);
	istringstream iss(sline);
	iss >> lambda;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> FWHM;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> delay;
	nu = C / lambda;
	T = 1.0 / nu;
	nT = T / dt;
	pinudt = 2.0 * pi * nu * dt;
	nw = 2 * nT;
	nd = 4 * nw;
	nTrans = nd + nw;
	cout << " ..." << endl;
}

Source::~Source()
{
}

double Source::ContinuousWave(double Amplitude, double ntimestep, double location)
{
	retard = ntimestep - location / C / dt;
	if (retard<nd){
		return (Amplitude * exp(-(retard - nd)*(retard - nd) / (nw*nw)) * sin(pinudt*retard));
	}
	else{
		return (Amplitude * sin(pinudt*retard));
	}
}

double Source::Pulse(double Amplitude, double ntimestep, double location)
{
	retard = ntimestep - location / C / dt;
	return (Amplitude * exp(-double(retard - nd)*(retard - nd) / (nw*nw)) * sin(pinudt*(retard - nd)));
}