#include "stdafx.h"
#include "SemiFluorescence.h"
#include "Structure.h"
#include "Fields.h"
#include "Constants.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <time.h> 
#include "Fluorescence.h"
using namespace std;


SemiFluorescence::SemiFluorescence(Structure &sS, char *fileName)
{
	cout << "Reading emitters properties from file \"" << fileName << "\" ";
	ifstream InputFile;
	string sline;
	istringstream iss;
	InputFile.open(fileName);
	getline(InputFile, sline);
	iss.str(sline);
	iss >> lambda;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> EmissionAmp;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> EmissionThershold;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> NumSources;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> i1;
	iss >> i2;
	getline(InputFile, sline);
	iss.clear();
	iss.str(sline);
	iss >> j1;
	iss >> j2;

	dt = sS.dt;
	nu = C / lambda;
	T = 1.0 / nu;
	nT = T / dt;
	pinudt = 2.0 * pi * nu * dt;
	nw = 2 * nT;
	nd = 4 * nw;

	srand(time(NULL));
	idif = i2 - i1;
	jdif = j2 - j1;
	location = new int*[NumSources];
	Amps = new double[NumSources];
	for (i = 0; i < NumSources; ++i){
		location[i] = new int[3];
		location[i][0] = rand() % idif + i1;
		location[i][1] = rand() % jdif + j1;
		location[i][2] = 0;
		Amps[i] = 0.0;
	}
	cout << "..." << endl;

}

void SemiFluorescence::UpdateEmitters(Field &fF)
{
	for (i = 0; i < NumSources; ++i){
		if (location[i][2] > 2*nd){//fF.Ex[location[i][0]][location[i][1]] >= EmissionThershold && 
			location[i][2] = 0;
			Amps[i] = fF.Ex[location[i][0]][location[i][1]] * 14.142135623730950488016887242097;
		}
		fF.Ex[location[i][0]][location[i][1]] = fF.Ex[location[i][0]][location[i][1]] + Emission(Amps[i], location[i][2], nd, nw, pinudt);
		location[i][2] = location[i][2] + 1;
	}
}

//double Emission(double Amplitude, double retardedtime, int nd, int nw, double pinudt)
//{
//	//retard = ntimestep - location / C / dt;
//	return (Amplitude * exp(-double(retardedtime - nd)*(retardedtime - nd) / (nw*nw)) * sin(pinudt*(retardedtime - nd)));
//}


SemiFluorescence::~SemiFluorescence()
{
}
