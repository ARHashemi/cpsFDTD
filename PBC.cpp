#include "PBC.h"
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

PBC::PBC(Structure &sS, char *fileName)
{
    cout << BOLDRED << ">" << RESET << "Constructing " << BOLDBLUE << "CPML " << RESET << "boundaries...";
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
}

PBC::~PBC()
{
    //dtor
}
