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
