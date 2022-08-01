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
#ifndef OUTIMAGE_H
#define OUTIMAGE_H

#include "CImg.h"
#include <iostream>
using namespace cimg_library;

class OutImage
{
public:
	OutImage(int N_x, int N_y, int  save_modulus, double cmap_min, double cmap_max, const char *file_name_base);
	~OutImage();
	void plot_save(double **dD, int time_step);
	void structure_plot_save(short **dD);
	void plot();
	void plot(double **dD, int time_step);
	CImg<double> img;
	short **color_table;
	int savemodulus, Nx, Ny, color_table_length = 10;
	double cmapmin, cmapmax, cmapmid, cmapthird, cmap2third;
	const char *filenamebase;
	CImgDisplay main_disp;

private:
	int temp, ii;
	char filename[20];
};


#endif
