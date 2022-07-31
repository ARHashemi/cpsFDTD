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
