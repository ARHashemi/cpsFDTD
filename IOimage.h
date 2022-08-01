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
#ifndef IOIMAGE_H
#define IOIMAGE_H

#include <stdio.h>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

class IOimage
{
	public:
		IOimage(int width, int height, double cmap_min, double cmap_max, const char *file_prefix, int save_modulus);
		~IOimage();
		void plot_save(double **dD, int time_step);
		void structure_plot(short **dD);
		int savemodulus, Nx, Ny;
		double cmapmin, cmapmax, cmapmid;
		const char *fileprefix;

		Mat out2D, img, colmap;



	protected:


	private:
		char filename[20];
};

#endif // IOIMAGE_H
