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
