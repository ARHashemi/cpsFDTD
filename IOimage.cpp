#include "IOimage.h"
#include <stdio.h>
#include <iostream>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

IOimage::IOimage(int width, int height, double cmap_min, double cmap_max, const char *file_prefix, int save_modulus=1)
{
    Nx = width;
    Ny = height;
    cmapmin = cmap_min;
	cmapmax = cmap_max;
	cmapmid = (cmapmax+cmapmin)/2.0;
	fileprefix = file_prefix;
	out2D.create(Ny, Nx, CV_8UC3);
	img.create(Ny, Nx, CV_8UC3);

	colmap.create(256, 1, CV_8UC3);
	Vec3b bgr = colmap.at<Vec3b>(0,0);
    for (int i=0; i<128; ++i){
        bgr[0] = 2*i;
        bgr[1] = 2*i;
        bgr[2] = 255;
        colmap.at<Vec3b>(i,0) = bgr;
    }
    for (int i=128; i<256; ++i){
        bgr[0] = 255;
        bgr[1] = 255 - 2*(i-128);
        bgr[2] = 255 - 2*(i-128);
        colmap.at<Vec3b>(i,0) = bgr;
    }
    imwrite("colormap.png", colmap);

    namedWindow(fileprefix, WINDOW_AUTOSIZE );
}

IOimage::~IOimage()
{
	//dtor
}

void IOimage::plot_save(double **dD, int time_step)
{
	if ((time_step % savemodulus) == 0){
		//out2D = Mat(Ny,Nx,CV_32F,dD);
		//img = out2D;
		//LUT(out2D, colmap, img);
		//sprintf(filename, "%s%06d.png", "Ex_", time_step/savemodulus);//fileprefix
		//imshow(fileprefix, img);
		//imwrite(filename, img);
		cout << "my book?!"<< endl;
	}
}

void IOimage::structure_plot(short **dD)
{

}
