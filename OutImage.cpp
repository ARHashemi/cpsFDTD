#include "OutImage.h"
#include "CImg.h"
#include <iostream>
using namespace cimg_library;

OutImage::OutImage(int N_x, int N_y, int save_modulus, double cmap_min, double cmap_max, const char *file_name_base)
{
	//CImg<unsigned char> img(Nx,Ny,1,3,0);
	Nx = N_x;
	Ny = N_y;
	savemodulus = save_modulus;
	cmapmin = cmap_min;
	cmapmax = cmap_max;
	cmapmid = (cmapmax+cmapmin)/2.0;
	cmapthird = (cmapmax+cmapmin)/3.0;
	cmap2third = (cmapmax + cmapthird)/2.0;
	filenamebase = file_name_base;
	img.assign(Nx,Ny,1,3);
	//cimg_forXY(img,x,y){
		//img(x,y,0) = 254;
		//img(x,y,1) = 254;
		//img(x,y,2) = 254;
	//}
	main_disp.assign(img); // display it
	color_table = new short*[color_table_length];
	for (int i=0; i<color_table_length; ++i){
		color_table[i] = new short[3];
	}
	color_table[0][0] = 254;color_table[0][1] = 254;color_table[0][2] = 254;
	color_table[1][0] = 254;color_table[1][1] = 0;color_table[1][2] = 0;
	color_table[2][0] = 0;color_table[2][1] = 0;color_table[2][2] = 254;
	color_table[3][0] = 0;color_table[3][1] = 254;color_table[3][2] = 0;
	color_table[4][0] = 0;color_table[4][1] = 0;color_table[4][2] = 0;
	color_table[5][0] = 254;color_table[5][1] = 254;color_table[5][2] = 0;
	color_table[6][0] = 0;color_table[6][1] = 254;color_table[6][2] = 254;
	color_table[7][0] = 254;color_table[7][1] = 0;color_table[7][2] = 254;
	color_table[8][0] = 127;color_table[8][1] = 0;color_table[8][2] = 254;
	color_table[9][0] = 0;color_table[9][1] = 127;color_table[9][2] = 254;
}

OutImage::~OutImage()
{
}

void OutImage::plot()
{
	for (int i = 0; i<255; i++)
      {
	for (int j = 0; j<255; j++)
	  {
              img(i,j,1) = i; // different colour on a pixel - green to max
              img(i,j,2) = j;
              img(254-i,254-j,0) = sqrt(i*j);
	  }
      }
    main_disp.display(img); // display it
    img.save_png("output.png"); // write it
    std::cin.ignore();
}

void OutImage::plot(double **dD, int time_step)
{
	if ((time_step % savemodulus) == 0){

	for (int i = 0; i<Nx; i++)
      {
	for (int j = 0; j<Ny; j++)
	  {
		  temp = ((dD[i][j]-cmapmin)*254.0/(cmapmax-cmapmin));//sqrt(dD[i][j]*cmapmax)
		  if (dD[i][j]<=cmapmin){
            img(i,j,0) = 0;
            img(i,j,1) = 0;
            img(i,j,2) = 0;
		  }else if (dD[i][j]>=cmapmax){
            img(i,j,0) = 254;
			img(i,j,1) = 254;
			img(i,j,2) = 254;
		  }else{
		      if (dD[i][j]<=cmapthird){
                img(i,j,0) = temp;
                img(i,j,1) = 0;
                img(i,j,2) = 0;
		      }else if (dD[i][j]>cmapthird && dD[i][j]<=cmap2third){
                img(i,j,0) = 254;
                img(i,j,1) = temp;
                img(i,j,2) = 0;
		      }else if (dD[i][j]>cmap2third){
                img(i,j,0) = 254;
                img(i,j,1) = 254;
                img(i,j,2) = temp;
		      }
		  }
//		  if (dD[i][j]>=cmapmid && dD[i][j]<cmapmax){
//			  img(i,j,0) = 254;//temp;
//		      img(i,j,1) = 254 - temp;
//              img(i,j,2) = 254 - temp;
//		  }else if(dD[i][j]>=cmapmid && dD[i][j]>=cmapmax){
//			  img(i,j,0) = 254;
//		      img(i,j,1) = 254;
//			  img(i,j,2) = 254;
//		  }else if(dD[i][j]<cmapmid && dD[i][j]>cmapmin){
//			  temp = (int)(dD[i][j]*254.0/cmapmin);//sqrt(dD[i][j]*cmapmin)
//			  img(i,j,0) = 254;//temp;
//			  img(i,j,2) = temp;//254 -
//		      img(i,j,1) = temp;//254 -
//		  }else if(dD[i][j]<cmapmid && dD[i][j]<=cmapmin){
//			  img(i,j,0) = 0;//254;
//			  img(i,j,2) = 0;
//		      img(i,j,1) = 0;
//		  }else{
//			  img(i,j,0) = 0;//254;
//			  img(i,j,1) = 0;//254;
//			  img(i,j,2) = 0;//254;
//		  }
	  }
      }
    main_disp.display(img);
    main_disp.set_title("%s%d",filenamebase,time_step);
//    sprintf(filename, "%s%06d.png", filenamebase, time_step/savemodulus);
//    img.save_png(filename); // write it
    //std::cin.ignore();
}
}

void OutImage::plot_save(double **dD, int time_step)
{
	if ((time_step % savemodulus) == 0){

	for (int i = 0; i<Nx; i++)
      {
	for (int j = 0; j<Ny; j++)
	  {
		  temp = ((dD[i][j]-cmapmin)*254.0/(cmapmax-cmapmin));//sqrt(dD[i][j]*cmapmax)
		  if (dD[i][j]<=cmapmin){
            img(i,j,0) = 0;
            img(i,j,1) = 0;
            img(i,j,2) = 0;
		  }else if (dD[i][j]>=cmapmax){
            img(i,j,0) = 254;
			img(i,j,1) = 254;
			img(i,j,2) = 254;
		  }else{
		      if (dD[i][j]<=cmapthird){
                img(i,j,0) = temp;
                img(i,j,1) = 0;
                img(i,j,2) = 0;
		      }else if (dD[i][j]>cmapthird && dD[i][j]<=cmap2third){
                img(i,j,0) = 254;
                img(i,j,1) = temp;
                img(i,j,2) = 0;
		      }else if (dD[i][j]>cmap2third){
                img(i,j,0) = 254;
                img(i,j,1) = 254;
                img(i,j,2) = temp;
		      }
		  }
//		  if (dD[i][j]>=cmapmid && dD[i][j]<cmapmax){
//			  img(i,j,0) = 254;//temp;
//		      img(i,j,1) = 254 - temp;
//              img(i,j,2) = 254 - temp;
//		  }else if(dD[i][j]>=cmapmid && dD[i][j]>=cmapmax){
//			  img(i,j,0) = 254;
//		      img(i,j,1) = 254;
//			  img(i,j,2) = 254;
//		  }else if(dD[i][j]<cmapmid && dD[i][j]>cmapmin){
//			  temp = (int)(dD[i][j]*254.0/cmapmin);//sqrt(dD[i][j]*cmapmin)
//			  img(i,j,0) = 254;//temp;
//			  img(i,j,2) = temp;//254 -
//		      img(i,j,1) = temp;//254 -
//		  }else if(dD[i][j]<cmapmid && dD[i][j]<=cmapmin){
//			  img(i,j,0) = 0;//254;
//			  img(i,j,2) = 0;
//		      img(i,j,1) = 0;
//		  }else{
//			  img(i,j,0) = 0;//254;
//			  img(i,j,1) = 0;//254;
//			  img(i,j,2) = 0;//254;
//		  }
	  }
      }
    main_disp.display(img);
    main_disp.set_title("%s%d",filenamebase,time_step);
    sprintf(filename, "%s%06d.png", filenamebase, time_step/savemodulus);
    img.save_png(filename); // write it
    //std::cin.ignore();
}
}

void OutImage::structure_plot_save(short **dD)
{
	for (int i = 0; i<Nx; i++){
		for (int j = 0; j<Ny; j++){
			for (int k=0; k<3; ++k){
				img(i,j,k) = color_table[dD[i][j]][k];
			}
		}
	}
    main_disp.display(img);
    sprintf(filename, "%s_structure.png", filenamebase);
    img.save_png(filename); // write it
    //std::cin.ignore();
}
