#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <iomanip>
#include <istream>
#include <fstream>
#include<string.h>
//#include<stdlib.h>
#include <time.h>
//#include <process.h>
#include <omp.h>
#include <ctime>
//#include <chrono>
//#include <fcntl.h>

using namespace std;
//#include <opencv2/opencv.hpp>
//#include "opencv2/core/core_c.h"
//#include "opencv2/core/core.hpp"
//#include "opencv2/flann/miniflann.hpp"
//#include "opencv2/imgproc/imgproc_c.h"
//#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/video/video.hpp"
//#include "opencv2/features2d/features2d.hpp"
//#include "opencv2/objdetect/objdetect.hpp"
//#include "opencv2/calib3d/calib3d.hpp"
//#include "opencv2/ml/ml.hpp"
//#include "opencv2/highgui/highgui_c.h"
//#include "opencv2/highgui/highgui.hpp"
////#include "opencv2/contrib/contrib.hpp"

//using namespace cv;

#include "CImg.h"
using namespace cimg_library;

#include "Constants.h"
#include "Fields.h"
#include "Structure.h"
#include "BasicFDTD.h"
#include "CPML.h"
#include "OutImage.h"
//#include "IOimage.h"
#include "Ports.h"
#include "TF_SF.h"
#include "DL_model.h"
#include "Fluorescence.h"

const char material_file[20] = "material.txt";
const char mesh_file[30] = "structure.txt";
const char plainport1_file[20] = "plainport1.txt";
const char plainport2_file[20] = "plainport2.txt";
const char plainport3_file[20] = "plainport3.txt";
const char plainport4_file[20] = "plainport4.txt";
const char lineport1_file[20] = "lineportx.txt";
const char lineport4_file[20] = "lineportx2.txt";
const char lineport2_file[20] = "lineporty.txt";
const char lineport3_file[20] = "lineportz.txt";//"lineportx2.txt";
const char sphericalport1_file[20] = "spherical_port1.txt";
const char sphericalport2_file[20] = "spherical_port2.txt";
const char sphericalportin1_file[22] = "spherical_portin1.txt";
const char sphericalportin2_file[22] = "spherical_portin2.txt";
char cpml_file[20] = "CFS_PML.txt";
char tfsf_file[20] = "TF_SF.txt";
char DL_file[20] = "DrudeLorentz.txt";
char flourescence_file[20] = "Fluorescence.txt";

double PoutTransient(double time, double T, double pinu)
{
	return exp(-(time-5.0*T)*(time-5*T)/(4.0*T*T)) * sin(pinu*time);//
}

double PoutCW(double time, double pinu)
{
	return sin(pinu*time);//
}

double GaussianSpatialDistr(int r2, int R2)
{
	return exp(-r2/R2);// 1;
}

int main()
{
	string doubleline = "\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550";
	string singleline = "\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500";
	cout << endl << "\t" << BOLDGREEN << "\u2554" << doubleline << doubleline << doubleline << "\u2557" << endl;
	cout << "\t" << "\u2551" << RESET << "\t ARH 3D FDTD \t" << BOLDGREEN << " \u2551" << endl;
	cout << "\t" << "\u255A" << doubleline << doubleline << doubleline << "\u255D" << RESET << endl;

	ofstream logOut;
	logOut.open("logFDTD.txt");
	logOut << endl << "\t" << "\u2554" << doubleline << doubleline << doubleline << "\u2557" << endl;
	logOut << "\t" << "\u2551" << "\t ARH 3D FDTD \t" << " \u2551" << endl;
	logOut << "\t" << "\u255A" << doubleline << doubleline << doubleline << "\u255D" << endl;

    Structure S(material_file);
	S.generate_mesh(mesh_file);
	S.InitializeBasicFDTDCoefficients();

	Field F(S, "TEM");


	CPML pml(S, cpml_file);
	pml.InitializeCPML(S);

	DL_model DL(S, DL_file);
	cout << DL.epsInf << "\t" << DL.omega_D << "\t" << DL.gamma_D << "\t" << DL.omega_L << "\t" << DL.gamma_L << "\t" << DL.delta_eps_L << endl;

	Fluorescence FL(S, flourescence_file);


	TF_SF tfsf(S, tfsf_file);


	SphericalPort sport1(S,sphericalport1_file);
	SphericalPort sport2(S,sphericalport2_file);
	ofstream sOut1, sOut2, sOut3;
	sOut1.open("Ex_xy_arc.txt");
	sOut2.open("Ey_xy_arc.txt");
	sOut3.open("Ez_xy_arc.txt");
	ofstream sOut4, sOut5, sOut6;
	sOut4.open("Ex_xy2_arc.txt");
	sOut5.open("Ey_xy2_arc.txt");
	sOut6.open("Ez_xy2_arc.txt");
    sport1.Export_Port("SPort1.txt");
	sport2.Export_Port("SPort2.txt");

	PlainPort port1(S, plainport1_file);
	PlainPort port2(S, plainport2_file);
	OutImage img1(port1.Nw, port1.Nh, 1, 0, 0.02, "Ex_yz_");
	OutImage img2(port2.Nw, port2.Nh, 1, 0, 0.02, "Ex_xz_");
	PlainPort port3(S, plainport3_file);
	OutImage img3(port3.Nw, port3.Nh, 1, 0, 0.05, "Ex_xy_");
//	PlainPort port4(S, plainport4_file);
//	OutImage img4(port4.Nw, port4.Nh, 1, 0, 10.0, "Ex_xy2_");

//	//IOimage image1(port1.Nh, port1.Nw, -1.5, 1.5, "Out_yz_Ex_", 1);
//	//img1.structure_plot_save(port1.Export_ID());
	img1.structure_plot_save(port1.Export_ID());
	img2.structure_plot_save(port2.Export_ID());
	img3.structure_plot_save(port3.Export_ID());
//	img4.structure_plot_save(port4.Export_ID());




	cout << "\u250C" << singleline << "\u252C" << singleline << "\u252C" << singleline << "\u2500\u2500\u2500\u2510" << endl;
	cout << "\u2502" << "dx = " << S.dx << "\t \u2502" << "dy =  " << S.dy << " \t  \u2502" << "dt = " << std::setw(13) << std::left << S.dt << " \u2502" << endl;
	cout << "\u251C" << singleline << "\u253C" << singleline << "\u253C" << singleline << "\u2500\u2500\u2500\u2524" << endl;
	cout << "\u2502" << "Nx = " << S.Nx << "\t \u2502" << "Ny =  " << S.Ny << "\t  \u2502" << "Nz = " << std::setw(13) << std::left << S.Nz << " \u2502" << endl;
	cout << "\u2514" << singleline << "\u2534" << singleline << "\u2534" << singleline << "\u2500\u2500\u2500\u2518" << endl;
	cout << "\u250C" << singleline << "\u2500\u2500\u2500\u252C" << singleline << "\u252C" << singleline << "\u2510" << endl;
	cout << "\u2502" << "lamda = " << tfsf.lambda << "   \u2502" << "nu = " << std::setw(11) << tfsf.nu << "\u2502" << "Nt = " << std::setw(11) << std::left << tfsf.Nt << "\u2502" << endl;
	cout << "\u2514" << singleline << "\u2500\u2500\u2500\u2534" << singleline << "\u2534" << singleline << "\u2518" << endl;

	logOut << "\u250C" << singleline << "\u252C" << singleline << "\u252C" << singleline << "\u2500\u2500\u2500\u2510" << endl;
	logOut << "\u2502" << "dx = " << S.dx << "\t \u2502" << "dy =  " << S.dy << " \t  \u2502" << "dt = " << std::setw(13) << std::left << S.dt << " \u2502" << endl;
	logOut << "\u251C" << singleline << "\u253C" << singleline << "\u253C" << singleline << "\u2500\u2500\u2500\u2524" << endl;
	logOut << "\u2502" << "Nx = " << S.Nx << "\t \u2502" << "Ny =  " << S.Ny << "\t  \u2502" << "Nz = " << std::setw(13) << std::left << S.Nz << " \u2502" << endl;
	logOut << "\u2514" << singleline << "\u2534" << singleline << "\u2534" << singleline << "\u2500\u2500\u2500\u2518" << endl;
	logOut << "\u250C" << singleline << "\u2500\u2500\u2500\u252C" << singleline << "\u252C" << singleline << "\u2510" << endl;
	logOut << "\u2502" << "lamda = " << tfsf.lambda << "   \u2502" << "nu = " << std::setw(11) << tfsf.nu << "\u2502" << "Nt = " << std::setw(11) << std::left << tfsf.Nt << "\u2502" << endl;
	logOut << "\u2514" << singleline << "\u2500\u2500\u2500\u2534" << singleline << "\u2534" << singleline << "\u2518" << endl;



    cout << BOLDRED << "Time-stepping begins..." << RESET;
    time_t now, newnow;
    time(&now);
    logOut << "Time-stepping begined at " << asctime(localtime(&now)) << endl;


    clock_t begin_time = clock();

	BasicHFieldUpdate(F,S);
//	tfsf.y1UpdateHField(0,F,S);//z1
	pml.UpdateCPMLHFields(F);

	BasicEFieldUpdate(F,S);
	DL.EFieldUpdate(F, S);
	FL.UpdateStatesPopulations(F);
//	FL.UpdatePolarizationDensities(F);
//	FL.UpdateElectricField(F);
//	tfsf.y1UpdateEField(0,F,S);//z1
	pml.UpdateCPMLEFields(F);

//	sport1.Update_Fields_on_Port(F);
//    sport2.Update_Fields_on_Port(F);
//    sport1.Export_Ex_onPort(sOut1);
//    sport1.Export_Ey_onPort(sOut2);
//    sport1.Export_Ez_onPort(sOut3);
//    sport2.Export_Ex_onPort(sOut4);
//    sport2.Export_Ey_onPort(sOut5);
//    sport2.Export_Ez_onPort(sOut6);

//    lport1.Update_Fields_on_Port(F);
//    lport2.Update_Fields_on_Port(F);
////	lport3.Update_Fields_on_Port(F);
////	lport4.Update_Fields_on_Port(F);
//	lport1.Export_Ex_onPort(lOut1);
//	lport1.Export_Ey_onPort(lOut2);
//	lport1.Export_Ez_onPort(lOut3);
//	lport2.Export_Ex_onPort(lOut4);
//	lport2.Export_Ey_onPort(lOut5);
//	lport2.Export_Ez_onPort(lOut6);

	port1.Update_Fields_on_Port(F);
	port2.Update_Fields_on_Port(F);
	img1.plot_save(port1.Export_EIntensity(),0);
	img2.plot_save(port2.Export_EIntensity(),0);
	port3.Update_Fields_on_Port(F);
	img3.plot_save(port3.Export_EIntensity(),0);
//	port4.Update_Fields_on_Port(F);
//	img4.plot_save(port4.Export_EIntensity(),0);
	clock_t end_time = clock();
	double elapsed_secs = double(end_time - begin_time) / CLOCKS_PER_SEC;

	cout << endl << "First time-step elapsed time is " << elapsed_secs << " sec." << endl;
	logOut << "First time-step elapsed time is " << elapsed_secs << " sec." << endl;

	int hour, minut, second;
	hour = (int)(tfsf.Nt*elapsed_secs / 3600.0);
	minut = ((int)(tfsf.Nt*elapsed_secs) % 3600) / 60.0;
	second = ((int)(tfsf.Nt*elapsed_secs) % 3600) % 60;
	cout << "Estimated total time for simulation is " << hour << " hour(s) and " << minut << " minut(s) and " << second << " seconds." << endl;
	logOut << "Estimated total time for simulation is " << hour << " hour(s) and " << minut << " minut(s) and " << second << " seconds." << endl;

	int Nx_2 = S.Nx/2.0, Ny_2 = S.Ny/2.0, Nz_2 = S.Nz/2.0;
	int i, j, k, N2, N1, N0;
	k = 25;// S.Nz -

	int i1TFSF = 25, i2TFSF = S.Nx - 25;
	int j1TFSF = 25, j2TFSF = S.Ny - 25;
	int R = 25;

	double eta0 = sqrt(mu0/eps0), Einc, Hinc, P;

	ofstream Out1, Out2, Out3, Out4;
	Out1.open("P.txt");
	Out2.open("N.txt");
	Out3.open("Coh.txt");
	Out4.open("Transi.txt");
//    ofstream lOut1, lOut2, lOut3;
//	lOut1.open("Ex_xz.txt");
//	lOut2.open("Ey_xz.txt");
//	lOut3.open("Ez_xz.txt");
	int tra20, tra21, tra02, tra10r, tra10nr, tra10stim;
    for (int n = 1; n < tfsf.Nt; ++n){
		BasicHFieldUpdate(F,S);
//		for (i=i1TFSF-1; i<i2TFSF; ++i){
//			for (j=j1TFSF; j<j2TFSF; ++j){
//				F.Hy[i][j][k-1] = F.Hy[i][j][k-1] + S.bH[0] * GaussianSpatialDistr(((i-Nx_2)*(i-Nx_2)+(j-Ny_2)*(j-Ny_2)), R*R) * Einc / S.dz;//10.0 ;
//			}
//		}
//		Hinc = 10.0 * PoutTransient((n*S.dt+S.dt), T, pinu) / eta0;
//		tfsf.y1UpdateHField(n,F,S);//z1
		pml.UpdateCPMLHFields(F);

		BasicEFieldUpdate(F,S);
		DL.EFieldUpdate(F, S);
//		FL.UpdateStatesPopulations(F);
//		FL.UpdateElectricField(F);
//		FL.UpdatePolarizationDensities(F);
//		tfsf.y1UpdateEField(n,F,S);//z1
//        F.Ez[Nx_2][Ny_2][Nz_2] = F.Ez[Nx_2][Ny_2][60] + 10.0*tfsf.Einc[n];
//        F.Ey[Nx_2][Ny_2][Nz_2] = F.Ey[Nx_2][Ny_2][60] + 10.0*tfsf.Einc[n];
//        F.Ex[Nx_2][Ny_2][Nz_2] = F.Ex[Nx_2][Ny_2][60] + 10.0*tfsf.Einc[n];
		pml.UpdateCPMLEFields(F);
		//F.Ez[Nx_2][Ny_2][Nz_2] = F.Ez[Nx_2][Ny_2][Nz_2] + 10.0*sin(2.0*pi*tfsf.nu*S.dt);
//		for (i=i1TFSF; i<(i2TFSF-1); ++i){
//			for (j=j1TFSF; j<j2TFSF; ++j){
//				F.Ex[i][j][k] = F.Ex[i][j][k] + S.bE[0] * GaussianSpatialDistr(((i-Nx_2)*(i-Nx_2)+(j-Ny_2)*(j-Ny_2)), R*R) * Hinc / S.dz;//10.0 ;
//			}
//		}
//		Einc = 10.0 * PoutTransient((n*S.dt+S.dt/2.0), T, pinu);

//		port1.Update_Fields_on_Port(F);
//		port2.Update_Fields_on_Port(F);
//		img1.plot_save(port1.Export_EIntensity(),n);
////		image1.plot_save(port1.Export_Ex(),n);
//		img2.plot_save(port2.Export_EIntensity(),n);
//		port3.Update_Fields_on_Port(F);
//		img3.plot_save(port3.Export_EIntensity(),n);
//		port3.Export_Ex_onPort(lOut1);
//		port3.Export_Ey_onPort(lOut2);
//		port3.Export_Ez_onPort(lOut3);
//		port4.Update_Fields_on_Port(F);
//		img4.plot_save(port4.Export_EIntensity(),n);
		cout << "\r" << n << " of " << tfsf.Nt;
		N2 = 0;
		N1 = 0;
		N0 = 0;
		tra20 = 0;
		tra21 = 0;
		tra02 = 0;
		tra10r = 0;
		tra10nr = 0;
		tra10stim = 0;
//		P = 0.0;
		for(i = 0; i < FL.NumSources; ++i){
			//P = P + FL.FluoMole[i].P;
			N2 = N2 + FL.FluoMole[i].N2;
			N1 = N1 + FL.FluoMole[i].N1;
			N0 = N0 + FL.FluoMole[i].N0;
			tra02 = tra02 + FL.FluoMole[i].trans02;
			tra21 = tra21 + FL.FluoMole[i].trans21;
			tra20 = tra20 + FL.FluoMole[i].trans20;
			tra10r = tra10r + FL.FluoMole[i].trans10r;
			tra10nr = tra10nr + FL.FluoMole[i].trans10nr;
			tra10stim = tra10stim + FL.FluoMole[i].trans10stim;
		}
		Out1 << FL.FluoMole[30].Pe << "\t" << FL.FluoMole[100].Pe << "\t" << FL.FluoMole[1000].Pe;// << "\t" << FL.FluoMole[10000].P << "\t" << FL.FluoMole[40000].P;// << "\t" << FL.FluoMole[400000].P;
//		//Out1 << P;// << "\t";
		Out1 << endl;
		Out2 << N0 << "\t" << N1 << "\t" << N2 << endl;
        Out4 << tra02 << "\t" << tra20 << "\t" << tra21 << "\t" << tra10r << "\t" << tra10nr << "\t" << tra10stim << endl;
		sport1.Update_Fields_on_Port(F);
		sport2.Update_Fields_on_Port(F);
		sport1.Export_Ex_onPort(sOut1);
		sport1.Export_Ey_onPort(sOut2);
		sport1.Export_Ez_onPort(sOut3);
		sport2.Export_Ex_onPort(sOut4);
		sport2.Export_Ey_onPort(sOut5);
		sport2.Export_Ez_onPort(sOut6);
////		sportin1.Update_Fields_on_Port(F);
////		sportin2.Update_Fields_on_Port(F);
////		sportin1.Export_EIntensity_onPort(sOutin1);
////		sportin2.Export_EIntensity_onPort(sOutin2);
//		lport1.Update_Fields_on_Port(F);
//		lport2.Update_Fields_on_Port(F);
//////		lport3.Update_Fields_on_Port(F);
//////		lport4.Update_Fields_on_Port(F);
//		lport1.Export_Ex_onPort(lOut1);
//        lport1.Export_Ey_onPort(lOut2);
//        lport1.Export_Ez_onPort(lOut3);
//        lport2.Export_Ex_onPort(lOut4);
//        lport2.Export_Ey_onPort(lOut5);
//        lport2.Export_Ez_onPort(lOut6);
//		for(int k = 10000; k <= FL.NumSources; k=k+10000){
//			Out3 << FL.PhazeSum(k) << "\t";
//		}
//		Out3 << endl;
    }
    Out1.close();
    Out2.close();
    Out3.close();
    Out4.close();
    sOut1.close();
    sOut2.close();
    sOut3.close();
    sOut4.close();
    sOut5.close();
    sOut6.close();
//	sOutin1.close();
//    sOutin2.close();
//    lOut1.close();
//    lOut2.close();
//    lOut3.close();
//    lOut4.close();
//    lOut5.close();
//    lOut6.close();

//    for (int n = nd; n < Nt; ++n){
//		BasicEFieldUpdate(F,S);
//		//DL.EFieldUpdate(F, S);
//		pml.UpdateCPMLEFields(F);
//		//tfsf.UpdateFields(n,F,S);
//		for (i=i1TFSF; i<(i2TFSF-1); ++i){
//			for (j=j1TFSF; j<j2TFSF; ++j){
//				F.Ex[i][j][k] = F.Ex[i][j][k] + S.bE[0] * GaussianSpatialDistr(((i-Nx_2)*(i-Nx_2)+(j-Ny_2)*(j-Ny_2)), R*R) * Hinc / S.dz;//10.0 ;
//			}
//		}
//		Einc = 10.0 * PoutCW(n*S.dt+S.dt/2.0, pinu);
//		BasicHFieldUpdate(F,S);
//		for (i=i1TFSF-1; i<i2TFSF; ++i){
//			for (j=j1TFSF; j<j2TFSF; ++j){
//				F.Hy[i][j][k-1] = F.Hy[i][j][k-1] + S.bH[0] * GaussianSpatialDistr(((i-Nx_2)*(i-Nx_2)+(j-Ny_2)*(j-Ny_2)), R*R) * Einc / S.dz;//10.0 ;
//			}
//		}
//		Hinc = 10.0 * PoutCW(n*S.dt+S.dt, pinu) / eta0;
//		pml.UpdateCPMLHFields(F);
//		port1.Update_Fields_on_Port(F);
//		//port2.Update_Fields_on_Port(F);
//		img1.plot_save(port1.Export_Ex(),n);
//		//image1.plot_save(port1.Export_Ex(),n);
//		//img2.plot_save(port2.Export_Ex(),n);
//		cout << "\r" << n << " of " << Nt;
//    }
	time(&newnow);
    logOut << "Time-stepping finished at " << asctime(localtime(&newnow)) << endl;
    elapsed_secs = difftime(now,newnow);
    hour = (int)(elapsed_secs / 3600.0);
	minut = ((int)(elapsed_secs) % 3600) / 60.0;
	second = ((int)(elapsed_secs) % 3600) % 60;
    logOut << "Total simulation time-stepping took " << hour << " hour(s) and " << minut << " minut(s) and " << second << " seconds." << endl;
	PressEnterToContinue();

}



