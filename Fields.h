#ifndef FIELDS_H
#define FIELDS_H

#include "Structure.h"
void DynamicMemoryAllocate(double ***D, int sizex, int sizey, int sizez, double initialValue);
void DynamicMemoryDeAllocate(double ***D, int sizex, int sizey);

class Field
{
public:
	Field(Structure &sS, char *smode);
	~Field();
	double ***Ex;
	double ***Ey;
	double ***Ez;
	double ***Hx;
	double ***Hy;
	double ***Hz;
	char *mode;
	char const *TE_mode = "TE";
	char const *TM_mode = "TM";

private:
	int xN;
	int yN;
	int zN;

	//void setSpaceDimensions(int nNx, int nNy);

};



#endif
