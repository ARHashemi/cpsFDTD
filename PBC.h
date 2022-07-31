#ifndef PBC_H
#define PBC_H
#include "Structure.h"
#include "Fields.h"

class PBC
{
    public:
        PBC(Structure &sS, char *fileName);
        virtual ~PBC();

        double kx, ky, kz, k;
        double rotx, roty, rotz;
        double kstep;
        double a1x, a1y, a1z;
        double a2x, a2y, a2z;
        double b1x, b1y, b1z;
        double b2x, b2y, b2z;

    protected:

    private:
};

#endif // PBC_H
