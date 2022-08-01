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
