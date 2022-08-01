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
#ifndef BASICFDTD_H
#define BASICFDTD_H

#include "Fields.h"
#include "Structure.h"


void BasicEFieldUpdate(Field &fF, Structure &sS);
void BasicHFieldUpdate(Field &fF, Structure &sS);

//void BasicEFieldUpdateTE(Field &fF, Structure &sS);
//void BasicHFieldUpdateTE(Field &fF, Structure &sS);

#endif
