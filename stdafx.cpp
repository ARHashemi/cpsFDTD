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

// stdafx.cpp : source file that includes just the standard includes
// TopLevelFDTD_OO.pch will be the pre-compiled header
// stdafx.obj will contain the pre-compiled type information

#include "stdafx.h"
#include <iostream>
#include <limits>
void PressEnterToContinue()
  {
  std::cout << std::endl << "Press ENTER to continue... " << std::flush;
  std::cin.ignore( std::numeric_limits <std::streamsize> ::max(), '\n' );
  }

// TODO: reference any additional headers you need in STDAFX.H
// and not in this file

