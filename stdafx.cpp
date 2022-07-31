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

