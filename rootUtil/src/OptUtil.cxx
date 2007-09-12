// -*- Mode: c++ -*-
#ifndef OptUtil_cxx
#define OptUtil_cxx
/*
* Project: GLAST
* Package: rootUtil
*    File: $Id$
* Authors:
*   EC, Eric Charles,    SLAC              echarles@slac.stanford.edu
*
* Copyright (c) 2007
*                   Regents of Stanford University. All rights reserved.
*
*
*
*/

// This Class's header
#include "rootUtil/OptUtil.h"

// c++/stl headers
#include <string>

bool OptUtil::has_option(const char* options, char opt) {
  // Checks to see if the character 'opt' occurs in the string 'options'
  //
  std::string os(options);
  return (os.find(opt) != os.npos);
}  

#endif
