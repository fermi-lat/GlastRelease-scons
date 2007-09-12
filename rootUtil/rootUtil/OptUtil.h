// -*- Mode: c++ -*-
#ifndef OptUtil_h
#define OptUtil_h
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


namespace OptUtil {

  // Checks to see if the character 'opt' appears in the string 'options'
  bool has_option(const char* options, char opt);
}

#endif
