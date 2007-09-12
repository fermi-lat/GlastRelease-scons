// -*- Mode: c++ -*-
#ifndef PointerUtil_h
#define PointerUtil_h
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

// ROOT Headers
#include <Rtypes.h>

// forward declares
class TCollection;
class CompositeEventList;

namespace PointerUtil {
  
  // Merge a Collection of input files into a single pointer skim
  CompositeEventList* mergeSkimsFromFiles(TCollection& skimFiles, const char* fileName=0, const char* option = "CREATE");  
  
  // Merge a Collection of pointer skims
  CompositeEventList* mergeSkims(TCollection& skims, const char* fileName=0, const char* option = "CREATE");

}

#endif
