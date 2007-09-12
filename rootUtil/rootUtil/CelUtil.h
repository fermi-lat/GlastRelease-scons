// -*- Mode: c++ -*-
#ifndef CelUtil_h
#define CelUtil_h
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

namespace CelUtil {
  
  // Merge a Collection of input files into a single pointer skim
  CompositeEventList* mergeCelFiles(TCollection&, const char* ouputfileName=0, const char* option = "CREATE");  
  
  // Merge a Collection of pointer skims
  CompositeEventList* mergeCompositeEventLists(TCollection&, const char* outputfileName=0, const char* option = "CREATE");

}

#endif
