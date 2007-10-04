
#ifndef CelEventIDs_h
#define CelEventIDs_h

/*
* Project: GLAST
* Package: rootUtil
*    File: $Id$
* Authors:
*   DC, David Chamont, LLR, chamont@llr.in2p3.fr
*
* Copyright (c) 2007
*                   Regents of Stanford University. All rights reserved.
*
*/

#include "BranchGroup.h"
#include "BgDataHandle.h"

class CelEventIDs : public BranchGroup
 {

  public :

    CelEventIDs() ; 
    CelEventIDs( const CelEventIDs & ) ;
    ~CelEventIDs() ;
    CelEventIDs & operator=( const CelEventIDs & ) ;
    void printInfo( const char * /* options */ ) const ;
    
  private :
  
	BgDataHandle<char *> _productionID ;
	BgDataHandle<Long64_t> _runID ;
	BgDataHandle<Long64_t> _eventID ;

    ClassDef(CelEventIDs,0)
 } ;
 
#endif

