
#ifndef CelFileAndTreeOffset_h
#define CelFileAndTreeOffset_h

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
*/



#include "BranchGroup.h"
#include "BgDataHandle.h"

#include <TString.h>
class TChain ;
class TTree ;
class TFile ;
class TObjArray ;
//class TArrayL64 ;

class CelFileAndTreeOffset  : public BranchGroup
 {

  public :

    // construction
	CelFileAndTreeOffset() ; // Needed for ROOT
	CelFileAndTreeOffset( const TString & componentName ) ;
	CelFileAndTreeOffset( const CelFileAndTreeOffset & ) ;
	CelFileAndTreeOffset & operator=( const CelFileAndTreeOffset & ) ;
    virtual ~CelFileAndTreeOffset() ;

    // Access
    inline void reset() { _offset = 0 ; }
    inline void increment( Long64_t offset ) { _offset = _offset + offset ; }
    inline Long64_t getOffset() const { return _offset ; }
    inline const TString & componentName() const { return _componentName ; }

    // print the info about the current set
    void printInfo( const char * options ) const ;  

  private :

    // data
    TString _componentName ;
    BgDataHandle<Long64_t> _offset ;

    ClassDef(CelFileAndTreeOffset,0)
 } ;


#endif

