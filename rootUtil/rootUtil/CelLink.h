
#ifndef CelLink_h
#define CelLink_h

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
#include "DataHandle.h"
#include <TObject.h>
#include <map>
#include <vector>
class TTree ;
class TFile ;
class TChain ;
class TVirtualIndex ;
class TCollection ;
class CelComponent ;
class TObjArray ;

//
// CelLink is made for reading/writing entries of the
// link tree. It associates a given composite event with
// a given entry in the file/tree tree, and also store the
// corresponding offset. The branches are :
// 
// Link_EventIndex  -> Absolute index of current composite event
// Link_MetaIndex   -> Index of the associated entry in the File Tree
// Link_MetaOffset  -> Offset in events for this entry in the File Tree
// 

class CelLink : public BranchGroup
 {

  public :

    // ctor's and d'tor
    CelLink() ; 
    CelLink( const CelLink & ) ;
    ~CelLink() ;

    // Operators
    CelLink & operator=( const CelLink & ) ;

    // Access
    Long64_t eventIndex() const { return _eventIndex ; }
    Long64_t metaIndex() const { return _metaIndex ; }
    Long64_t metaOffset() const { return _metaOffset ; }

    // print the info about the current event
    void printEventInfo( const char * options ) const ;  

  private :
  
    DataHandle<Long64_t> _eventIndex ; //! Absolute index of current composite event
    DataHandle<Long64_t> _metaIndex ;  //! Index of the associated entry in the File Tree
    DataHandle<Long64_t> _metaOffset ; //! Offset in events for this entry in the File Tree

    ClassDef(CelLink,0)
 } ;
 
#endif

