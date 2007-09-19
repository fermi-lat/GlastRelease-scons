
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
// a given entry in the file-tree-names tree, and also store the
// corresponding offset. The branches are :
// 
// Link_EventIndex  -> Absolute index of current composite event
// Link_SetIndex   -> Index of the associated entry in the file-tree-names Tree
// Link_SetOffset  -> Offset in events for this entry in the file-tree-names Tree
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
    
    // Accessors
    Long64_t eventIndex() const { return _eventIndex ; }
    Long64_t setIndex() const { return _setIndex ; }
    Long64_t setOffset() const { return _setOffset ; }
    void setEventIndex( Long64_t index )
     { _eventIndex = index ; }
    void incrementEventIndex()
     { _eventIndex = _eventIndex+1 ; }
    void incrementSetIndex()
     { _setOffset = _eventIndex ; _setIndex= _setIndex+1 ; }

    // print the info about the current event
    void printEventInfo( const char * options ) const ;  

  private :
  
    DataHandle<Long64_t> _eventIndex ; //! Absolute index of current composite event
    DataHandle<Long64_t> _setIndex ;  //! Index of the associated entry in the File Tree
    DataHandle<Long64_t> _setOffset ; //! Offset in events for this entry in the File Tree

    ClassDef(CelLink,0)
 } ;
 
#endif

