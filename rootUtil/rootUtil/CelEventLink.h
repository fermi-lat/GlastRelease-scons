
#ifndef CelEventLink_h
#define CelEventLink_h

/*
* Project: GLAST
* Package: rootUtil
*    File: $Id$
* Authors:
*   EC, Eric Charles, SLAC, echarles@slac.stanford.edu
*   DC, David Chamont, LLR, chamont@llr.in2p3.fr
*
* Copyright (c) 2007
*   Regents of Stanford University. All rights reserved.
*
*/

#include "BranchGroup.h"
#include "BgDataHandle.h"
#include <TObject.h>

//
// CelEventLink associates a given composite event with
// a given entry in the tree of file names, and also store the
// corresponding offset. The branches are :
// 
// Event_Index  -> Absolute index of current composite event
// Event_SetIndex   -> Index of the associated entry in the file-tree-names Tree
// Event_SetOffset  -> Offset in events for this entry in the file-tree-names Tree
// 

class CelEventLink : public BranchGroup
 {

  public :

    // ctor's and d'tor
    CelEventLink() ; 
    CelEventLink( const CelEventLink & ) ;
    ~CelEventLink() ;

    // Operators
    CelEventLink & operator=( const CelEventLink & ) ;
    
    // Accessors
    Long64_t eventIndex() const { return _eventIndex ; }
    Long64_t fileSetIndex() const { return _fileSetIndex ; }
    Long64_t fileSetOffset() const { return _fileSetOffset ; }
    void setEventIndex( Long64_t index )
     { _eventIndex = index ; }
    void incrementEventIndex()
     { _eventIndex = _eventIndex+1 ; }
    void incrementFileSetIndex()
     { _fileSetOffset = _eventIndex ; _fileSetIndex= _fileSetIndex+1 ; }

    // print the info about the current event
    void printInfo( const char * options ) const ;  

  private :
  
	BgDataHandle<Long64_t> _eventIndex ; //! Absolute index of current composite event
    BgDataHandle<Long64_t> _fileSetIndex ;   //! Index of the associated entry in the File Tree
    BgDataHandle<Long64_t> _fileSetOffset ;  //! Offset in events for this entry in the File Tree

    ClassDef(CelEventLink,0)
 } ;
 
#endif

