
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

//
// CelEventLink associates the current composite event with
// a given set of file names, and also store the
// corresponding offset.
// 

class CelEventLink : public BranchGroup
 {

  public :

    // ctor's and d'tor
    CelEventLink() ; 
    CelEventLink( const CelEventLink & ) ;
    ~CelEventLink() ;
    CelEventLink & operator=( const CelEventLink & ) ;
    
    // Accessors
    Long64_t eventIndex() const { return _eventIndex ; }
    Long64_t setIndex() const { return _setIndex ; }
    void changeEventIndex( Long64_t index ) { _eventIndex = index ; }
    void incrementEventIndex() { _eventIndex = _eventIndex+1 ; }
    void incrementFileSetIndex() { _setIndex = _setIndex+1 ; }

    // print the info about the current event
    void printInfo( const char* options ) const ;

  private :
  
	BgDataHandle<Long64_t> _eventIndex ;     //! Absolute index of current composite event
    BgDataHandle<Long64_t> _setIndex ;   //! Index of the associated entry in the File Tree

    ClassDef(CelEventLink,0)
 } ;
 
#endif

