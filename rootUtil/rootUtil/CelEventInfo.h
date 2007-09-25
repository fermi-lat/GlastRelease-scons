
#ifndef CelEventInfo_h
#define CelEventInfo_h

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
class CelEventComponent ;
class TObjArray ;

//
// CelEventInfo is made for reading/writing entries of the
// events internal tree. It associates a given composite event with
// a given entry in the file-tree-names tree, and also store the
// corresponding offset. The ids are not yet used. The branches are :
// 
// Event_Index  -> Absolute index of current composite event
// Event_RunId/L    -> Unique? id of the corresponding run
// Event_EventId/L  -> Unique? id of the corresponding event
// Event_FileSetIndex   -> Index of the associated entry in the file-tree-names Tree
// Event_FileSetOffset  -> Offset in events for this entry in the file-tree-names Tree
// 

class CelEventInfo : public BranchGroup
 {

  public :

    // ctor's and d'tor
    CelEventInfo() ; 
    CelEventInfo( const CelEventInfo & ) ;
    ~CelEventInfo() ;

    // Operators
    CelEventInfo & operator=( const CelEventInfo & ) ;
    
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
    void printEventInfo( const char * options ) const ;  

  private :
  
	DataHandle<Long64_t> _eventIndex ; //! Absolute index of current composite event
	DataHandle<Long64_t> _runId ;      //! placeholder
	DataHandle<Long64_t> _eventId ;    //! placeholder
    DataHandle<Long64_t> _fileSetIndex ;   //! Index of the associated entry in the File Tree
    DataHandle<Long64_t> _fileSetOffset ;  //! Offset in events for this entry in the File Tree

    ClassDef(CelEventInfo,0)
 } ;
 
#endif

