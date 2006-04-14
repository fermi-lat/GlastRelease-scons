#include "Event/Recon/AcdRecon/AcdTkrGapPoca.h"
#include <iostream>

#include "GaudiKernel/MsgStream.h"

namespace Event {

  /// Default constructor.  Set everything to default values
  AcdTkrGapPoca::AcdTkrGapPoca() 
    :AcdTkrLocalCoords(),AcdPocaData(){
    ini();
  }

  /// Constructor for use in transient -> persistent conversion 
  /// Takes arguements as they are stored in ROOT
  AcdTkrGapPoca::AcdTkrGapPoca(const idents::AcdGapId& acdId, int trackIndex,
			       const Event::AcdTkrLocalCoords& local,
			       const Event::AcdPocaData& pocaData)
    :AcdTkrLocalCoords(),AcdPocaData()
  {
    set(acdId,trackIndex,local,pocaData);
  }
  
  /// Copy constructor
  AcdTkrGapPoca::AcdTkrGapPoca(const Event::AcdTkrGapPoca& other)
    :AcdTkrLocalCoords(),AcdPocaData()
  {
    set(other.getId(),other.trackIndex(),other,other);
  }

  /// Assignment operator
  AcdTkrGapPoca& AcdTkrGapPoca::operator=(const Event::AcdTkrGapPoca& other)
  {
    if ( this == &other ) return *this;
    set(other.getId(),other.trackIndex(),other,other);
    return *this;
  }


  /// set all the values
  void AcdTkrGapPoca::set(const idents::AcdGapId& acdId, int trackIndex,
			  const Event::AcdTkrLocalCoords& local,
			  const Event::AcdPocaData& pocaData)
  {
    m_id = acdId;
    m_trackIndex = trackIndex;
    AcdTkrLocalCoords::copy(local);
    AcdPocaData::set(pocaData);
  }  
  
  
  /// reset all the values to their default
  void AcdTkrGapPoca::ini()
  {
    m_id = idents::AcdGapId();
    m_trackIndex = -1;
    AcdTkrLocalCoords::ini();
    AcdPocaData::ini();
  }
  
  /// Print out this structure
  void AcdTkrGapPoca::writeOut(MsgStream& stream) const 
  {
    stream << MSG::DEBUG
	   << "AcdTrkGapPoca.  Gap: " << (int)m_id.asShort()
	   << ".  Track: " << (int)m_trackIndex
	   << ".  ";
    AcdTkrLocalCoords::writeOut(stream);
    AcdPocaData::writeOut(stream);
    stream << endreq;
  }

  /// Copy c'tor
  AcdTkrGapPocaCol::AcdTkrGapPocaCol(const std::vector<AcdTkrGapPoca*>& other) {
    for ( std::vector<AcdTkrGapPoca*>::const_iterator itr = other.begin();
	  itr != other.end(); itr++ ) {
      AcdTkrGapPoca* onePoca = const_cast<AcdTkrGapPoca*>(*itr);
      add(onePoca);
    }
  }
  
  void AcdTkrGapPocaCol::del() 
    //Purpose: delete all AcdTkrIntersection object from memory    
  {
    int n = num();
    for (int i = 0; i < n; i++) {
      delete operator[](i);
    }
    clear();
  }
  
  void AcdTkrGapPocaCol::ini()
    
    //Purpose:  delete all pointers to clusters
    // from collection 
    
  {
    clear();
  }
  
  
  void AcdTkrGapPocaCol::writeOut(MsgStream& stream) const
    
    // Purpose: provide symbolic output of some data members
    //          of all clusters in collection for debugging purposes
    //
    // Input:
    //        stream - Gaudi message stream
  {
    
    // loop over all clusters
    for (unsigned int i = 0; i < size();i++) {
      
      // call the writeOut() method for each cluster
      (operator[](i))->writeOut(stream);
    }
    
  }
}
