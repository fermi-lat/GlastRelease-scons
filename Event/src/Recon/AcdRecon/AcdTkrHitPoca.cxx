#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"

#include "GaudiKernel/MsgStream.h"
#include "CLHEP/Matrix/Matrix.h"

namespace Event {

  /// Default constructor.  Set everything to default values
  AcdTkrHitPoca::AcdTkrHitPoca() 
    :AcdTkrLocalCoords(),AcdPocaData(){
    ini();
  }

  /// Constructor for use in transient -> persistent conversion 
  /// Takes arguements as they are stored in ROOT
  AcdTkrHitPoca::AcdTkrHitPoca(const idents::AcdId& acdId, int trackIndex,
			       const Event::AcdTkrLocalCoords& local,
			       const Event::AcdPocaData& pocaData)
    :AcdTkrLocalCoords(),AcdPocaData()
  {
    set(acdId,trackIndex,local,pocaData);
  }
  
  /// Copy constructor
  AcdTkrHitPoca::AcdTkrHitPoca(const Event::AcdTkrHitPoca& other)
    :AcdTkrLocalCoords(),AcdPocaData()
  {
    set(other.getId(),other.trackIndex(),other,other);
  }

  /// Assignment operator
  AcdTkrHitPoca& AcdTkrHitPoca::operator=(const Event::AcdTkrHitPoca& other)
  {
    if ( this == &other ) return *this;
    set(other.getId(),other.trackIndex(),other,other);
    return *this;
  }


  /// set all the values
  void AcdTkrHitPoca::set(const idents::AcdId& acdId, int trackIndex,
			  const Event::AcdTkrLocalCoords& local,
			  const Event::AcdPocaData& pocaData)
  {
    m_id = acdId;
    m_trackIndex = trackIndex;
    AcdTkrLocalCoords::copy(local);
    AcdPocaData::set(pocaData);
  }  
  
  
  /// reset all the values to their default
  void AcdTkrHitPoca::ini()
  {
    m_id = idents::AcdId();
    m_trackIndex = -1;
    AcdTkrLocalCoords::ini();
    AcdPocaData::ini();
  }
  
  /// Print out this structure
  void AcdTkrHitPoca::writeOut(MsgStream& stream) const 
  {
    stream << MSG::DEBUG
	   << "AcdTkrHitPoca.  Tile: " << m_id.id() 
	   << ".  Track: " << (int)m_trackIndex
	   << ".  ";
    AcdTkrLocalCoords::writeOut(stream);
    AcdPocaData::writeOut(stream);
    stream << endreq;
  }

  /// Copy c'tor
  AcdTkrHitPocaCol::AcdTkrHitPocaCol(const std::vector<AcdTkrHitPoca*>& other) {
    for ( std::vector<AcdTkrHitPoca*>::const_iterator itr = other.begin();
	  itr != other.end(); itr++ ) {
      AcdTkrHitPoca* onePoca = const_cast<AcdTkrHitPoca*>(*itr);
      add(onePoca);
    }
  }
  
  void AcdTkrHitPocaCol::delTkrHitPocas() 
    //Purpose: delete all AcdTkrIntersection object from memory    
  {
    int n = num();
    for (int i = 0; i < n; i++) {
      delete operator[](i);
    }
    clear();
  }
  
  void AcdTkrHitPocaCol::ini()
    
    //Purpose:  delete all pointers to clusters
    // from collection 
    
  {
    clear();
  }
  
  
  void AcdTkrHitPocaCol::writeOut(MsgStream& stream) const
    
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
