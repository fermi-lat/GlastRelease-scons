#include "Event/Recon/AcdRecon/AcdSplashVars.h"

#include "GaudiKernel/MsgStream.h"
#include "CLHEP/Matrix/Matrix.h"

namespace Event {

  /// Default constructor.  Set everything to default values
  AcdSplashVars::AcdSplashVars() 
  {
    ini();
  }

  /// Constructor for use in transient -> persistent conversion 
  /// Takes arguements as they are stored in ROOT
  AcdSplashVars::AcdSplashVars(const idents::AcdId& acdId, int trackIndex,
			   const Point& calEntryPoint, const Vector& calEntryVector,
			   const float& tileSolidAngle, const float& weightedTrackAngle,
			   const float& weightedPathlength)
  {
    set(acdId,trackIndex,calEntryPoint,calEntryVector,tileSolidAngle,weightedTrackAngle,weightedPathlength);
  }
  
  /// Copy constructor
  AcdSplashVars::AcdSplashVars(const Event::AcdSplashVars& other)
  {
    set(other.getId(),other.getTrackIndex(),
	other.calEntryPoint(),other.calEntryVector(),
	other.tileSolidAngle(),other.weightedTrackAngle(),other.weightedPathlength());
  }

  /// Assignment operator
  AcdSplashVars& AcdSplashVars::operator=(const Event::AcdSplashVars& other)
  {
    if ( this == &other ) return *this;
    set(other.getId(),other.getTrackIndex(),
	other.calEntryPoint(),other.calEntryVector(),
	other.tileSolidAngle(),other.weightedTrackAngle(),other.weightedPathlength());
    return *this;
  }


  /// set all the values
  void AcdSplashVars::set(const idents::AcdId& acdId, int trackIndex, 
			const Point& calEntryPoint, const Vector& calEntryVector,
			const float& tileSolidAngle, const float& weightedTrackAngle,
			const float& weightedPathlength)
  {
    m_id = acdId;
    m_trackIndex = trackIndex;
    m_calEntryPoint = calEntryPoint;
    m_calEntryVector = calEntryVector;
    m_tileSolidAngle = tileSolidAngle;
    m_weightedTrackAngle = weightedTrackAngle;
    m_weightedPathlength = weightedPathlength;
  }  
  
  
  /// reset all the values to their default
  void AcdSplashVars::ini()
  {
    m_id = idents::AcdId();
    m_trackIndex = -1;
    m_calEntryPoint = Point();
    m_calEntryVector = Vector();
    m_tileSolidAngle = -1.;
    m_weightedTrackAngle = -99.;
    m_weightedPathlength = -1.;
  }
  
  /// Print out this structure
  void AcdSplashVars::writeOut(MsgStream& stream) const 
  {
    stream << MSG::DEBUG
	   << "AcdSplashVars.  Tile: " << m_id.id() 
	   << ".  Track: " << (int)m_trackIndex
	   << ".  ";
    stream << m_calEntryPoint << std::endl;
    stream << m_calEntryVector << std::endl;
    stream << "Solid Angle: " << m_tileSolidAngle << std::endl;
    stream << "Weigted Track Angle: " << m_weightedTrackAngle << std::endl;
    stream << "Weigted Path Length: " << m_weightedPathlength << std::endl;
    stream << endreq;
  }

  /// Copy c'tor
  AcdSplashVarsCol::AcdSplashVarsCol(const std::vector<AcdSplashVars*>& other) {
    for ( std::vector<AcdSplashVars*>::const_iterator itr = other.begin();
	  itr != other.end(); itr++ ) {
      AcdSplashVars* one = const_cast<AcdSplashVars*>(*itr);
      add(one);
    }
  }
  
  void AcdSplashVarsCol::del()
    //Purpose: delete all AcdTkrIntersection object from memory    
  {
    int n = num();
    for (int i = 0; i < n; i++) {
      delete operator[](i);
    }
    clear();
  }
  
  void AcdSplashVarsCol::ini()
    
    //Purpose:  delete all pointers to clusters
    // from collection 
    
  {
    clear();
  }
  
  
  void AcdSplashVarsCol::writeOut(MsgStream& stream) const
    
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
