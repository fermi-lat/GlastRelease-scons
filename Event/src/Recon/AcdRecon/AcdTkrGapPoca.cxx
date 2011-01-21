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
  AcdTkrGapPoca::AcdTkrGapPoca(const idents::AcdGapId& gapId, int trackIndex, 
                               const float active[2],
                               float vetoSigmaHit, float vetoSigmaProj, float vetoSigmaProp,
                               int volumePlane, float arcLengthToPlane, float cosTheta, 
                               const HepPoint3D& global, const float localPosition[2], 
                               const HepSymMatrix& localCovProj, const HepSymMatrix& localCovProp,
                               int volume, int region, float arcLength, 
                               float doca, float docaErrProj, float docaErrProp,
                               const Point& poca, const Vector& voca)
    : AcdTkrLocalCoords(volumePlane,arcLengthToPlane,cosTheta,
                        global,localPosition,active,
                        localCovProj,localCovProp),
      AcdPocaData(volume,region,arcLength,doca,docaErrProj,docaErrProp,poca,voca),
      m_id(gapId),
      m_trackIndex(trackIndex),
      m_vetoSigmaHit(vetoSigmaHit),
      m_vetoSigmaProj(vetoSigmaProj),
      m_vetoSigmaProp(vetoSigmaProp){
  }
  
  /// Old Constructor for backwards compatiblity  
  AcdTkrGapPoca::AcdTkrGapPoca(const idents::AcdGapId& gapId, int trackIndex, 
                               const Event::AcdTkrLocalCoords& local, const Event::AcdPocaData& pocaData )
    : AcdTkrLocalCoords(local),
      AcdPocaData(pocaData),
      m_id(gapId),
      m_trackIndex(trackIndex),
      m_vetoSigmaHit(-1),
      m_vetoSigmaProj(-1),
      m_vetoSigmaProp(-1){
  }
  

  /// Copy constructor
  AcdTkrGapPoca::AcdTkrGapPoca(const Event::AcdTkrGapPoca& other)
    :AcdTkrLocalCoords(other),AcdPocaData(other)
  {
    set(other.getId(),other.trackIndex(),
        other.vetoSigmaHit(),other.vetoSigmaProj(),other.vetoSigmaProp());
  }

  /// Assignment operator
  AcdTkrGapPoca& AcdTkrGapPoca::operator=(const Event::AcdTkrGapPoca& other)
  {
    if ( this == &other ) return *this;
    set(other.getId(),other.trackIndex(),
        other.vetoSigmaHit(),other.vetoSigmaProj(),other.vetoSigmaProp());
    AcdTkrLocalCoords::copy(other);
    AcdPocaData::setPocaData(other);
    return *this;
  }

  bool AcdTkrGapPoca::operator<(const Event::AcdTkrGapPoca& other) const {
    if ( this == &other ) return false;
    float vs2 = vetoSigma2();
    float ovs2 = other.vetoSigma2();      
    if ( vs2 < ovs2 ) return true;
    if ( vs2 > ovs2 ) return false;

    float d = getDoca();
    float od = other.getDoca();
    if ( d < od ) return true;
    if ( d > od ) return false;

    if ( m_id.asShort() < other.getId().asShort() ) return true;
    if ( m_id.asShort() > other.getId().asShort() ) return false;
    
    return false;    
  }


  /// reset all the values to their default
  void AcdTkrGapPoca::ini()
  {
    m_id = idents::AcdGapId();
    m_trackIndex = -1;
    m_vetoSigmaHit = 0.;
    m_vetoSigmaProj = 0.;
    m_vetoSigmaProp = 0.;
    AcdTkrLocalCoords::ini();
    AcdPocaData::ini();
  }
  
  /// Print out this structure
  void AcdTkrGapPoca::writeOut(MsgStream& stream) const 
  {
    stream << MSG::DEBUG
           << "AcdTrkGapPoca.  Gap: " << (int)m_id.asShort()
           << ".  Track: " << (int)m_trackIndex << std::endl
           << "Sigma " << m_vetoSigmaHit << ' ' << m_vetoSigmaProj << ' ' << m_vetoSigmaProp << std::endl;
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
