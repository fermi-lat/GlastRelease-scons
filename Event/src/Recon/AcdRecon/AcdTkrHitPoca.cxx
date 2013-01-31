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
                               const float active2d[2], const float mips[2], 
                               float vetoSigmaHit, float vetoSigmaProj, float vetoSigmaProp,
                               int volumePlane, float arcLengthToPlane, float cosTheta, 
                               const HepPoint3D& global, const float localPosition[2], 
                               const CLHEP::HepSymMatrix& localCovProj, const CLHEP::HepSymMatrix& localCovProp,
                               int volume, int region, float arcLength, 
                               float doca, float docaErrProj, float docaErrProp,
                               const Point& poca, const Vector& voca, const unsigned short flags[2])
    :AcdTkrLocalCoords(volumePlane,arcLengthToPlane,cosTheta,
                       global,localPosition,active2d,
                       localCovProj,localCovProp),
     AcdPocaData(volume,region,arcLength,doca,docaErrProj,docaErrProp,poca,voca),
     m_id(acdId),
     m_trackIndex(trackIndex),
     m_vetoSigmaHit(vetoSigmaHit),
     m_vetoSigmaProj(vetoSigmaProj),
     m_vetoSigmaProp(vetoSigmaProp){
    m_mips[0] = mips[0];
    m_mips[1] = mips[1];
    m_flags[0] = flags[0];
    m_flags[1] = flags[1];
  }
  
  AcdTkrHitPoca::AcdTkrHitPoca( const idents::AcdId& acdId, int trackIndex, 
                                const Event::AcdTkrLocalCoords& local, const Event::AcdPocaData& pocaData ) 
    :AcdTkrLocalCoords(local),
     AcdPocaData(pocaData),
     m_id(acdId),
     m_trackIndex(trackIndex),     
     m_vetoSigmaHit(-1.),
     m_vetoSigmaProj(-1.),
     m_vetoSigmaProp(-1.){
    m_mips[0] = -1.;
    m_mips[1] = -1.;
    m_flags[0] = 0;
    m_flags[1] = 0;
  }


  

  /// Copy constructor
  AcdTkrHitPoca::AcdTkrHitPoca(const Event::AcdTkrHitPoca& other)
    :AcdTkrLocalCoords(other),AcdPocaData(other)
  {
    set(other.getId(),other.trackIndex(),other.m_mips,
        other.vetoSigmaHit(),other.vetoSigmaProj(),other.vetoSigmaProp(),other.m_flags);
  }

  /// Assignment operator
  AcdTkrHitPoca& AcdTkrHitPoca::operator=(const Event::AcdTkrHitPoca& other)
  {
    if ( this == &other ) return *this;
    set(other.getId(),other.trackIndex(),other.m_mips,
        other.vetoSigmaHit(),other.vetoSigmaProj(),other.vetoSigmaProp(),other.m_flags);
    AcdTkrLocalCoords::copy(other);
    AcdPocaData::setPocaData(other);
    return *this;
  }


  /// Comparison operator
  bool AcdTkrHitPoca::operator<(const Event::AcdTkrHitPoca& other) const {
    // Check identity
    if ( this == &other ) return false;
    // First compare vetoSigma2
    float myVS2 = vetoSigma2();
    float otherVS2 = other.vetoSigma2();
    if ( myVS2 < otherVS2 ) return true;
    if ( myVS2 > otherVS2 ) return false;
    // They might be equal.  
    // This is the case when a track extapolate into two tiles with MIP-like hits.
    // myVS2 == otherVS2 == 0. 
    //   or two tiles without hits
    // myVS2 == otherVS2 == 1.e4
    // Larger hit wins.  
    if ( vetoSigmaHit() < other.vetoSigmaHit()) return true;
    if ( vetoSigmaHit() > other.vetoSigmaHit()) return false;
    // Still equal.  This is the case when a track extrapolates inside two tiles without hits    
    // Larger scaled active distance wins
    if ( vetoSigmaProj() < other.vetoSigmaProj()) return true;
    if ( vetoSigmaProj() > other.vetoSigmaProj()) return false;    
    // Try the acd id
    if ( getId().id() < other.getId().id() ) return true;
    if ( getId().id() > other.getId().id() ) return false;
    // use the pointer address.  
    return this < &other;    
  }


  /// set all the values
  void AcdTkrHitPoca::set(const idents::AcdId& acdId, int trackIndex,
                          const float mips[2], 
                          float vetoSigmaHit, float vetoSigmaProj, float vetoSigmaProp,
                          const unsigned short flags[2])
  {
    m_id = acdId;
    m_trackIndex = trackIndex;
    m_mips[0] = mips[0];
    m_mips[1] = mips[1];
    m_vetoSigmaHit = vetoSigmaHit;
    m_vetoSigmaProj = vetoSigmaProj;
    m_vetoSigmaProp = vetoSigmaProp;
    m_flags[0] = flags[0];
    m_flags[1] = flags[1];
  }  
  

  /// combine the sigma from the hit with the sigma from the track
  float AcdTkrHitPoca::vetoSigma2() const
  {
    // Only add positive values.  
    float retVal = m_vetoSigmaHit >= 0. ? (m_vetoSigmaHit*m_vetoSigmaHit) : 0.;
    retVal += m_vetoSigmaProj >= 0. ? (m_vetoSigmaProj*m_vetoSigmaProj) : 0.;
    return retVal;
  }

  
  /// reset all the values to their default
  void AcdTkrHitPoca::ini()
  {
    m_id = idents::AcdId();
    m_trackIndex = -1;
    m_mips[0] = 0.;
    m_mips[1] = 0.;
    m_vetoSigmaHit = 0.;
    m_vetoSigmaProj = 0.;
    m_vetoSigmaProp = 0.;
    m_flags[0] = 0;
    m_flags[1] = 0;
    AcdTkrLocalCoords::ini();
    AcdPocaData::ini();
  }
  
  /// Print out this structure
  void AcdTkrHitPoca::writeOut(MsgStream& stream) const 
  {
    stream << MSG::DEBUG
           << "AcdTkrHitPoca.  Tile: " << m_id.id()  
           << ".  Track: " << ((int)m_trackIndex) 
           << ".  Mips:  " << m_mips[0] << ',' << m_mips[1] << std::endl
           << "Sigma " << m_vetoSigmaHit << ' ' << m_vetoSigmaProj << ' ' << m_vetoSigmaProp << std::endl;
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
