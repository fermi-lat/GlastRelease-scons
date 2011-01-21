// File and Version information:
// $Header$
//
//  Implementation file of AcdTkrPoca and AcdTkrPocaCol classes
//  
// Authors:
//
//    Eric Charles
//
//

#include "Event/Recon/AcdRecon/AcdPocaData.h"

#include "GaudiKernel/MsgStream.h"


namespace Event {

  AcdPocaData::AcdPocaData() {
    ini();
  }
  
  AcdPocaData::AcdPocaData(int volume, int region, float arcLength, 
                           float doca, float docaErrProj, float docaErrProp,
                           const Point& poca, const Vector& voca) {
    setPocaData(volume,region,arcLength,doca,docaErrProj,docaErrProp,poca,voca);
  }
  
  AcdPocaData::AcdPocaData(const AcdPocaData& other) {
    setPocaData(other.m_volume,other.m_region,other.m_arcLength,
                other.m_doca,other.m_docaErr_proj,other.m_docaErr_prop,
                other.m_poca,other.m_voca);
  }
  
  AcdPocaData& AcdPocaData::operator=(const AcdPocaData& other) {
    setPocaData(other);
    return *this;
  }

  void AcdPocaData::setPocaData(int volume, int region, float arcLength, 
                                float doca, float docaErrProj, float docaErrProp,
                                const Point& poca, const Vector& pocaVector) {

    m_volume = volume;
    m_region = region;
    m_arcLength = arcLength;
    m_doca = doca;
    m_docaErr_proj = docaErrProj;
    m_docaErr_prop = docaErrProp;
    m_poca = poca;
    m_voca = pocaVector;
  }

  void AcdPocaData::setPocaData(float arcLength, float doca, float docaErr, 
                                const Point& poca, const Vector& pocaVector) {
    m_volume = -1;
    m_region = -1;
    m_arcLength = arcLength;
    m_doca = doca;
    m_docaErr_proj = docaErr;
    m_docaErr_prop = docaErr;
    m_poca = poca;
    m_voca = pocaVector;
  }
  
  void AcdPocaData::setPocaData(const Event::AcdPocaData& other) {
    setPocaData(other.m_volume,other.m_region,other.m_arcLength,
                other.m_doca,other.m_docaErr_proj,other.m_docaErr_prop,
                other.m_poca,other.m_voca);
  }
  
  void AcdPocaData::writeOut(MsgStream& stream ) const
    // Purpose: provide ascii output of some data members for
    //          debugging purposes
    // Input:
    //        stream - Gaudi message stream
  {
    stream << "Vol:     " << m_volume << ':' << m_region 
           << ".  Arc: " << m_arcLength
           << ".  Doca: " << m_doca 
           << " +- " << m_docaErr_proj << "(proj)"
           << " +- " << m_docaErr_prop << "(prop)"
           << ".  Poca: (" << m_poca.x() << ',' << m_poca.y() << ',' << m_poca.z()
           << ".  PocaDir: (" << m_voca.x() << ',' << m_voca.y() << ',' << m_voca.z()
           << ").  ";
  }

  
  
  void AcdPocaData::ini()
    // Purpose: reset all data members to 0
    //
  {
    m_volume = -1;
    m_region = -1;
    m_arcLength = 0.;
    m_doca = 0.;
    m_docaErr_proj = 0.;
    m_docaErr_prop = 0.;
    m_poca = Point();
    m_voca = Vector();
  }

}
