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
  
  AcdPocaData::AcdPocaData(float arcLength, float doca, float docaErr, 
			   const Point& poca, const Vector& pocaVector) {
    set(arcLength,doca,docaErr,poca,pocaVector);
  }
  
  AcdPocaData::AcdPocaData(const AcdPocaData& other) {
    set(other.m_arcLength,other.m_doca,other.m_docaErr,other.m_poca,other.m_pocaVector);
  }
  
  AcdPocaData& AcdPocaData::operator=(const AcdPocaData& other) {
    set(other);
    return *this;
  }

  void AcdPocaData::set(float arcLength, float doca, float docaErr, 
			const Point& poca, const Vector& pocaVector) {
    m_arcLength = arcLength;
    m_doca = doca;
    m_docaErr = docaErr;
    
    m_poca = poca;
    m_pocaVector = pocaVector;
  }
  
  void AcdPocaData::set(const Event::AcdPocaData& other) {
    set(other.m_arcLength,other.m_doca,other.m_docaErr,other.m_poca,other.m_pocaVector);
  }
  
  void AcdPocaData::writeOut(MsgStream& stream ) const
    // Purpose: provide ascii output of some data members for
    //          debugging purposes
    // Input:
    //        stream - Gaudi message stream
  {
    stream << "Arc: " << m_arcLength
	   << ".  Doca: " << m_doca << " +- " << m_docaErr
	   << ".  Poca: (" << m_poca.x() << ',' << m_poca.y() << ',' << m_poca.z()
	   << ".  PocaDir: (" << m_pocaVector.x() << ',' << m_pocaVector.y() << ',' << m_pocaVector.z()
	   << ").  ";
  }

  
  
  void AcdPocaData::ini()
    // Purpose: reset all data members to 0
    //
  {
    m_arcLength = 0.;
    m_doca = 0.;
    m_docaErr = 0.;
    
    m_poca = Point();
    m_pocaVector = Vector();
  }

}
