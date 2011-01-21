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

#include "Event/Recon/AcdRecon/AcdTkrPoca.h"

#include "GaudiKernel/MsgStream.h"


using namespace Event;

AcdTkrPoca::AcdTkrPoca() {
  ini();
}

AcdTkrPoca::AcdTkrPoca(const idents::AcdId& acdId, int trackIndex,
                       double doca, double docaErr, unsigned docaRegion,
                       const Point& poca, const Event::TkrTrackParams& paramsAtPoca) {
  set(acdId,trackIndex,doca,docaErr,docaRegion,poca,paramsAtPoca);
}

void AcdTkrPoca::writeOut(MsgStream& /* stream */ ) const
// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

}

void AcdTkrPoca::set(const idents::AcdId& acdId, int trackIndex,
                     double doca, double docaErr, unsigned docaRegion,
                     const Point& poca, const Event::TkrTrackParams& paramsAtPoca)
// Purpose: set all data members at once
//
{
  m_acdId = acdId;
  m_trackIdx = trackIndex;
  m_doca = doca;
  m_docaErr = docaErr;
  m_docaRegion = docaRegion;
  
  m_poca = poca;
  m_paramsAtPoca = paramsAtPoca;
}


void AcdTkrPoca::ini()
// Purpose: reset all data members to 0
//
{
  m_acdId = idents::AcdId();
  m_trackIdx = 0;
  m_doca = 0.;
  m_docaErr = 0.;
  m_docaRegion = NONE_TILE;
  
  m_poca = Point();
  m_paramsAtPoca = TkrTrackParams();
}

AcdTkrPocaCol::AcdTkrPocaCol(const std::vector<AcdTkrPoca*>& acdTkrPocas) {
//Purpose: take ownership of pocas from another vector
  for ( std::vector<AcdTkrPoca*>::const_iterator itr = acdTkrPocas.begin();
        itr != acdTkrPocas.end(); itr++ ) {
    AcdTkrPoca* poca = const_cast<AcdTkrPoca*>(*itr);
    add(poca);
  }
}

void AcdTkrPocaCol::delTkrPocas() 

//Purpose: delete all AcdTkrPoca object from memory

{
  int nInter = num();
  for (int iIn = 0; iIn < nInter; iIn++) {
    delete operator[](iIn);
  }
  clear();
}

void AcdTkrPocaCol::ini()

//Purpose:  delete all pointers from collection 

{
  clear();
}


void AcdTkrPocaCol::writeOut(MsgStream& stream) const

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
