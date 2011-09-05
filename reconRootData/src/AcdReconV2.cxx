#include "reconRootData/AcdReconV2.h"
#include "reconRootData/AcdHit.h"
#include "reconRootData/AcdTkrAssoc.h"
#include "reconRootData/AcdCalAssoc.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"
#include "TROOT.h"
#include "TClass.h"
#include "TStreamerInfo.h"

ClassImp(AcdReconV2)



AcdReconV2::AcdReconV2()
  :TObject(),
   m_acdHitCol( new TClonesArray( AcdHit::Class() ) ),
   m_acdTkrAssocCol( new TClonesArray( AcdTkrAssoc::Class() ) ),
   m_acdCalAssocCol( new TClonesArray( AcdCalAssoc::Class() ) ){  
  Clear();
}

AcdReconV2::~AcdReconV2() {
  if ( m_acdHitCol != 0 ) {
    m_acdHitCol->Delete();
  }
  delete m_acdHitCol;
  m_acdHitCol = 0;
  if ( m_acdTkrAssocCol != 0 ) {
    m_acdTkrAssocCol->Delete();
  }
  delete m_acdTkrAssocCol;
  m_acdTkrAssocCol = 0;
  if ( m_acdCalAssocCol != 0 ) {
    m_acdCalAssocCol->Delete();
  }
  delete m_acdCalAssocCol;
  m_acdCalAssocCol = 0;
  Clear();
}


void AcdReconV2::Clear(Option_t* option) {
  m_topology.Clear(option);
  if ( m_acdHitCol ) {
    m_acdHitCol->Clear();
  }
  if ( m_acdTkrAssocCol ) {
    m_acdTkrAssocCol->Clear();
  }
  if ( m_acdCalAssocCol ) {
    m_acdCalAssocCol->Clear();
  }
}

void AcdReconV2::Print(Option_t *option) const {
  TObject::Print(option);
  m_topology.Print(option);
}

AcdHit* AcdReconV2::addAcdHit(const AcdHit& toAdd) {
  UInt_t n = getHitCol().GetEntriesFast();
  TClonesArray& hits = *m_acdHitCol;
  AcdHit* newHit = (AcdHit*)(new (hits[n]) AcdHit(toAdd));
  return newHit;
}

AcdTkrAssoc* AcdReconV2::addAcdTkrAssoc(const AcdTkrAssoc& toAdd) {
  UInt_t n = getTkrAssocCol().GetEntriesFast();
  TClonesArray& assocs = *m_acdTkrAssocCol; 
  AcdTkrAssoc* newAssoc =  (AcdTkrAssoc*)(new (assocs[n]) AcdTkrAssoc(toAdd));
  return newAssoc;
}

AcdCalAssoc* AcdReconV2::addAcdCalAssoc(const AcdCalAssoc& toAdd) {
  UInt_t n = getCalAssocCol().GetEntriesFast();
  TClonesArray& assocs = *m_acdCalAssocCol; 
  AcdCalAssoc* newAssoc =  (AcdCalAssoc*)(new (assocs[n]) AcdCalAssoc(toAdd));
  return newAssoc;
}

//======================================================
// For Unit Tests
//======================================================

void AcdReconV2::Fake( Int_t /* ievent */, Float_t /* randNum */ ) {
  
}

#define UTIL_COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),a.get ## att(),#att)
#define ID_COMPARE_IN_RANGE(att) get ## att().CompareInRange(a.get ## att(),#att)

Bool_t AcdReconV2::CompareInRange( const AcdReconV2 & /* a */, const std::string & /* name */ ) const {

  bool result = true ;
  return result ;
}
