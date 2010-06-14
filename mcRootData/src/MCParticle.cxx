///////////////////////////////////////////////////////////////////////////
//                                                                       
// The MCParticle class provides an interface to a particle for use in
// GLAST Monte Carlo applications
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "TClass.h"

#include "mcRootData/MCParticle.h"
#include "mcRootData/MCVertex.h"

ClassImp(MCParticle)

/*! Default Constructor
 *  Reserved for use by ROOT I/O.  Do not use this constructor, as no
 *  initialization is performed.
 */
MCParticle::MCParticle() {
}
//________________________________________________________________________
MCParticle::MCParticle(TLorentzVector *pMom, Int_t id, MCParticle *pParent, MCVertex *pVertex) {
  m_pMom = pMom;
  m_nPDGId = id;
  m_pParent = pParent;
  m_pVert = pVertex;
}
//________________________________________________________________________
MCParticle::MCParticle(Double_t momX, Double_t momY, Double_t momZ, Double_t energy,
                       Int_t id, MCParticle *pParent, MCVertex *pVertex) {
  m_pMom = new TLorentzVector(momX, momY, momZ, energy);
  m_nPDGId = id;
  m_pParent = pParent;
  m_pVert = pVertex;
}
//________________________________________________________________________
MCParticle::~MCParticle() {
  // destructor
  if (m_pMom)
    delete m_pMom;
  if (m_pVert)
    delete m_pVert;
}
//________________________________________________________________________
Double_t MCParticle::getMomX() {
  return m_pMom->X();
}
//________________________________________________________________________
Double_t MCParticle::getMomY() {
  return m_pMom->Y();
}
//________________________________________________________________________
Double_t MCParticle::getMomZ() {
  return m_pMom->Z();
}
//________________________________________________________________________
Double_t MCParticle::getMom() {
  Double_t  x = m_pMom->X(),
            y = m_pMom->Y(),
            z = m_pMom->Z();

  return (Double_t)sqrt(x*x + y*y + z*z);
}
//________________________________________________________________________
Double_t MCParticle::getEnergy() {
  return m_pMom->T();
}
//________________________________________________________________________
TLorentzVector *MCParticle::getMomLV() {
  return m_pMom;
}
//________________________________________________________________________
MCVertex *MCParticle::getVertex() {
  return m_pVert;
}
//________________________________________________________________________
MCParticle *MCParticle::getParent() {
  return m_pParent;
}
//________________________________________________________________________
Int_t MCParticle::getPDGId() {
  return m_nPDGId;
}
//________________________________________________________________________
void MCParticle::setMomX(Double_t momX) {
  m_pMom->SetX(momX);
}
//________________________________________________________________________
void MCParticle::setMomY(Double_t momY) {
  m_pMom->SetY(momY);
}
//________________________________________________________________________
void MCParticle::setMomZ(Double_t momZ) {
  m_pMom->SetZ(momZ);
}
//________________________________________________________________________
void MCParticle::setEnergy(Double_t energy) {
  m_pMom->SetT(energy);
}
//________________________________________________________________________
