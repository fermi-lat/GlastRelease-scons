///////////////////////////////////////////////////////////////////////////
//                                                                       
// The MCVertex class provides a 4-Space (x,y,z,t) position for the
// MCParticle class to use
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "TClass.h"

#include "mcRootData/MCVertex.h"

ClassImp(MCVertex)

/*! Default Constructor
 *  Reserved for use by ROOT I/O.  Do not use this constructor, as no
 *  initialization is performed.
 */
MCVertex::MCVertex() {
}
//________________________________________________________________________
MCVertex::MCVertex(TLorentzVector *pos) {
  m_pPos = pos;
}
//________________________________________________________________________
MCVertex::MCVertex(Double_t posX, Double_t posY, Double_t posZ, Double_t time) {
  m_pPos = new TLorentzVector(posX, posY, posZ, time);
}
//________________________________________________________________________
MCVertex::~MCVertex() {
  if (m_pPos)
    delete m_pPos;
}
//________________________________________________________________________
Double_t MCVertex::getPosX() {
  return m_pPos->X();
}
//________________________________________________________________________
Double_t MCVertex::getPosY() {
  return m_pPos->Y();
}
//________________________________________________________________________
Double_t MCVertex::getPosZ() {
  return m_pPos->Z();
}
//________________________________________________________________________
Double_t MCVertex::getTime() {
  return m_pPos->T();
}
//________________________________________________________________________
TLorentzVector *MCVertex::getPosLV() {
  return m_pPos;
}
//________________________________________________________________________
void MCVertex::setPosX(Double_t posX) {
  m_pPos->SetX(posX);
}
//________________________________________________________________________
void MCVertex::setPosY(Double_t posY) {
  m_pPos->SetY(posY);
}
//________________________________________________________________________
void MCVertex::setPosZ(Double_t posZ) {
  m_pPos->SetZ(posZ);
}
//________________________________________________________________________
void MCVertex::setTime(Double_t time) {
  m_pPos->SetT(time);
}
//________________________________________________________________________