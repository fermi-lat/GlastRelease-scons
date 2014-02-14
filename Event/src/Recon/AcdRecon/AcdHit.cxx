// File and Version information:
// $Header$
//
//  Implementation file of AcdHit and AcdHitCol classes
//  
// Authors:
//
//    Eric Charles
//
//

#include "Event/Recon/AcdRecon/AcdHit.h"
#include "Event/Digi/AcdDigi.h"

#include "GaudiKernel/MsgStream.h"

using namespace Event;

/// Constructor for use in reconstruction, takes digi and calibrated values
AcdHit::AcdHit(const Event::AcdDigi& digi, float mipsPmtA, float mipsPmtB) {
  // null values
  ini();
  // get the raw PHA values
  m_pha[A] = digi.getPulseHeight(AcdDigi::A);
  m_pha[B] = digi.getPulseHeight(AcdDigi::B);
  
  // get the calibrated values
  m_mipsPmt[A] = mipsPmtA;
  m_mipsPmt[B] = mipsPmtB;

  // get the flags from the digi
  setFlags(digi);
}

/// Constructor for use in persistent -> transient conversion, 
/// Takes arguements as they are stored in ROOT
AcdHit::AcdHit(const idents::AcdId& id, unsigned short flagsA, unsigned short flagsB, 
               unsigned short phaA, unsigned short phaB,
               float mipsPmtA, float mipsPmtB) {
  set(id,flagsA,flagsB,phaA,phaB,mipsPmtA,mipsPmtB);
}


void AcdHit::set(const idents::AcdId& id, 
                 unsigned short flagsA, unsigned short flagsB, 
                 unsigned short phaA, unsigned short phaB,
                 float mipsPmtA, float mipsPmtB) {
  // just copy everything
  m_acdId = id;
  m_flags[A] = flagsA;
  m_flags[B] = flagsB;
  m_pha[A] = phaA;
  m_pha[B] = phaB;
  m_mipsPmt[A] = mipsPmtA;
  m_mipsPmt[B] = mipsPmtB;
}

void AcdHit::writeOut(MsgStream& stream) const
// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  stream << MSG::DEBUG
         << "Tile: " << m_acdId.id()
         << "Phas: (" << m_pha[A] << ',' << m_pha[B] 
         << "). Flags: (" << m_flags[A] << ',' << m_flags[B] 
         << "). Mips: (" << m_mipsPmt[A] << ',' <<  m_mipsPmt[B]
         << ")."
         << endreq;
}



void AcdHit::ini()
// Purpose: reset all data members to 0
//
{ 
  m_acdId = idents::AcdId();
  m_flags[A] = m_flags[B] = 0;
  m_pha[A] = m_pha[B] = 0;
  m_mipsPmt[A] = m_mipsPmt[B] = 0.;
}


void AcdHit::setFlags(const Event::AcdDigi& digi)
// Purpose: to pull various flags out of the digis into the hits
//
{
  m_acdId = digi.getId();
  m_flags[A] = 0; m_flags[B] = 0;

  // FIXME  -> Do all the flags
  m_flags[A] |= digi.getAcceptMapBit(Event::AcdDigi::A) << PMT_ACCEPT_BIT;
  m_flags[B] |= digi.getAcceptMapBit(Event::AcdDigi::B) << PMT_ACCEPT_BIT;

  m_flags[A] |= digi.getHitMapBit(Event::AcdDigi::A) << PMT_VETO_BIT;
  m_flags[B] |= digi.getHitMapBit(Event::AcdDigi::B) << PMT_VETO_BIT;

  m_flags[A] |= (digi.getRange(Event::AcdDigi::A)==Event::AcdDigi::HIGH) << PMT_RANGE_BIT;
  m_flags[B] |= (digi.getRange(Event::AcdDigi::B)==Event::AcdDigi::HIGH) << PMT_RANGE_BIT;

  m_flags[A] |= digi.getCno(Event::AcdDigi::A) << PMT_CNO_BIT;
  m_flags[B] |= digi.getCno(Event::AcdDigi::B) << PMT_CNO_BIT;

  m_flags[A] |= (digi.getOddParityError(Event::AcdDigi::A)==Event::AcdDigi::ERROR) << PMT_ODD_PARITY_ERROR_BIT;
  m_flags[B] |= (digi.getOddParityError(Event::AcdDigi::B)==Event::AcdDigi::ERROR) << PMT_ODD_PARITY_ERROR_BIT;

  m_flags[A] |= (digi.getHeaderParityError(Event::AcdDigi::A)==Event::AcdDigi::ERROR) << PMT_HEADER_PARITY_ERROR_BIT;
  m_flags[B] |= (digi.getHeaderParityError(Event::AcdDigi::B)==Event::AcdDigi::ERROR) << PMT_HEADER_PARITY_ERROR_BIT;
  
}


void AcdHit::correctAcceptMapBits(bool acceptA, bool acceptB) 
/// this is to allow us to kill accept map bits for periodic triggers in the overlays
{
  // This masks out the accept bit
  static unsigned short clearAcceptMask(0xFFFE);
  m_flags[A] &= clearAcceptMask;
  m_flags[B] &= clearAcceptMask;
  
  if (acceptA)  m_flags[A] |= PMT_ACCEPT_MASK;
  if (acceptB)  m_flags[B] |= PMT_ACCEPT_MASK;
}


AcdHitCol::AcdHitCol(const std::vector<AcdHit*>& acdhits) {
//Purpose: take ownership of hits from a vector
  for ( std::vector<AcdHit*>::const_iterator itr = acdhits.begin();
        itr != acdhits.end(); itr++ ) {
    AcdHit* hit = const_cast<AcdHit*>(*itr);
    add(hit);
  }
}

void AcdHitCol::delHits()
//Purpose: delete all AcdHit object from memory
{
  int nHit = num();
  for (int iHit = 0; iHit < nHit; iHit++) {
    delete operator[](iHit);
  }
  clear();
}

void AcdHitCol::ini()
//Purpose:  delete all pointers to clusters
// from collection 
{
  clear();
}


void AcdHitCol::writeOut(MsgStream& stream) const
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
