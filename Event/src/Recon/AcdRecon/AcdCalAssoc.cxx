// File and Version information:
// $Header$
//
//  Implementation file of AcdCalAssoc and AcdCalAssocCol classes
//  
// Authors:
//
//    Alex Drlica-Wagner
//
//

#include "Event/Recon/AcdRecon/AcdCalAssoc.h"

#include "GaudiKernel/MsgStream.h"

using namespace Event;

AcdCalAssoc::AcdCalAssoc()
  :m_index(-1),
   m_upward(true),
   m_energy(0.),
   m_start(0.,0.,0.),
   m_dir(0.,0.,0.),
   m_arcLength(0.),
   m_cov_start(5,0),
   m_cov_end(5,0),
   m_tkrSSDVeto(0),
   m_cornerDoca(0){
}

/// Copy constructor
AcdCalAssoc::AcdCalAssoc(const AcdCalAssoc& other)
  :m_index(other.m_index),
   m_upward(other.m_upward),
   m_energy(other.m_energy),
   m_start(other.m_start),
   m_dir(other.m_dir),
   m_arcLength(other.m_arcLength),
   m_cov_start(other.m_cov_start),
   m_cov_end(other.m_cov_end),
   m_tkrSSDVeto(other.m_tkrSSDVeto),
   m_cornerDoca(other.m_cornerDoca){
}


/// Constructor for use in reconstruction, 
AcdCalAssoc::AcdCalAssoc(int index, bool up, float energy, 
                         const HepPoint3D& start, const HepVector3D& dir, float arcLength,
                         const CLHEP::HepSymMatrix& covStart, const CLHEP::HepSymMatrix& covEnd,
                         int tkrSSDVeto, float cornerDoca)
  :m_index(index),
   m_upward(up),
   m_energy(energy),
   m_start(start),
   m_dir(dir),
   m_arcLength(arcLength),
   m_cov_start(covStart),
   m_cov_end(covEnd),
   m_tkrSSDVeto(tkrSSDVeto),
   m_cornerDoca(cornerDoca){
}


void AcdCalAssoc::set(int index, bool up, float energy, 
                      const HepPoint3D& start, const HepVector3D& dir, float arcLength,
                      const CLHEP::HepSymMatrix& covStart, const CLHEP::HepSymMatrix& covEnd,
                      int tkrSSDVeto, float cornerDoca){
  // just in copy everything
  m_index = index;
  m_upward = up;
  m_energy = energy;
  m_start = start;
  m_dir = dir;
  m_arcLength = arcLength;
  m_cov_start = covStart;
  m_cov_end = covEnd;
  m_tkrSSDVeto = tkrSSDVeto;
  m_cornerDoca = cornerDoca;
}

void AcdCalAssoc::writeOut(MsgStream& stream) const
// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  stream << MSG::DEBUG
         << "AcdCalAssoc " << m_index << ' ' << (m_upward ? "up" : "down") << " E=" << m_energy << ' ' 
         << m_start << ' ' << m_dir << ' ' << " s= " << m_arcLength << " SSDVeto = " << m_tkrSSDVeto 
         << " CornerDoca = " << m_cornerDoca
         << std::endl
         << m_cov_start
         << m_cov_end
         << endreq;
}



void AcdCalAssoc::ini()
// Purpose: reset all data members to 0
//
{ 
  m_index = -1;
  m_upward = true;
  m_energy = 0.;
  m_start.set(0.,0.,0.);
  m_dir.set(0.,0.,0.);
  m_arcLength = 0.;
  m_cov_start = CLHEP::HepSymMatrix(5,0);
  m_cov_end = CLHEP::HepSymMatrix(5,0);
  m_tkrSSDVeto  = 0;
  m_cornerDoca = 0.;
  m_hitPocae.clear();
  m_gapPocae.clear();
}
 
AcdCalAssocCol::AcdCalAssocCol(const std::vector<AcdCalAssoc*>& acdCalAssocs) {
//Purpose: take ownership of CalAssocs from a vector
  for ( std::vector<AcdCalAssoc*>::const_iterator itr = acdCalAssocs.begin();
        itr != acdCalAssocs.end(); itr++ ) {
    AcdCalAssoc* CalAssoc = const_cast<AcdCalAssoc*>(*itr);
    add(CalAssoc);
  }
}

void AcdCalAssocCol::delCalAssocs()
//Purpose: delete all AcdCalAssoc object from memory
{
  int nCalAssoc = num();
  for (int iCalAssoc = 0; iCalAssoc < nCalAssoc; iCalAssoc++) {
    delete operator[](iCalAssoc);
  }
  clear();
}

void AcdCalAssocCol::ini()
//Purpose:  delete all pointers to clusters
// from collection 
{
  clear();
}


void AcdCalAssocCol::writeOut(MsgStream& stream) const
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
