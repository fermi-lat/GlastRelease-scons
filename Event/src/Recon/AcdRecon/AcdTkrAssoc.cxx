// File and Version information:
// $Header$
//
//  Implementation file of AcdTkrAssoc and AcdTkrAssocCol classes
//  
// Authors:
//
//    Eric Charles
//
//

#include "Event/Recon/AcdRecon/AcdTkrAssoc.h"

#include "GaudiKernel/MsgStream.h"

using namespace Event;

AcdTkrAssoc::AcdTkrAssoc()
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
AcdTkrAssoc::AcdTkrAssoc(const AcdTkrAssoc& other)
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
AcdTkrAssoc::AcdTkrAssoc(int index, bool up, float energy, 
                         const HepPoint3D& start, const HepVector3D& dir, float arcLength,
                         const HepSymMatrix& covStart, const HepSymMatrix& covEnd,
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


void AcdTkrAssoc::set(int index, bool up, float energy, 
                      const HepPoint3D& start, const HepVector3D& dir, float arcLength,
                      const HepSymMatrix& covStart, const HepSymMatrix& covEnd,
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

void AcdTkrAssoc::writeOut(MsgStream& stream) const
// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  stream << MSG::DEBUG
         << "AcdTkrAssoc " << m_index << ' ' << (m_upward ? "up" : "down") << " E=" << m_energy << ' ' 
         << m_start << ' ' << m_dir << ' ' << " s= " << m_arcLength << " SSDVeto = " << m_tkrSSDVeto 
         << " CornerDoca = " << m_cornerDoca
         << std::endl
         << m_cov_start
         << m_cov_end
         << endreq;
}



void AcdTkrAssoc::ini()
// Purpose: reset all data members to 0
//
{ 
  m_index = -1;
  m_upward = true;
  m_energy = 0.;
  m_start.set(0.,0.,0.);
  m_dir.set(0.,0.,0.);
  m_arcLength = 0.;
  m_cov_start = HepSymMatrix(5,0);
  m_cov_end = HepSymMatrix(5,0);
  m_tkrSSDVeto  = 0;
  m_cornerDoca = 0.;
  m_hitPocae.clear();
  m_gapPocae.clear();
}
 
AcdTkrAssocCol::AcdTkrAssocCol(const std::vector<AcdTkrAssoc*>& acdTkrAssocs) {
//Purpose: take ownership of TkrAssocs from a vector
  for ( std::vector<AcdTkrAssoc*>::const_iterator itr = acdTkrAssocs.begin();
        itr != acdTkrAssocs.end(); itr++ ) {
    AcdTkrAssoc* TkrAssoc = const_cast<AcdTkrAssoc*>(*itr);
    add(TkrAssoc);
  }
}

void AcdTkrAssocCol::delTkrAssocs()
//Purpose: delete all AcdTkrAssoc object from memory
{
  int nTkrAssoc = num();
  for (int iTkrAssoc = 0; iTkrAssoc < nTkrAssoc; iTkrAssoc++) {
    delete operator[](iTkrAssoc);
  }
  clear();
}

void AcdTkrAssocCol::ini()
//Purpose:  delete all pointers to clusters
// from collection 
{
  clear();
}


void AcdTkrAssocCol::writeOut(MsgStream& stream) const
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
