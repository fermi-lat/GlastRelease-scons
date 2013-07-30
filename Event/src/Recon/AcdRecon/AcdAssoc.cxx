// File and Version information:
// $Header$
//
//  Implementation file of AcdAssoc,AcdAssocCol, and AcdCalAssocCol classes
//  
// Authors:
//
//    Eric Charles
//    Alex Drlica-Wagner
//

#include "Event/Recon/AcdRecon/AcdAssoc.h"

#include "GaudiKernel/MsgStream.h"

using namespace Event;

AcdAssoc::AcdAssoc()
  :m_index(-1),
   m_upward(true),
   m_energy(0.),
   m_start(0.,0.,0.),
   m_dir(0.,0.,0.),
   m_arcLength(0.),
   m_cov_start(5,0),
   m_cov_end(5,0),
   m_tkrSSDVeto(0),
   m_cornerDoca(0),
   m_point(0),
   m_energy15(0.),
   m_energy30(0.),
   m_energy45(0.),
   m_triggerEnergy15(0.),
   m_triggerEnergy30(0.),   
   m_triggerEnergy45(0.){
}


/// Constructor for use in reconstruction, 
AcdAssoc::AcdAssoc(int index, bool up, float energy, 
		   const HepPoint3D& start, const HepVector3D& dir, float arcLength,
		   const CLHEP::HepSymMatrix& covStart, const CLHEP::HepSymMatrix& covEnd,
		   int tkrSSDVeto, float cornerDoca,
		   float energy15, float energy30, float energy45,
		   float triggerEnergy15, float triggerEnergy30, float triggerEnergy45)		   
  :m_index(index),
   m_upward(up),
   m_energy(energy),
   m_start(start),
   m_dir(dir),
   m_arcLength(arcLength),
   m_cov_start(covStart),
   m_cov_end(covEnd),
   m_tkrSSDVeto(tkrSSDVeto),
   m_cornerDoca(cornerDoca),
   m_point(0),
   m_energy15(energy15),
   m_energy30(energy30),
   m_energy45(energy45),
   m_triggerEnergy15(triggerEnergy15),
   m_triggerEnergy30(triggerEnergy30),   
   m_triggerEnergy45(triggerEnergy45){      
}


AcdAssoc::~AcdAssoc(){
  if ( m_point != 0 ) {
    delete m_point;
  }
}

void AcdAssoc::set(int index, bool up, float energy, 
		   const HepPoint3D& start, const HepVector3D& dir, float arcLength,
		   const CLHEP::HepSymMatrix& covStart, const CLHEP::HepSymMatrix& covEnd,
		   int tkrSSDVeto, float cornerDoca,
		   float energy15, float energy30, float energy45,
		   float triggerEnergy15, float triggerEnergy30, float triggerEnergy45){		   
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
  m_energy15 = energy15;
  m_energy30 = energy30;
  m_energy45 = energy45;
  m_triggerEnergy15 = triggerEnergy15;
  m_triggerEnergy30 = triggerEnergy30;   
  m_triggerEnergy45 = triggerEnergy45;      
}

void AcdAssoc::writeOut(MsgStream& stream) const
// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{

  stream << MSG::DEBUG
         << "AcdAssoc " << m_index << ' ' << (m_upward ? "up" : "down") << " E=" << m_energy << ' ' 
         << m_start << ' ' << m_dir << ' ' << " s= " << m_arcLength << " SSDVeto = " << m_tkrSSDVeto 
         << " CornerDoca = " << m_cornerDoca
         << std::endl
         << m_cov_start
         << m_cov_end
	 << "Cone Energies: " << m_energy15 << ' ' << m_energy30 << ' ' << m_energy45 
	 << std::endl
	 << "Cone Eneriges (Trig): " << m_triggerEnergy15 << ' ' << m_triggerEnergy30 << ' ' << m_triggerEnergy45 
         << endreq;
}



void AcdAssoc::ini()
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
  if ( m_point != 0 ) {
    delete m_point;
  }
  m_point = 0;
  m_energy15 = 0.;
  m_energy30 = 0.;
  m_energy45 = 0.;
  m_triggerEnergy15 = 0.;
  m_triggerEnergy30 = 0.;
  m_triggerEnergy45 = 0.;
}
 
AcdAssocCol::AcdAssocCol(const std::vector<AcdAssoc*>& acdAssocs) {
//Purpose: take ownership of AcdAssocs from a vector
  for ( std::vector<AcdAssoc*>::const_iterator itr = acdAssocs.begin();
        itr != acdAssocs.end(); itr++ ) {
    AcdAssoc* acdAssoc = const_cast<AcdAssoc*>(*itr);
    add(acdAssoc);
  }
}

void AcdAssocCol::delAssocs()
//Purpose: delete all AcdAssoc object from memory
{
  int nAssoc = num();
  for (int iAssoc = 0; iAssoc < nAssoc; iAssoc++) {
    delete operator[](iAssoc);
  }
  clear();
}

void AcdAssocCol::ini()
//Purpose:  delete all pointers to clusters
// from collection 
{
  clear();
}


void AcdAssocCol::writeOut(MsgStream& stream) const
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
