///////////////////////////////////////////////////////////////////////////
//                                                                       
// The Tagger class contains the information about the ESA beam tagger  
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/Tagger.h"
#include "digiRootData/stuff.h"
//#include "TClass.h"

ClassImp(Tagger)

///________________________________________________________________________
Tagger::Tagger() 
:m_PID(0)  // protect destructor
{
  // Default constructor
  m_nBeamEnergy = 1;  // ESA code doesn't like energies of zero
  m_nMagCurrent = 0;
}
///________________________________________________________________________
Tagger::Tagger(ESAPID* pid, UInt_t eBeam, UShort_t iMag) : 
  m_PID(pid),
  m_nBeamEnergy(eBeam),
  m_nMagCurrent(iMag) {
  // Create a Tagger object with user specified pid, beamEnergy, and magCurrent
}
//_________________________________________________________________________
Tagger::~Tagger() {
  // Destructor, free up reserved memory
    if (m_PID) {
        delete m_PID;
        m_PID = 0;
    }
}
//_________________________________________________________________________
Float_t Tagger::getBeamEnergy() const { 
  // Returns the Beam Energy in GeV
  return (Float_t)m_nBeamEnergy / 1000.0;
}
//________________________________________________________________________
Float_t Tagger::getMagCurrent() const { 
  // Returns the Tagger Magnet Current in Amps
  return (Float_t)m_nMagCurrent / 10.0;
}
//________________________________________________________________________
void Tagger::setBeamEnergy(UInt_t eBeam) { 
  // Sets the Beam Energy to value (MeV) specified in iMag
  m_nBeamEnergy = eBeam;
}
//________________________________________________________________________
void Tagger::setMagCurrent(UShort_t iMag) { 
  // Sets the Tagger Magnet Current to value (deciAmps) specified in iMag
  m_nMagCurrent = iMag;
}
//________________________________________________________________________
void Tagger::setBeamEnergy(Float_t eBeam) {
  m_nBeamEnergy = stuff::Round(eBeam * 1000.0);
}
//________________________________________________________________________
void Tagger::setMagCurrent(Float_t iMag) {
  m_nMagCurrent = stuff::Round(iMag * 10.0);
}
