///////////////////////////////////////////////////////////////////////////
//                                                                       
// The Tagger class contains the information about the ESA beam tagger  
//                                                                       
///////////////////////////////////////////////////////////////////////////

#include "digiRootData/Tagger.h"
#include "digiRootData/stuff.h"
#include "TClass.h"

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
  delete m_PID;
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
//________________________________________________________________________
/// Implement the streamer ourselves for now...to take advantage
/// of schema evolution - and to allow our TBEvent class to handle
/// both old (<= Root v2.25) Root files, and new (>= Root v3.00) files
void Tagger::Streamer(TBuffer &R__b)
{
   // Stream an object of class Tagger.

   if (R__b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); 
      if (R__v > 1) 
      { 
          Tagger::Class()->ReadBuffer(R__b, this, R__v, R__s, R__c);
          return;
      }
      /// Old Versions
      TObject::Streamer(R__b);
      R__b >> m_PID;
      R__b >> m_nBeamEnergy;
      R__b >> m_nMagCurrent;
   } else {
       Tagger::Class()->WriteBuffer(R__b, this);
   }
}