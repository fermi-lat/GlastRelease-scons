#ifndef TAGGER_H
#define TAGGER_H

#include "ESAPID.h"
#include "TObject.h"

/*! \class Tagger
\brief
 
\li Jan 1999 Daniel Flath - ROOT HTML Documentation added
\li Dec 1999 Daniel Flath - Rewrite for GLAST
\li Dec 04,1999 Richard Dubois Clone from CalHit
*/

class Tagger : public TObject {		// so far it does not contain the raw data
 private:
  ESAPID* m_PID;            // Particle ID Info
  UInt_t m_nBeamEnergy;     // Beam Energy (MeV)
  UShort_t m_nMagCurrent;   // Tagger Magnet Current (deciAmps)
  
 protected:
  
 public:
  
  Tagger();
  Tagger(ESAPID* pid, UInt_t eBeam = 1, UShort_t iMag = 0);
  virtual ~Tagger();
  ESAPID* getPID() { return m_PID; };
    
  Float_t getBeamEnergy() const;
  Float_t getMagCurrent() const;
  void setBeamEnergy(UInt_t eBeam);  // (MeV)
  void setMagCurrent(UShort_t iMag); // (deciAmps)
  void setBeamEnergy(Float_t eBeam); // (GeV)
  void setMagCurrent(Float_t iMag);  // (Amps)

  ClassDef(Tagger,3)      // Tagger information
};

#endif
