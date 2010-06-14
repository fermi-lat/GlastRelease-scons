// --------------------------------------------------------------
//      GEANT 4 -  2000
//---------------------------------------------------------------

#ifndef GeneratorAction_h
#define GeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class GeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
   GeneratorAction();
    ~GeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4ParticleGun* particleGun;
};

#endif


