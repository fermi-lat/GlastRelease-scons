#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

    void setParticle(std::string pname);
    void setMomentum(G4ThreeVector pmom){particleGun->SetParticleMomentumDirection(pmom);}
    void setPosition(G4ThreeVector ppos){particleGun->SetParticlePosition(ppos);}
    // Set energy in MeV
    void setEnergy(G4double pen){particleGun->SetParticleEnergy(pen*MeV);}

  private:
    G4ParticleGun* particleGun;
};

#endif


