// $Header$
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "CLHEP/Vector/LorentzVector.h"

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
    void setPosition(G4ThreeVector ppos){particleGun->SetParticlePosition(ppos*mm);}
    // Set energy in MeV
    void setEnergy(G4double pen){particleGun->SetParticleEnergy(pen*MeV);}

    G4ParticleDefinition* GetParticleDefinition() {
        return particleGun->GetParticleDefinition();
    }
    HepLorentzVector GetFourMomentum(){
        double mass = GetParticleDefinition()->GetPDGMass(),
            e = particleGun->GetParticleEnergy()+mass,
            p = sqrt(e*e-mass*mass);
        HepLorentzVector p4(e, 
            p*particleGun->GetParticleMomentumDirection());
        return p4;
    }

  private:
    G4ParticleGun* particleGun;
};

#endif


