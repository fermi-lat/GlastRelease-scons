#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "CLHEP/Vector/LorentzVector.h"

class G4Event;
class IParticlePropertySvc;
namespace Event{
class McParticle;
}

/** 
 * @class PrimaryGeneratorAction
 *
 * @brief A class to deal with primary particle generation
 *
 * This class implements the ordinary Geant4 mechanism to produce a primary
 * particle and inject it in the detector. This is a quite simple class and its
 * public interface is used primarly by G4Generator algorithm to set the
 * particle properties values coming from FluxSvc
 *  
 * @author R.Giannitrapani
 *    
 * $Header$
 */
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

  public:
    /// This implements a Geant4 superclass method (note the name convention);
    /// this is the main trigger to particle production in the detector and is
    /// called automatically by Geant4
    void GeneratePrimaries(G4Event*);

    /// This method init the generator with an McParticle coming from the
    /// FluxSvc.
    /// @param part The pointer to the Event::McParticle
    /// @param ppsvc The pointer to the ParticlePropertySvc of GAUDI
    void init(Event::McParticle* part, IParticlePropertySvc* ppsvc);

    /// @param pname The name of the particle
    void setParticle(std::string pname);
    /// @param pmom The 3d momentum of the particle
    void setMomentum(G4ThreeVector pmom){
      particleGun->SetParticleMomentumDirection(pmom);}
    /// @param ppos The 3d position (initial) of the particle
    void setPosition(G4ThreeVector ppos){
      particleGun->SetParticlePosition(ppos*mm);}
    /// @param pen Energy of the particle in MeV
    void setEnergy(G4double pen){particleGun->SetParticleEnergy(pen*MeV);}

    /// @return The PDG particle definition
    G4ParticleDefinition* GetParticleDefinition() {
      return particleGun->GetParticleDefinition();
    }
    
    /// @return The four momentum as a Lorentz vector
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


