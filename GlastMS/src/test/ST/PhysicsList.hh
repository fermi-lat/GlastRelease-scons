// --------------------------------------------------------------
//      GEANT 4 -  2000
//
// Comments
//
//  Construct/define particles and physics processes
//
//  Particle used
//
//    Leptons: e+ e- mu+ mu- nue nue* numu numu*

//
//  Processes 
//
//    energy loss
//    multiple scattering
//    bremsstrahlung
//    compton
//    gamma conversion
//    photoelectric
//    ionization
//
//---------------------------------------------------------------

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class PhysicsList: public G4VUserPhysicsList
{
  public:
    PhysicsList();
    ~PhysicsList();

  protected:
    // Construct particle and physics process
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

 protected:
    // these methods Construct particles 
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();

  protected:
  // these methods Construct physics processes and register them
    void ConstructGeneral();
    void ConstructEM();

  private:
  // Data members to store the cuts
    G4double cutForGamma;
    G4double cutForElectron; 
    G4double cutForProton;

};
#endif







