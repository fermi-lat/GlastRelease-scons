// File and Version Information:
// $Header$
//
// Description: this class is called by Geant4 to generate the primary particle
// during the event run
//
// Author(s):
//      R.Giannitrapani

#include "PrimaryGeneratorAction.h"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  
  // A new Geant4 particle gun is created (with just one particle to be
  // produced)
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  // By default we use an electron starting in the origin, directed along the
  // negative z axis with 30 MeV of energy. 
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  particleGun->SetParticleEnergy(30.*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*mm,0.*mm,0.*mm));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::setParticle(std::string pname)
{
  // Purpose and Method: this method set the particle by searching a name in the
  // particle table (the name should be recognized as a valid name by Geant4)
  // Inputs: pname is the name of the particle
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
    = particleTable->FindParticle(pname);
  
  particleGun->SetParticleDefinition(particle);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Purpose and Method: this method is automatically called by the internal
  // mechanisms of Geant4 and it generates a primary particle in the simulation
  // Inputs: the G4Event pointer anEvent that represent the actual event
  particleGun->GeneratePrimaryVertex(anEvent);
}


