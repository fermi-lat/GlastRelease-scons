#include "PrimaryGeneratorAction.h"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  particleGun->SetParticleEnergy(30.*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::setParticle(std::string pname)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
    = particleTable->FindParticle(pname);
  
  particleGun->SetParticleDefinition(particle);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->GeneratePrimaryVertex(anEvent);
}


