// File and Version Information:
// $Header$
//
// Description: this class is called by Geant4 to generate the primary particle
// during the event run
//
// Author(s):
//      R.Giannitrapani

#include "PrimaryGeneratorAction.h"

#include "Event/MonteCarlo/McParticle.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

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

void PrimaryGeneratorAction::init(Event::McParticle* part, IParticlePropertySvc* ppsvc)
{
  // Purpose and Method: this method set the particle by passing an
  // Event::McParticle pointer and the ParticlePropertySvc; 
  // Inputs: part is the pointer to the Event::McParticle
  // Inputs: ppsvc is the pointer to the IParticlePropertySvc

  Event::McParticle::StdHepId hepid= part->particleProperty();
  ParticleProperty* ppty = ppsvc->findByStdHepID( hepid );

  const HepLorentzVector& pfinal = part->finalFourMomentum();
  Hep3Vector dir=    pfinal.vect().unit();
  HepPoint3D p =   part->finalPosition();
  // note possibility of truncation error here! especially with MeV.
  double ke =   pfinal.e() - pfinal.m(); 
  
  // Set the G4 primary generator
  // the position has to be expressed in mm
  // while the energy in MeV

  // note the conversion from mass in Gev/c^2 to AMU in case of the ion.
  if (isIon(hepid))
    setIon((int)ppty->charge() , (int)ppty->mass()/0.93149);
  else
    setParticle(ppty->particle());
  
  setMomentum(dir);
  setPosition(p);
  setEnergy(ke);  
}

bool PrimaryGeneratorAction::isIon(int id)
{
  // Purpose and Method: this method return true if the id correspond to an ion
  // Inputs: id is the IDHEP identifier
  
  if ((id>=80000) && (id<90000))
    return 1;
  else return 0;
}

void PrimaryGeneratorAction::setIon(int pNatomic, int pMatomic, double pElevel)
{
  // Purpose and Method: this method set the ion by extracting it from the
  // particle table by using the atomic mass and number
  // Inputs: pNatomic is the atomic number
  // Inputs: pMatomic is the atomic mass
  // Inputs: pElevel is the excitation level of the ion (by default the ground state)

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* ion
    = particleTable->GetIon(pNatomic, pMatomic, pElevel);
  
  particleGun->SetParticleDefinition(ion);
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


