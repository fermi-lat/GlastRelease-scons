// File and Version Information:
// $Header$
//
// Description: this class is called by Geant4 to generate the primary particle
// during the event run
//
// Author(s):
//      R.Giannitrapani

#include "PrimaryGeneratorAction.h"

//#include "Event/MonteCarlo/McParticle.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4Generator/G4GenException.h"

// TU: Hack for CLHEP 1.9.2.2
typedef HepGeom::Point3D<double>  HepPoint3D;
typedef HepGeom::Vector3D<double> HepVector3D;

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
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  //std::cout << particle->GetParticleName() << std::endl;
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  particleGun->SetParticleEnergy(30.*MeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*mm,0.*mm,0.*mm));

  m_primaryVertex = 0;
  m_secondaryVertexVec.clear();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::init(Event::McParticleCol* pcol, IParticlePropertySvc* ppsvc, double dz)
{
    // Purpose and Method: this method sets the particle by passing an
    // Event::McParticle pointer and the ParticlePropertySvc; 
    // Inputs: part is the pointer to the Event::McParticle
    // Inputs: ppsvc is the pointer to the IParticlePropertySvc
    // Inputs: dz is an optional z offset of the particle, to match LAT offset
        
    // Count the number of particles input

    // First particle is the "generator" particle... 
    // This particle defines the event vertex
    Event::McParticle* primary    = pcol->front();
    HepPoint3D         primVtxPos = primary->initialPosition();

    // Look for sources with energy over the limit
    double maxAllowedEnergy = 20001000.; // Temporary until I remember how to extract this from G4!
    if (primary->initialFourMomentum().e() > maxAllowedEnergy)
    {
        std::stringstream errorStr(" ");
        errorStr << "PrimaryGeneratorAction found energy over maximum allowed: " 
            << primary->initialFourMomentum().e() << ", source: "  
            << primary->getProcess() << std::endl;
        throw G4GenException(errorStr.str());
    }

    // Adjust by the delta z input
    primVtxPos = CLHEP::Hep3Vector(primVtxPos.x(), primVtxPos.y(), primVtxPos.z() + dz);

    // Create the vertex for this primary particle
    m_primaryVertex = new G4PrimaryVertex(primVtxPos, 0.);

    // Is the first particle to be added to the list of particles to be tracked?
    if (!(primary->statusFlags() & Event::McParticle::NOTTRACK)) 
    {
        m_primaryVertex->SetPrimary(convertToG4Primary(primary,ppsvc));
    }
    
    // Retrieve the vector of daughter particles
    const SmartRefVector<Event::McParticle>& daughterVec = primary->daughterList();
    
    SmartRefVector<Event::McParticle>::const_iterator dIter = daughterVec.begin();

    m_secondaryVertexVec.clear();

    //for(Event::McParticleCol::iterator colIter = pcol->begin(); colIter != pcol->end(); colIter++)
    for(dIter = daughterVec.begin(); dIter != daughterVec.end(); dIter++)
    {
        const Event::McParticle* mcPart  = *dIter;

        HepPoint3D scndVtxPos = mcPart->initialPosition();

        // Adjust by the delta z input
        scndVtxPos = CLHEP::Hep3Vector(scndVtxPos.x(), scndVtxPos.y(), scndVtxPos.z() + dz);

        G4PrimaryVertex* scndParticle = new G4PrimaryVertex(scndVtxPos, 0.);

        scndParticle->SetPrimary(convertToG4Primary(mcPart, ppsvc));

        m_secondaryVertexVec.push_back(scndParticle);

        //m_primaryVertex->SetPrimary(convertToG4Primary(mcPart, ppsvc));
    }

    return;
}
    
G4PrimaryParticle* PrimaryGeneratorAction::convertToG4Primary(const Event::McParticle* mcPart, IParticlePropertySvc* ppsvc)
{
    // Retrieve and convert the particle definition
    Event::McParticle::StdHepId hepid   = mcPart->particleProperty();
    G4ParticleDefinition*       partDef = 0;

    if (isIon(hepid))
    {
        //      std::cout << " Ion " << (int)ppty->charge() << " " << (int)ppty->mass()/931.49 << std::endl;
        int z=0;
        int a=0;
        a = (hepid - 1e9)/1e6;
        z = (hepid - 1e9 -a*1e6)/1e3;

        partDef = G4ParticleTable::GetParticleTable()->GetIon(z, a, 0.0);
    }
    else
    {
        ParticleProperty* ppty = ppsvc->findByStdHepID( hepid );
        //partDef = G4ParticleTable::GetParticleTable()->FindParticle(ppty->particle());
        partDef = G4ParticleTable::GetParticleTable()->FindParticle(ppty->pdgID());
    }

    // Look for a null pointer
    if (!partDef)
    {
        std::stringstream errorStr(" ");
        errorStr << "PrimaryGeneratorAction failed to find particle id: " << hepid 
            << " in the G4ParticleTable " << std::endl;
        throw G4GenException(errorStr.str());
    }

    // Position and momentum
    const CLHEP::HepLorentzVector& pinitial = mcPart->initialFourMomentum();

    // New Primary Particle
    G4PrimaryParticle* primPart = new G4PrimaryParticle(partDef, pinitial.px(), pinitial.py(), pinitial.pz());
  
    primPart->SetMass(partDef->GetPDGMass());
    primPart->SetCharge(partDef->GetPDGCharge());

    return primPart;
}

void PrimaryGeneratorAction::init(Event::McParticle* part, IParticlePropertySvc* ppsvc, double dz)
{
  // Purpose and Method: this method sets the particle by passing an
  // Event::McParticle pointer and the ParticlePropertySvc; 
  // Inputs: part is the pointer to the Event::McParticle
  // Inputs: ppsvc is the pointer to the IParticlePropertySvc
  // Inputs: dz is an optional z offset of the particle, to match LAT offset

  Event::McParticle::StdHepId hepid= part->particleProperty();
  ParticleProperty* ppty = ppsvc->findByStdHepID( hepid );

  const CLHEP::HepLorentzVector& pfinal = part->finalFourMomentum();
  CLHEP::Hep3Vector dir = pfinal.vect().unit();
  HepPoint3D        p   = part->finalPosition();
  p = CLHEP::Hep3Vector(p.x(), p.y(), p.z()+dz);
  // note possibility of truncation error here! especially with MeV.
  double ke =   pfinal.e() - pfinal.m(); 
  
  // Set the G4 primary generator
  // the position has to be expressed in mm
  // while the energy in MeV

  // note the conversion from mass in Mev/c^2 to AMU in case of the ion.
  if (isIon(hepid))
    {
      //      std::cout << " Ion " << (int)ppty->charge() << " " << (int)ppty->mass()/931.49 << std::endl;
      int z=0;
      int a=0;
      a = (hepid - 1e9)/1e6;
      z = (hepid - 1e9 -a*1e6)/1e3;
      
      //      std::cout<< " Ion new " << z << " " << a << std::endl;
      //      setIon((int)ppty->charge() , (int)ppty->mass()/931.49);
      setIon(z,a);
    }
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
  
  //  if ((id>=80000) && (id<90000))
  if (id>1e9) // new convention
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
  //if (m_primaryVertexVec.empty()) particleGun->GeneratePrimaryVertex(anEvent);
  //else
  //{
  //    for(std::vector<G4PrimaryVertex*>::iterator vtxIter = m_primaryVertexVec.begin();
  //        vtxIter != m_primaryVertexVec.end(); vtxIter++)
  //    {
  //        anEvent->AddPrimaryVertex(*vtxIter);
  //    }
    anEvent->AddPrimaryVertex(m_primaryVertex);

    // If secondaries then add those as well
    for(std::vector<G4PrimaryVertex*>::iterator vtxIter = m_secondaryVertexVec.begin();
        vtxIter != m_secondaryVertexVec.end(); vtxIter++)
    {
        anEvent->AddPrimaryVertex(*vtxIter);
    }
}


