//
////////////////////////////////////////////////////////////////
//                                                            //
//  Title:  Decay physics for a Linear Collider Detector      //
//  Date:   16 June 2004                                      //
//  Author: D.H. Wright (SLAC)                                //
//                                                            //
////////////////////////////////////////////////////////////////
//

#include "G4HadronSim/LCDecayPhysics.hh"

// #include "globals.hh"
// #include "G4ios.hh"
// #include <iomanip>   

LCDecayPhysics::LCDecayPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{
}

LCDecayPhysics::~LCDecayPhysics()
{
}

// #include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

void LCDecayPhysics::ConstructParticle()
{
}

void LCDecayPhysics::ConstructProcess()
{
  // Add Decay Process
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (fDecayProcess.IsApplicable(*particle)) { 
      pmanager ->AddProcess(&fDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(&fDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(&fDecayProcess, idxAtRest);
    }
  }
}


