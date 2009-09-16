//
////////////////////////////////////////////////////////////////
//                                                            //
//  Title:  Decay physics for Space Electronics Applications  //
//  Date:   20 March 2005                                     //
//  Author: D.H. Wright (SLAC)                                //
//                                                            //
////////////////////////////////////////////////////////////////
//

#include "G4HadronSim/SEDecayPhysics.hh"
#include "G4ProcessManager.hh"


SEDecayPhysics::SEDecayPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{;}

SEDecayPhysics::~SEDecayPhysics()
{;}

void SEDecayPhysics::ConstructParticle()
{;}

void SEDecayPhysics::ConstructProcess()
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


