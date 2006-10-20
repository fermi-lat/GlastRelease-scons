//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Ion physics for a Linear Collider Detector                     //
//  Date:   6 July 2004                                                    //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef LCIonPhysics_h
#define LCIonPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "GlastMS/MultipleScatteringFactory.h"

class LCIonPhysics : public G4VPhysicsConstructor
{
  public: 
  //    LCIonPhysics(const G4String& name ="ion");

    LCIonPhysics(const G4String& name,GlastMS::MultipleScatteringFactory& msFactory );
    virtual ~LCIonPhysics();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

private: 

  GlastMS::MultipleScatteringFactory& m_msFactory;

};

#endif





