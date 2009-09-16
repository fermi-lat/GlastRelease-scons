//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Ion physics for Space Electronics Applications                 //
//  Date:   22 March 2005                                                  //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef SEIonPhysics_h
#define SEIonPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "GlastMS/MultipleScatteringFactory.h"

class SEIonPhysics : public G4VPhysicsConstructor
{
  public: 
  //    SEIonPhysics(const G4String& name ="ion");

    SEIonPhysics(const G4String& name,GlastMS::MultipleScatteringFactory& msFactory );
    virtual ~SEIonPhysics();

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





