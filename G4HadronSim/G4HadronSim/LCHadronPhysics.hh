//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Hadron physics for a Linear Collider Detector                  //
//  Date:   23 June 2004                                                   //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef LCHadronPhysics_h
#define LCHadronPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "GlastMS/MultipleScatteringFactory.h"


class LCHadronPhysics : public G4VPhysicsConstructor
{
public: 
  //LCHadronPhysics(const G4String& name="hadron");
  LCHadronPhysics(const G4String& name, GlastMS::MultipleScatteringFactory& msFactory);
  virtual ~LCHadronPhysics();

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





