//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Hadron physics for Space Electronics Applications              //
//  Date:   20 March 2005                                                  //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef SEHadronPhysics_h
#define SEHadronPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "GlastMS/MultipleScatteringFactory.h"


class SEHadronPhysics : public G4VPhysicsConstructor
{
public: 
  //SEHadronPhysics(const G4String& name="hadron");
  SEHadronPhysics(const G4String& name, GlastMS::MultipleScatteringFactory& msFactory);
  virtual ~SEHadronPhysics();

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





