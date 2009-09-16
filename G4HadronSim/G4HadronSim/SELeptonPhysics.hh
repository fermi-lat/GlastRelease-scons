//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//  Title:  Lepton physics for Space Electronics Applications              //
//  Date:   22 March 2005                                                  //
//  Author: D.H. Wright (SLAC)                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

#ifndef SELeptonPhysics_h
#define SELeptonPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "GlastMS/MultipleScatteringFactory.h"
#include "GlastMS/EnergyLossFactory.h"

class SELeptonPhysics : public G4VPhysicsConstructor
{
  public: 
  //   SELeptonPhysics(const G4String& name ="lepton");
  SELeptonPhysics(const G4String& name, GlastMS::MultipleScatteringFactory& msfactory,
		  GlastMS::EnergyLossFactory& eLossFactory); 
  virtual ~SELeptonPhysics();
  
    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
    virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
    virtual void ConstructProcess();

private:
  
  GlastMS::MultipleScatteringFactory& m_msFactory;
  GlastMS::EnergyLossFactory&         m_eLossFactory;

};

#endif





