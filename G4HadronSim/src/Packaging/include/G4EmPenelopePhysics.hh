//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id$
// GEANT4 tag $Name$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmPenelopePhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 05.12.2005 V.Ivanchenko add controlled verbosity
//
//----------------------------------------------------------------------------
//

#ifndef G4EmPenelopePhysics_h
#define G4EmPenelopePhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

// for lib list detection....
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4EmProcessOptions.hh"

#include "GlastMS/MultipleScatteringFactory.h"
#include "GlastMS/EnergyLossFactory.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmPenelopePhysics : public G4VPhysicsConstructor
{
public:
  //  G4EmPenelopePhysics(const G4String& name = "EMpenelope", G4int ver = 1,
  //	      G4bool msc=true);
  // to add GLAST MCS & Ionisation
  
  G4EmPenelopePhysics(const G4String& name, G4int ver,
		      GlastMS::MultipleScatteringFactory& msfactory,
		      GlastMS::EnergyLossFactory& eLossFactory );
		      //		      G4bool msc=true, 

  virtual ~G4EmPenelopePhysics();

public:
  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int  verbose;
  //  G4bool mscStepLimit;
  
  GlastMS::MultipleScatteringFactory& m_msFactory;
  GlastMS::EnergyLossFactory&         m_eLossFactory;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif






