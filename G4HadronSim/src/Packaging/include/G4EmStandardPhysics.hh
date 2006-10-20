//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
// GEANT4 tag $Name$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 05.12.2005 V.Ivanchenko add controlled verbosity
//
//----------------------------------------------------------------------------
//

#ifndef G4EmStandardPhysics_h
#define G4EmStandardPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

// for lib list detection....
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4EmProcessOptions.hh"

#include "GlastMS/MultipleScatteringFactory.h"
#include "GlastMS/EnergyLossFactory.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmStandardPhysics : public G4VPhysicsConstructor
{
public:
  //  G4EmStandardPhysics(const G4String& name = "EMstandard", G4int ver = 1,
  //	      G4bool msc=true);
  // to add GLAST MCS & Ionisation
  
  G4EmStandardPhysics(const G4String& name, G4int ver,
		      GlastMS::MultipleScatteringFactory& msfactory,
		      GlastMS::EnergyLossFactory& eLossFactory );
		      //		      G4bool msc=true, 

  virtual ~G4EmStandardPhysics();

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






