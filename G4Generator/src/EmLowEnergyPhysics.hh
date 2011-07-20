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
// ClassName:   EmLowEnergyPhysics
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 05.12.2005 V.Ivanchenko add controlled verbosity
//
//----------------------------------------------------------------------------
//

#ifndef EmLowEnergyPhysics_h
#define EmLowEnergyPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

// for lib list detection....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EmLowEnergyPhysics : public G4VPhysicsConstructor
{
public:
  //  EmLowEnergyPhysics(const G4String& name = "EMlow", G4int ver = 1,
  //	      G4bool msc=true);
  // to add GLAST MCS & Ionisation
  
  EmLowEnergyPhysics(const G4String& name, G4int ver);
  virtual ~EmLowEnergyPhysics();

public:
  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int  verbose;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif






