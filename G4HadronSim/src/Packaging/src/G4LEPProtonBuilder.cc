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
#include "G4HadronSim/G4LEPProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LEPProtonBuilder::
G4LEPProtonBuilder() 
{
  theMin = 0;
  theMax=55*GeV;
}

G4LEPProtonBuilder::
~G4LEPProtonBuilder() {}

void G4LEPProtonBuilder::
Build(G4HadronElasticProcess * aP)
{
  theElasticModel = new G4LElastic();
  aP->RegisterMe(theElasticModel);
}

void G4LEPProtonBuilder::
Build(G4ProtonInelasticProcess * aP)
{
// G4cout << "adding inelastic Proton in LHEP" << G4endl;
  theLEProtonModel = new G4LEProtonInelastic();
  theLEProtonModel->SetMinEnergy(theMin);
  theLEProtonModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theLEProtonModel);
}

// 2002 by J.P. Wellisch
