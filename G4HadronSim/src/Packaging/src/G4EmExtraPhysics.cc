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
// $Id$
// GEANT4 tag $Name$
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmExtraPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
//
//----------------------------------------------------------------------------
//

#include "G4EmExtraPhysics.hh"

#include "G4SynchrotronRadiation.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"

G4EmExtraPhysics::G4EmExtraPhysics(const G4String& name): G4VPhysicsConstructor(name),
  wasActivated(false), synchOn(false),
  gammNucOn(true), muNucOn(false), theElectronSynch(0), thePositronSynch(0),
  theGNPhysics(0),theMuNuc1(0),theMuNuc2(0)
{
  theMessenger = new G4EmMessenger(this);
}

G4EmExtraPhysics::~G4EmExtraPhysics()
{
  delete theMessenger;
  if(theElectronSynch) delete theElectronSynch;
  if(thePositronSynch) delete thePositronSynch;
  if(theGNPhysics)     delete theGNPhysics;
  if(theMuNuc1)        delete theMuNuc1;
  if(theMuNuc2)        delete theMuNuc2;
}

void G4EmExtraPhysics::Synch(G4String & newState)
{
  if(newState == "on") synchOn = true;
  else                 synchOn = false;
}

void G4EmExtraPhysics::GammaNuclear(G4String & newState)
{
  if(newState == "on") gammNucOn = true;
  else                 gammNucOn = false;
}

void G4EmExtraPhysics::MuonNuclear(G4String & newState)
{
  if(newState == "on") muNucOn = true;
  else                 muNucOn = false;
}

void G4EmExtraPhysics::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
}

void G4EmExtraPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;
  wasActivated = true;

  if(synchOn) {
    pManager = G4Electron::Electron()->GetProcessManager();
    theElectronSynch = new G4SynchrotronRadiation();
    pManager->AddDiscreteProcess(theElectronSynch);

    pManager = G4Positron::Positron()->GetProcessManager();
    thePositronSynch = new G4SynchrotronRadiation();
    pManager->AddDiscreteProcess(thePositronSynch);
  }
  if (gammNucOn) {
    theGNPhysics = new G4ElectroNuclearBuilder();
    theGNPhysics->Build();
  }
  if(muNucOn) {
    pManager  = G4MuonPlus::MuonPlus()->GetProcessManager();
    theMuNuc1 = new G4MuNuclearInteraction("muNucl");
    pManager->AddDiscreteProcess(theMuNuc1);

    pManager  = G4MuonMinus::MuonMinus()->GetProcessManager();
    theMuNuc2 = new G4MuNuclearInteraction("muNucl");
    pManager->AddDiscreteProcess(theMuNuc2);
  }
}
