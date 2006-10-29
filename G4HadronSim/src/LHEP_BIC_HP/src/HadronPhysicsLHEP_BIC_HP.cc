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
// $Id$
// GEANT4 tag $Name$
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsLHEP_BIC_HP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
//  1.12.2005 G.Folger: migration to non static particles
//
//----------------------------------------------------------------------------
//
#include "HadronPhysicsLHEP_BIC_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsLHEP_BIC_HP::HadronPhysicsLHEP_BIC_HP(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{}

void HadronPhysicsLHEP_BIC_HP::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theLHEPNeutron=new G4LHEPNeutronBuilder);
  theNeutrons->RegisterMe(theBinaryNeutron=new G4BinaryNeutronBuilder);
  theNeutrons->RegisterMe(theHPNeutron=new G4NeutronHPBuilder);
  theLHEPNeutron->SetMinEnergy(19.9*MeV);
  theLHEPNeutron->SetMinInelasticEnergy(9.5*GeV);
  theBinaryNeutron->SetMaxEnergy(9.9*GeV);
  theBinaryNeutron->SetMinEnergy(19.9*MeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theLHEPPro=new G4LHEPProtonBuilder);
  thePro->RegisterMe(theBinaryPro=new G4BinaryProtonBuilder);
  theLHEPPro->SetMinEnergy(9.5*GeV);
  theBinaryPro->SetMaxEnergy(9.9*GeV);
  
  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theLHEPPiK=new G4LHEPPiKBuilder);

  theMiscLHEP=new G4MiscLHEPBuilder;
  theStoppingHadron=new G4StoppingHadronBuilder;
  
}

HadronPhysicsLHEP_BIC_HP::~HadronPhysicsLHEP_BIC_HP()
{
    delete theNeutrons;
    delete theLHEPNeutron;
    delete theBinaryNeutron;
    delete theHPNeutron;
    
    delete thePiK;
    delete theLHEPPiK;
    
    delete thePro;
    delete theLHEPPro;
    delete theBinaryPro;
    
    delete theMiscLHEP;
    delete theStoppingHadron;
}

void HadronPhysicsLHEP_BIC_HP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP_BIC_HP::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();
  theStoppingHadron->Build();
}
// 2002 by J.P. Wellisch
