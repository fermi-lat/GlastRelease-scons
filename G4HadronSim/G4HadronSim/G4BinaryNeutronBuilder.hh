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
#ifndef G4BinaryNeutronBuilder_h
#define G4BinaryNeutronBuilder_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4VNeutronBuilder.hh"

#include "G4BinaryCascade.hh"   
#include "G4NeutronInelasticCrossSection.hh"

class G4BinaryNeutronBuilder : public G4VNeutronBuilder
{
  public: 
    G4BinaryNeutronBuilder();
    virtual ~G4BinaryNeutronBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess * aP);
    virtual void Build(G4HadronFissionProcess * aP);
    virtual void Build(G4HadronCaptureProcess * aP);
    virtual void Build(G4NeutronInelasticProcess * aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}
    void SetMaxEnergy(G4double aM) {theMax = aM;}

  private:
    G4BinaryCascade * theModel;    
    G4NeutronInelasticCrossSection theXSec;
    G4double theMin;
    G4double theMax;

};

// 2002 by J.P. Wellisch

#endif

