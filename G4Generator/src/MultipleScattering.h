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
// $Id:
// --------------------------------------------------------------
//    GEANT 4 class header file
//
// ***** This is the Geant4 v3.2 release version of MSC code ****
//       resurrected for use with G4 5.1 for testing purposes
//       6/3/03 Tracy Usher
//  
//    History: based on object model of
//    2nd December 1995, G.Cosmo
//    --------- MultipleScattering physics process --------
//               by Laszlo Urban, October 1997
// **************************************************************
//      UNIVERSAL: for arbitrary single charged particle
//  09/12/98:      charge can be != +- 1 !!!!   L.Urban
//  30/09/99:    nuclear size effect correction L.Urban 
// --------------------------------------------------------------
// 22/10/98: cleanup , L.Urban

#ifndef GLAST_MultipleScattering_h
#define GLAST_MultipleScattering_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4EnergyLossTables.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4PionPlus.hh"
#include "G4Proton.hh"
#include "G4PhysicsLogVector.hh"
#include "G4GPILSelection.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Material.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4UnitsTable.hh"

class MultipleScattering : public G4VContinuousDiscreteProcess
{
 public:

   MultipleScattering(const G4String& processName="msc") ;

   ~MultipleScattering() ;
          
   G4bool IsApplicable ( const G4ParticleDefinition& ) ;

   void SetPhysicsTableBining(G4double lowE,G4double highE,G4int nBins);

   void BuildPhysicsTable(const G4ParticleDefinition& aParticleType) ;

   void PrintInfoDefinition();

   G4double GetContinuousStepLimit(const G4Track& aTrack,
                                   G4double previousStepSize,
                                   G4double currentMinimumStep,
                                   G4double& currentSafety) ; 

   G4double GetMeanFreePath(const G4Track& aTrack,
                            G4double previousStepSize,
                            G4ForceCondition* condition) ;

   G4VParticleChange* AlongStepDoIt(const G4Track& aTrack,const G4Step& aStep);

   G4VParticleChange* PostStepDoIt(const G4Track& aTrack,const G4Step& aStep) ;

   G4double GetLambda(G4double KineticEnergy,G4Material* material);

   void SetScatteringParameter(G4double value)
           { scatteringparameter = value ; } ;
   void SetTuning(G4double value) { tuning = value ; };
   void SetCpar  (G4double value) { cpar   = value ; };
   void SetTlimitmsc  (G4double value) { Tlimit = value ; };
   void SetLateralDisplacementFlag(G4bool flag) {fLatDisplFlag = flag;};
   
   void SetNuclCorrPar(G4double val) { NuclCorrPar = val; } ;
   void SetFactPar(G4double val) { FactPar = val ; } ;

 protected:

   G4double ComputeTransportCrossSection(
                             const G4ParticleDefinition& aParticleType,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber,
                                   G4double AtomicWeight) ;

   G4double TrueToGeomTransformation(const G4DynamicParticle* aParticle,
                                     G4Material* aMaterial,
                                     G4double truePathLength) ;

 private:

 //  hide assignment operator as  private
   MultipleScattering & operator = (const MultipleScattering &right) ;
   MultipleScattering ( const MultipleScattering &) ;


 // data members ...................................................
 private:

   G4PhysicsTable* theTransportMeanFreePathTable ;

   G4double fTransportMeanFreePath ;
   G4double range,alpha1 ;
   G4int stepFlag ;

   G4double biglambda ;

   G4double LowestKineticEnergy ;
   G4double HighestKineticEnergy ;
   G4int TotBin ;

   const G4Electron* theElectron ;
   const G4Positron* thePositron ;

   G4Material* lastMaterial;
   G4double lastKineticEnergy;
   G4int materialIndex ;
  
   G4double tLast ;
   G4double zLast ;

   G4double Tlimit ;

   // model parameters
   G4double scatteringparameter;
   G4double tuning;
   G4double cpar;

   // with/without lateral displacement
   G4bool fLatDisplFlag ;

   // nuclear size effect correction
   G4double NuclCorrPar ;
   G4double FactPar ;

   //New ParticleChange
   G4ParticleChangeForMSC fParticleChange ;
   
};

#include "MultipleScattering.icc"

#endif
 

 
