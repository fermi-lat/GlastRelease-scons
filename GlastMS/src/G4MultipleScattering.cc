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
//    GEANT 4 class implementation file
//
//    ***** This is the Geant4 v3.2 release version of MSC ******
//   Incorporated into the GLAST simulation for testing 6/3/03
//   Slight modifications needed (as noted) to handle updates to 
//   Geant4 version 5.1
//   Tracy Usher (blame him)
//
//    History: based on object model of
//    2nd December 1995, G.Cosmo
//   -------- G4MultipleScattering physics process ------------
//             by Laszlo Urban, October 1997
// **************************************************************
// 24/10/97  correction in PostStepDoIt for tau << 1. L.Urban
// 09/12/98: charge can be different from +-1 !!!! L.Urban
// 29/07/99: corr. for low energy , L.Urban
// 17/09/99: corr. for high energy and/or small step , L.Urban
// 30/09/99: nuclear size effect correction,  L.Urban
// 22/03/00: value of member cpar has changed! , L.Urban
// 20/06/00: nuclear size correction for particles other than e+/e- only ,  L.Urban
// 10/08/00 values of some data members has been changed, L.Urban
// 09/11/00 bug corrected in sigma computation, L.Urban
// 16/05/01 value of cpar changed back to the old value, L.Urban
// 18/05/01 V.Ivanchenko Clean up againist Linux ANSI compilation 
// --------------------------------------------------------------
// note changed to .h
#include "GlastMS/G4MultipleScattering.h"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
// GLAST change: in namespace
namespace GlastMS{
  G4MultipleScattering::G4MultipleScattering(const G4String& processName)
     : G4VContinuousDiscreteProcess(processName),
       theTransportMeanFreePathTable(0),
       fTransportMeanFreePath (1.e12),
       range(1.e10*mm),
       alpha1(5.),
       stepFlag(0),
       biglambda (1.e10*mm), 
       LowestKineticEnergy(0.1*keV),
       HighestKineticEnergy(100.*TeV),
       TotBin(100),
       theElectron(G4Electron::Electron()),
       thePositron(G4Positron::Positron()),
       lastMaterial(0),
       lastKineticEnergy(0.),
       materialIndex(0),
       tLast (0.0),
       zLast (0.0),
       Tlimit(1.*keV),
       scatteringparameter(0.9),
       tuning (1.00),
       cpar (1.50),
       fLatDisplFlag(true), 
       NuclCorrPar (0.0615),
       FactPar(0.40)
  { }

  G4MultipleScattering::~G4MultipleScattering()
  {
    if(theTransportMeanFreePathTable)
    {
      theTransportMeanFreePathTable->clearAndDestroy() ;
      delete theTransportMeanFreePathTable ;
    }
  }

//  .........methods..............................

  void G4MultipleScattering::BuildPhysicsTable(
                              const G4ParticleDefinition& aParticleType)

  //   tables are built for MATERIALS
  {
    const G4double sigmafactor = twopi*classic_electr_radius*
                                       classic_electr_radius ;
    G4double KineticEnergy,AtomicNumber,AtomicWeight,
             sigma,lambda ;
    G4double density ;

  // destroy old tables if any
    if(theTransportMeanFreePathTable)
    {
      theTransportMeanFreePathTable->clearAndDestroy() ;
      delete theTransportMeanFreePathTable ;
    }

  //  create table
    const G4MaterialTable* theMaterialTable=
                                   G4Material::GetMaterialTable();
//**** G4MaterialTable Definition changed - TU
//    G4int numOfMaterials = theMaterialTable->length() ;
    G4int numOfMaterials = theMaterialTable->size() ;

    theTransportMeanFreePathTable = new G4PhysicsTable(numOfMaterials) ;

  //  loop for materials
//****  Since G4MaterialTable is now typedef'd as a vector, change the loop
//      for (G4int J=0; J<numOfMaterials; J++)
    for(G4MaterialTable::const_iterator matIter = theMaterialTable->begin(); 
        matIter != theMaterialTable->end(); matIter++)
    {

  //  create physics vector and fill it
      G4PhysicsLogVector* aVector = new G4PhysicsLogVector(
                          LowestKineticEnergy,HighestKineticEnergy,TotBin) ;

  //  get elements in the material
//****      const G4Material* material = (*theMaterialTable)(J) ;
//will this change work with the gcc compiler?
      const G4Material* material = *matIter;
      const G4ElementVector* theElementVector =
                             material->GetElementVector() ;
      const G4double* theAtomicNumDensityVector =
                             material->GetAtomicNumDensityVector() ;       
      const G4int NumberOfElements =
                             material->GetNumberOfElements() ;
      density = material->GetDensity() ;
 
  //  loop for kinetic energy values
      for (G4int i=0; i<TotBin; i++)
      {
          KineticEnergy = aVector->GetLowEdgeEnergy(i) ;
          sigma = 0. ;

  //   loop for element in the material
          for (G4int iel=0; iel<NumberOfElements; iel++)
          {
  //**** Change () in the following to [] to compile with vc7
            AtomicNumber = (*theElementVector)[iel]->GetZ() ;
            AtomicWeight = (*theElementVector)[iel]->GetA() ;
            sigma += theAtomicNumDensityVector[iel]*
                     ComputeTransportCrossSection(aParticleType,
                                     KineticEnergy,
                                     AtomicNumber,AtomicWeight) ;
          }
          sigma *= sigmafactor ;
 
          lambda = 1./sigma ;

        aVector->PutValue(i,lambda) ;
      }

     theTransportMeanFreePathTable->insert(aVector) ;

    }

    if(   (&aParticleType == G4Electron::Electron())  ||
          (&aParticleType == G4MuonPlus::MuonPlus())  ||
          (&aParticleType == G4Proton::Proton())         ) 
    {
      PrintInfoDefinition() ;
    }
  }


  G4double G4MultipleScattering::ComputeTransportCrossSection(
                   const G4ParticleDefinition& aParticleType,
                         G4double KineticEnergy,
                         G4double AtomicNumber,G4double AtomicWeight)
  {
    const G4double epsfactor = 2.*electron_mass_c2*electron_mass_c2*
                               Bohr_radius*Bohr_radius/(hbarc*hbarc) ;
    const G4double epsmin = 1.e-4 , epsmax = 1.e10 ;


    const G4double Zdat[15] = {4.,6.,13.,20.,26.,29.,32.,38.,47.,
                               50.,56.,64.,74.,79.,82. } ;

    const G4double Tdat[22] = 
                      { 0.0001*MeV,0.0002*MeV,0.0004*MeV,0.0007*MeV,
                       0.001*MeV,0.002*MeV,0.004*MeV,0.007*MeV,0.01*MeV,
                       0.02*MeV,0.04*MeV,0.07*MeV,0.1*MeV,0.2*MeV,
                       0.4*MeV,0.7*MeV,1.*MeV,2.*MeV,4.*MeV,
                       7.*MeV,10.*MeV,20.*MeV} ;

  // corr. factors for e-/e+ lambda

    const G4double celectron[15][22] =
          {{1.125,1.072,1.051,1.047,1.047,1.050,1.052,1.054,
            1.054,1.057,1.062,1.069,1.075,1.090,1.105,1.111,
            1.112,1.108,1.100,1.093,1.089,1.087            },
           {1.408,1.246,1.143,1.096,1.077,1.059,1.053,1.051,
            1.052,1.053,1.058,1.065,1.072,1.087,1.101,1.108,
            1.109,1.105,1.097,1.090,1.086,1.082            },
           {2.833,2.268,1.861,1.612,1.486,1.309,1.204,1.156,
            1.136,1.114,1.106,1.106,1.109,1.119,1.129,1.132,
            1.131,1.124,1.113,1.104,1.099,1.098            },
           {3.879,3.016,2.380,2.007,1.818,1.535,1.340,1.236,
            1.190,1.133,1.107,1.099,1.098,1.103,1.110,1.113,
            1.112,1.105,1.096,1.089,1.085,1.098            },
           {6.937,4.330,2.886,2.256,1.987,1.628,1.395,1.265,
            1.203,1.122,1.080,1.065,1.061,1.063,1.070,1.073,
            1.073,1.070,1.064,1.059,1.056,1.056            },
           {9.616,5.708,3.424,2.551,2.204,1.762,1.485,1.330,
            1.256,1.155,1.099,1.077,1.070,1.068,1.072,1.074,
            1.074,1.070,1.063,1.059,1.056,1.052            },
           {11.72,6.364,3.811,2.806,2.401,1.884,1.564,1.386,
            1.300,1.180,1.112,1.082,1.073,1.066,1.068,1.069,
            1.068,1.064,1.059,1.054,1.051,1.050            },
           {18.08,8.601,4.569,3.183,2.662,2.025,1.646,1.439,
            1.339,1.195,1.108,1.068,1.053,1.040,1.039,1.039,
            1.039,1.037,1.034,1.031,1.030,1.036            },
           {18.22,10.48,5.333,3.713,3.115,2.367,1.898,1.631,
            1.498,1.301,1.171,1.105,1.077,1.048,1.036,1.033,
            1.031,1.028,1.024,1.022,1.021,1.024            },
           {14.14,10.65,5.710,3.929,3.266,2.453,1.951,1.669,
            1.528,1.319,1.178,1.106,1.075,1.040,1.027,1.022,
            1.020,1.017,1.015,1.013,1.013,1.020            },
           {14.11,11.73,6.312,4.240,3.478,2.566,2.022,1.720,
            1.569,1.342,1.186,1.102,1.065,1.022,1.003,0.997,
            0.995,0.993,0.993,0.993,0.993,1.011            },
           {22.76,20.01,8.835,5.287,4.144,2.901,2.219,1.855,
            1.677,1.410,1.224,1.121,1.073,1.014,0.986,0.976,
            0.974,0.972,0.973,0.974,0.975,0.987            },
           {50.77,40.85,14.13,7.184,5.284,3.435,2.520,2.059,
            1.837,1.512,1.283,1.153,1.091,1.010,0.969,0.954,
            0.950,0.947,0.949,0.952,0.954,0.963            },
           {65.87,59.06,15.87,7.570,5.567,3.650,2.682,2.182,
            1.939,1.579,1.325,1.178,1.108,1.014,0.965,0.947,
            0.941,0.938,0.940,0.944,0.946,0.954            },
          // {45.60,47.34,15.92,7.810,5.755,3.767,2.760,2.239, // paper.....
           {55.60,47.34,15.92,7.810,5.755,3.767,2.760,2.239, 
            1.985,1.609,1.343,1.188,1.113,1.013,0.960,0.939,
            0.933,0.930,0.933,0.936,0.939,0.949            }};

     const G4double cpositron[15][22] = {
           {2.589,2.044,1.658,1.446,1.347,1.217,1.144,1.110,
            1.097,1.083,1.080,1.086,1.092,1.108,1.123,1.131,
            1.131,1.126,1.117,1.108,1.103,1.100            },
           {3.904,2.794,2.079,1.710,1.543,1.325,1.202,1.145,
            1.122,1.096,1.089,1.092,1.098,1.114,1.130,1.137,
            1.138,1.132,1.122,1.113,1.108,1.102            },
           {7.970,6.080,4.442,3.398,2.872,2.127,1.672,1.451,
            1.357,1.246,1.194,1.179,1.178,1.188,1.201,1.205,
            1.203,1.190,1.173,1.159,1.151,1.145            },
           {9.714,7.607,5.747,4.493,3.815,2.777,2.079,1.715,
            1.553,1.353,1.253,1.219,1.211,1.214,1.225,1.228,
            1.225,1.210,1.191,1.175,1.166,1.174            },
           {17.97,12.95,8.628,6.065,4.849,3.222,2.275,1.820,
            1.624,1.382,1.259,1.214,1.202,1.202,1.214,1.219,
            1.217,1.203,1.184,1.169,1.160,1.151            },
           {24.83,17.06,10.84,7.355,5.767,3.707,2.546,1.996,
            1.759,1.465,1.311,1.252,1.234,1.228,1.238,1.241,
            1.237,1.222,1.201,1.184,1.174,1.159            },
           {23.26,17.15,11.52,8.049,6.375,4.114,2.792,2.155,
            1.880,1.535,1.353,1.281,1.258,1.247,1.254,1.256,
            1.252,1.234,1.212,1.194,1.183,1.170            },
           {22.33,18.01,12.86,9.212,7.336,4.702,3.117,2.348,
            2.015,1.602,1.385,1.297,1.268,1.251,1.256,1.258,
            1.254,1.237,1.214,1.195,1.185,1.179            },
           {33.91,24.13,15.71,10.80,8.507,5.467,3.692,2.808,
            2.407,1.873,1.564,1.425,1.374,1.330,1.324,1.320,
            1.312,1.288,1.258,1.235,1.221,1.205            },
           {32.14,24.11,16.30,11.40,9.015,5.782,3.868,2.917,
            2.490,1.925,1.596,1.447,1.391,1.342,1.332,1.327,
            1.320,1.294,1.264,1.240,1.226,1.214            },
           {29.51,24.07,17.19,12.28,9.766,6.238,4.112,3.066,
            2.602,1.995,1.641,1.477,1.414,1.356,1.342,1.336,
            1.328,1.302,1.270,1.245,1.231,1.233            },
           {38.19,30.85,21.76,15.35,12.07,7.521,4.812,3.498,
            2.926,2.188,1.763,1.563,1.484,1.405,1.382,1.371,
            1.361,1.330,1.294,1.267,1.251,1.239            },
           {49.71,39.80,27.96,19.63,15.36,9.407,5.863,4.155,
            3.417,2.478,1.944,1.692,1.589,1.480,1.441,1.423,
            1.409,1.372,1.330,1.298,1.280,1.258            },
           {59.25,45.08,30.36,20.83,16.15,9.834,6.166,4.407,
            3.641,2.648,2.064,1.779,1.661,1.531,1.482,1.459,
            1.442,1.400,1.354,1.319,1.299,1.272            },
           {56.38,44.29,30.50,21.18,16.51,10.11,6.354,4.542,
            3.752,2.724,2.116,1.817,1.692,1.554,1.499,1.474,
            1.456,1.412,1.364,1.328,1.307,1.282            }};

     G4double Z23,ParticleMass,rat2,Charge,TotalEnergy,beta2,bg2,
              eps,Z1,Z2,ratZ,T,E,b2small,b2big,ratb2,c1,c2,cc1,cc2,
              corr,sigma,corrfactor,ChargeSquare ;
     G4int iZ,iT ;

     Z23 = 2.*log(AtomicNumber)/3. ;
     Z23 = exp(Z23) ;

     ParticleMass = aParticleType.GetPDGMass() ;

     rat2 = ParticleMass/electron_mass_c2 ;
     rat2 = rat2*rat2 ;

     Charge = aParticleType.GetPDGCharge() ;
     ChargeSquare = Charge*Charge/(eplus*eplus) ;

     TotalEnergy = KineticEnergy + ParticleMass ;

     beta2 = KineticEnergy*(TotalEnergy+ParticleMass)/
             (TotalEnergy*TotalEnergy) ;
     bg2 = KineticEnergy*(TotalEnergy+ParticleMass)/
           (ParticleMass*ParticleMass) ;

     eps = rat2*epsfactor*bg2/Z23 ;

     if(eps<epsmin)
       sigma = 2.*eps*eps ;
     else if(eps<epsmax)
       sigma = log(1.+2.*eps)-2.*eps/(1.+2.*eps) ;
     else
       sigma = log(2.*eps)-1.+1./eps ;

     sigma *=ChargeSquare*AtomicNumber*AtomicNumber/rat2 ;
     sigma /= beta2*bg2 ;

  // nuclear size effect correction for high energy
  // ( a simple approximation at present)
     G4double corrnuclsize,a,x0,w1,w2,w ;

     x0 = 1. - NuclCorrPar*ParticleMass/(KineticEnergy*
               exp(log(AtomicWeight/(g/mole))/3.)) ;
     if((x0 < -1.) || (KineticEnergy  <= 10.*MeV))
     {
       x0=-1. ;
       corrnuclsize = 1. ;
     }
     else
     {
     a = 1.+1./eps ;
     if(eps > epsmax)
       w1=log(2.*eps)+1./eps-3./(8.*eps*eps) ;
     else
       w1=log((a+1.)/(a-1.))-2./(a+1.) ;

     w = 1./((1.-x0)*eps) ;
     if(w < epsmin)
       w2=-log(w)-1.+2.*w-1.5*w*w ;
     else
       w2 = log((a-x0)/(a-1.))-(1.-x0)/(a-x0) ;
 
     corrnuclsize = w1/w2 ;

     // ####################################################
     corrnuclsize = exp(-FactPar*proton_mass_c2/KineticEnergy)*
                    (corrnuclsize-1.)+1. ;
     }

  //  correct this value using the corrections computed for e+/e-
     KineticEnergy *= electron_mass_c2/ParticleMass ;

     //  interpolate in AtomicNumber and beta2

     //  get bin number in Z
       iZ = 14 ;

       while ((iZ>=0)&&(Zdat[iZ]>=AtomicNumber)) iZ -= 1 ;
       if (iZ==14)                               iZ = 13 ;
       if (iZ==-1)                               iZ = 0 ;

       Z1 = Zdat[iZ] ;
       Z2 = Zdat[iZ+1] ;
       ratZ = (AtomicNumber-Z1)/(Z2-Z1) ;

     // get bin number in T (beta2)
       iT = 21 ;
       while ((iT>=0)&&(Tdat[iT]>=KineticEnergy))
         iT -= 1 ;
       if(iT==21)
         iT = 20 ;
       if(iT==-1)
         iT = 0 ;

     //  calculate betasquare values
       T = Tdat[iT] ;
       E = T + electron_mass_c2 ;
       b2small = T*(E+electron_mass_c2)/(E*E) ;
       T = Tdat[iT+1] ;
       E = T + electron_mass_c2 ;
       b2big = T*(E+electron_mass_c2)/(E*E) ;
       ratb2 = (beta2-b2small)/(b2big-b2small) ;
 
       corrfactor = tuning*(1.+cpar)/(1.+cpar*beta2) ;

     if(Charge < 0.)
     {
       c1 = celectron[iZ][iT] ;
       c2 = celectron[iZ+1][iT] ;
       cc1 = c1+ratZ*(c2-c1) ;

       c1 = celectron[iZ][iT+1] ;
       c2 = celectron[iZ+1][iT+1] ;
       cc2 = c1+ratZ*(c2-c1) ;

       corr = cc1+ratb2*(cc2-cc1) ;

       sigma /= corr ;

     }

     if(Charge > 0.)
     {
       c1 = cpositron[iZ][iT] ;
       c2 = cpositron[iZ+1][iT] ;
       cc1 = c1+ratZ*(c2-c1) ;

       c1 = cpositron[iZ][iT+1] ;
       c2 = cpositron[iZ+1][iT+1] ;
       cc2 = c1+ratZ*(c2-c1) ;

       corr = cc1+ratb2*(cc2-cc1) ;

       sigma /= corr ;
     }

     sigma *= corrfactor ;

    //  nucl. size correction for particles other than e+/e- only at present !!!!
    if((&aParticleType != G4Electron::Electron()) &&
       (&aParticleType != G4Positron::Positron())  )    
       sigma /= corrnuclsize ;

     return sigma ;
  } 


  G4VParticleChange* G4MultipleScattering::PostStepDoIt(
                                               const G4Track& trackData,
                                               const G4Step& stepData)
  {
    static G4double taulim=1.e-6 , randlim = 0.25*taulim*taulim  ;
    static const G4double tausmall = 5.e-5,taubig =50.,
          kappa = 2.5, kappapl1 = kappa+1., kappami1 = kappa-1. ;
    const G4DynamicParticle* aParticle ;
    G4double KineticEnergy,truestep,tau,prob,cth,sth,phi,
             dirx,diry,dirz,w,w1,etau,rmean,safetyminustolerance;
    G4double rand ;
    G4bool isOut;

    fParticleChange.Initialize(trackData) ;

    truestep = stepData.GetStepLength() ;

    aParticle = trackData.GetDynamicParticle() ;

    KineticEnergy = aParticle->GetKineticEnergy() ;

    if(stepFlag == 0)
    {
       fTransportMeanFreePath = (*theTransportMeanFreePathTable)
                               (materialIndex)->GetValue(KineticEnergy,isOut);
    } 

  //  change direction first ( scattering ) ..........................
    if(stepFlag == 0)
    {
      tau = truestep/fTransportMeanFreePath ;

     if(tau > taulim) 
      prob = exp(-tau)*(1.+scatteringparameter*tau) ;
     else
      prob=1.-(1.-scatteringparameter)*tau ;
    }
    else
    {
      tau = truestep/range ;
      if(tau < taulim)
        prob = exp(-(alpha1-1.)*tau)*(1.+scatteringparameter*tau) ;
      else if(tau < 1.) 
        prob = exp((alpha1-1.)*log(1.-tau))*(1.+scatteringparameter*tau) ;
      else
        prob = 0. ;
    }

    if(G4UniformRand()<prob)
    {
      if(tau<taulim)
      {
        rand = G4UniformRand() ;
        if(rand > randlim)
        {
          cth = 1.-scatteringparameter*tau*(1./sqrt(rand)-1.)  ;
          if(cth < -1.) cth = -1. ;
        }
        else 
          cth = -1. ; 
      }
      else
      {
        w = 1.+scatteringparameter*tau ;
        w1 = w-1. ;
        cth = w-w1*(w+1.)/sqrt(w1*w1+4.*w*G4UniformRand()) ;
      }
    }
    else
      cth = -1.+2.*G4UniformRand() ;
                                                                 
    sth = sqrt(1.-cth*cth) ;

    phi = twopi*G4UniformRand() ;

    dirx = sth*cos(phi) ;
    diry = sth*sin(phi) ;
    dirz = cth ;

    G4ParticleMomentum ParticleDirection = aParticle->GetMomentumDirection();

    G4ThreeVector newDirection(dirx,diry,dirz) ;


    newDirection.rotateUz(ParticleDirection) ;

    fParticleChange.SetMomentumChange(newDirection.x(),
                                      newDirection.y(),
                                      newDirection.z()) ;

    if(fLatDisplFlag)
    {
      //  compute mean lateral displacement ...............
      //  only for safety > tolerance !!!!!!!!!
      safetyminustolerance = stepData.GetPostStepPoint()->GetSafety();

      if(safetyminustolerance > 0.)
      {
        if(tau<tausmall)
          rmean = 5.*tau*tau*tau/12. ;
        else
        {
          if(tau<taubig)
            etau = exp(-tau) ;
          else
            etau = 0. ;
          rmean = -kappa*tau ;
          rmean = -exp(rmean)/(kappa*kappami1) ;
          rmean += tau-kappapl1/kappa+kappa*etau/kappami1 ;
        }

        if(rmean>0.)
          rmean = 2.*fTransportMeanFreePath*sqrt(rmean/3.) ;
        else
          rmean = 0. ;

        // for rmean > 0) only
        if(rmean>0.)
        {
          if(rmean>safetyminustolerance)
            rmean = safetyminustolerance ;

          //  sample direction of lateral displacement
          phi = twopi*G4UniformRand() ;

          dirx = cos(phi) ;
          diry = sin(phi) ;
          dirz = 0. ;

          G4ThreeVector latDirection(dirx,diry,dirz);

          latDirection.rotateUz(ParticleDirection) ;

          // compute new endpoint of the Step
	  G4ThreeVector newPosition= 
			 stepData.GetPostStepPoint()->GetPosition()+
					     rmean*latDirection;
    
	  G4Navigator *navigator= 
	   G4TransportationManager::GetTransportationManager()
				   ->GetNavigatorForTracking();
	  navigator->LocateGlobalPointWithinVolume( newPosition );

	  fParticleChange.SetPositionChange(newPosition) ;
        }
      }
    }


    return &fParticleChange ;
 
  }

  
void G4MultipleScattering::PrintInfoDefinition()
{
  G4String comments = " Tables of transport mean free paths.";
        comments += "\n          New model of MSC , computes the lateral \n";
        comments += "          displacement of the particle , too.";

  G4cout << G4endl << GetProcessName() << ":  " << comments
         << "\n        PhysicsTables from " << G4BestUnit(LowestKineticEnergy,
                                                  "Energy")
         << " to " << G4BestUnit(HighestKineticEnergy,"Energy")
         << " in " << TotBin << " bins. \n";
}

}//namespace GlastMS