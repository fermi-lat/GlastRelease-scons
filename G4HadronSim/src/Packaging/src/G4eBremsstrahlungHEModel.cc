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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eBremsstrahlungHEModel
//
// Author:        Andreas Schaelicke 
//                extention of standard G4eBremsstrahlungModel
//
// Creation date: 28.03.2008
//
// Modifications:
//  31.03.2008    Added modified LPM suppression function according to
//                Klein, Rev.Mod.Phys.71 (1999) 1501  eqs.(72)-(77)
//                (ratio of BH xsec with and without LPM)
//  28.04.2008    improved low energy behaviour in LPM suppression
//  13.05.2008    1st prototype version
//                 - included Density & LPM combination ala Ter-Mikaelian
//                 - suppression function more smooth and numerical stable
//
// Main References:
//  S.Klein,  Rev. Mod. Phys. 71 (1999) 1501.
//  T.Stanev et.al., Phys. Rev. D25 (1982) 1291.
//  M.L.Ter-Mikaelian, High-energy Electromagnetic Processes in Condensed Media, Wiley, 1972.
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eBremsstrahlungHEModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForLoss.hh"

#include "G4LossTableManager.hh"

#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eBremsstrahlungHEModel::G4eBremsstrahlungHEModel(const G4ParticleDefinition* p,
                                               const G4String& nam)
  : G4VEmModel(nam),
    particle(0),
    isElectron(true),
    probsup(1.0),
    //    MigdalConstant(classic_electr_radius*electron_Compton_length*electron_Compton_length/pi),
    MigdalConstant(classic_electr_radius*electron_Compton_length*electron_Compton_length*4.0*pi),
    LPMconstant(fine_structure_const*electron_mass_c2*electron_mass_c2/(4.*pi*hbarc)),
    theLPMflag(true),
    isInitialised(false)
{
  if(p) SetParticle(p);
  theGamma = G4Gamma::Gamma();
  minThreshold = 1.0*keV;
  highKinEnergy= 100.*TeV;
  lowKinEnergy = 1.0*keV;
  highEnergyTh = DBL_MAX;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eBremsstrahlungHEModel::~G4eBremsstrahlungHEModel()
{
  size_t n = partialSumSigma.size();
  if(n > 0) {
    for(size_t i=0; i<n; i++) {
      delete partialSumSigma[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungHEModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  if(p == G4Electron::Electron()) isElectron = true;
  else                            isElectron = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungHEModel::MinEnergyCut(const G4ParticleDefinition*,
                                              const G4MaterialCutsCouple*)
{
  return minThreshold;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eBremsstrahlungHEModel::Initialise(const G4ParticleDefinition* p,
                                        const G4DataVector& cuts)
{
  // *** update flags from losstablemanager *** 
  G4LossTableManager* man = G4LossTableManager::Instance(); 
  //SetEnergyThreshold(man->BremsstrahlungTh());
  //SetLPMflag(man->LPMFlag());

  if(p) SetParticle(p);
  highKinEnergy = HighEnergyLimit();
  lowKinEnergy  = LowEnergyLimit();
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();

  if(theCoupleTable) {
    G4int numOfCouples = theCoupleTable->GetTableSize();

    G4int nn = partialSumSigma.size();
    G4int nc = cuts.size();
    if(nn > 0) {
      for (G4int ii=0; ii<nn; ii++){
	G4DataVector* a=partialSumSigma[ii];
	if ( a )  delete a;
      }
      partialSumSigma.clear();
    }
    if(numOfCouples>0) {
      for (G4int i=0; i<numOfCouples; i++) {
        G4double cute   = DBL_MAX;
        if(i < nc) cute = cuts[i];
	const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
	const G4Material* material = couple->GetMaterial();
	G4DataVector* dv = ComputePartialSumSigma(material, 0.5*highKinEnergy,
                             std::min(cute, 0.25*highKinEnergy));
	partialSumSigma.push_back(dv);
      }
    }
  }
  if(isInitialised) return;

  if(pParticleChange)
    fParticleChange = reinterpret_cast<G4ParticleChangeForLoss*>(pParticleChange);
  else
    fParticleChange = new G4ParticleChangeForLoss();

  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungHEModel::ComputeDEDXPerVolume(
					     const G4Material* material,
                                             const G4ParticleDefinition* p,
                                                   G4double kineticEnergy,
                                                   G4double cutEnergy)
{
  if(!particle) SetParticle(p);
  if(kineticEnergy < lowKinEnergy) return 0.0;

  const G4double thigh = 100.*GeV;

  G4double cut = std::min(cutEnergy, kineticEnergy);

  G4double rate, loss;
  const G4double factorHigh = 36./(1450.*GeV);
  const G4double coef1 = -0.5;
  const G4double coef2 = 2./9.;

  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomicNumDensityVector = material->GetAtomicNumDensityVector();

  G4double totalEnergy = kineticEnergy + electron_mass_c2;
  G4double dedx = 0.0;

  //  loop for elements in the material
  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    G4double Z     = (*theElementVector)[i]->GetZ();
    G4double natom = theAtomicNumDensityVector[i];

    // loss for MinKinEnergy<KineticEnergy<=100 GeV
    if (kineticEnergy <= thigh) {

      //      x = log(totalEnergy/electron_mass_c2);
      loss = ComputeBremLoss(Z, kineticEnergy, cut) ;
      if (!isElectron) loss *= PositronCorrFactorLoss(Z, kineticEnergy, cut);

    // extrapolation for KineticEnergy>100 GeV
    } else {

      //      G4double xhigh = log(thigh/electron_mass_c2);
      G4double cuthigh = thigh*0.5;

      if (cut < thigh) {

        loss = ComputeBremLoss(Z, thigh, cut) ;
        if (!isElectron) loss *= PositronCorrFactorLoss(Z, thigh, cut) ;
        rate = cut/totalEnergy;
        loss *= (1. + coef1*rate + coef2*rate*rate);
        rate = cut/thigh;
        loss /= (1.+coef1*rate+coef2*rate*rate);

      } else {

        loss = ComputeBremLoss(Z, thigh, cuthigh) ;
        if (!isElectron) loss *= PositronCorrFactorLoss(Z, thigh, cuthigh) ;
        rate = cut/totalEnergy;
        loss *= (1. + coef1*rate + coef2*rate*rate);
        loss *= cut*factorHigh;
      }
    }
    loss *= natom;

    G4double kp2 = MigdalConstant*totalEnergy*totalEnergy
                 * (material->GetElectronDensity()) ;

    // now compute the correction due to the supression(s)
    G4double kmin = 1.*eV;
    G4double kmax = cut;

    if (kmax > kmin) {

      G4double floss = 0.;
      G4int nmax = 100;

      G4double vmin=log(kmin);
      G4double vmax=log(kmax) ;
      G4int nn = (G4int)(nmax*(vmax-vmin)/(log(highKinEnergy)-vmin)) ;
      G4double u,fac,c,v,dv ;
      if(nn > 0) {

        dv = (vmax-vmin)/nn ;
        v  = vmin-dv ;

        for(G4int n=0; n<=nn; n++) {

          v += dv;
          u = exp(v);
          fac = u*SupressionFunction(material,kineticEnergy,u);
	  fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;
          if ((n==0)||(n==nn)) c=0.5;
          else                 c=1. ;
          fac   *= c ;
          floss += fac ;
        }
        floss *=dv/(kmax-kmin);

      } else {
        floss = 1.;
      }
      if(floss > 1.) floss = 1.;
      // correct the loss
      loss *= floss;
    }
    dedx += loss;
  }
  if(dedx < 0.) dedx = 0.;
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungHEModel::ComputeBremLoss(G4double Z, G4double T,
                                                 G4double Cut)

// compute loss due to soft brems
{
  static const G4double beta=1.0, ksi=2.0;
  static const G4double clossh = 0.254 , closslow = 1./3. , alosslow = 1. ;
  static const G4double Tlim= 10.*MeV ;

  static const G4double xlim = 1.2 ;
  static const G4int NZ = 8 ;
  static const G4int Nloss = 11 ;
  static const G4double ZZ[NZ] =
        {2.,4.,6.,14.,26.,50.,82.,92.};
  static const G4double coefloss[NZ][Nloss] = {
  // Z=2
 { 0.98916,        0.47564,        -0.2505,       -0.45186,        0.14462,
   0.21307,      -0.013738,      -0.045689,     -0.0042914,      0.0034429,
   0.00064189},

  // Z=4
 { 1.0626,        0.37662,       -0.23646,       -0.45188,        0.14295,
   0.22906,      -0.011041,      -0.051398,     -0.0055123,      0.0039919,
   0.00078003},
  // Z=6
 { 1.0954,          0.315,       -0.24011,       -0.43849,        0.15017,
   0.23001,      -0.012846,      -0.052555,     -0.0055114,      0.0041283,
   0.00080318},

  // Z=14
 { 1.1649,        0.18976,       -0.24972,       -0.30124,         0.1555,
   0.13565,      -0.024765,      -0.027047,    -0.00059821,      0.0019373,
   0.00027647},

  // Z=26
 { 1.2261,        0.14272,       -0.25672,       -0.28407,        0.13874,
   0.13586,      -0.020562,      -0.026722,    -0.00089557,      0.0018665,
   0.00026981},

  // Z=50
 { 1.3147,       0.020049,       -0.35543,       -0.13927,        0.17666,
   0.073746,      -0.036076,      -0.013407,      0.0025727,     0.00084005,
  -1.4082e-05},

  // Z=82
 { 1.3986,       -0.10586,       -0.49187,     -0.0048846,        0.23621,
   0.031652,      -0.052938,     -0.0076639,      0.0048181,     0.00056486,
  -0.00011995},

  // Z=92
 { 1.4217,         -0.116,       -0.55497,      -0.044075,        0.27506,
   0.081364,      -0.058143,      -0.023402,      0.0031322,      0.0020201,
   0.00017519}

    } ;
  static G4double aaa = 0.414;
  static G4double bbb = 0.345;
  static G4double ccc = 0.460;

  G4int iz = 0;
  G4double delz = 1.e6;
  for (G4int ii=0; ii<NZ; ii++)
    {
      G4double dz = std::abs(Z-ZZ[ii]);
      if(dz < delz)  {
        iz = ii;
        delz = dz;
      }
    }

  G4double xx = log10(T/MeV);
  G4double fl = 1.;

  if (xx <= xlim)
    {
      xx /= xlim;
      G4double yy = 1.0;
      fl = 0.0;
      for (G4int j=0; j<Nloss; j++) {
	fl += yy+coefloss[iz][j];
        yy *= xx;
      }
      if (fl < 0.00001) fl = 0.00001;
      else if (fl > 1.0) fl = 1.0;
    }

  G4double loss;
  G4double E = T+electron_mass_c2 ;

  loss = Z*(Z+ksi)*E*E/(T+E)*exp(beta*log(Cut/T))*(2.-clossh*exp(log(Z)/4.));
  if (T <= Tlim) loss /= exp(closslow*log(Tlim/T));
  if( T <= Cut)  loss *= exp(alosslow*log(T/Cut));
  //  correction
  loss *= (aaa+bbb*T/Tlim)/(1.+ccc*T/Tlim);
  loss *= fl;
  loss /= Avogadro;

  return loss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungHEModel::PositronCorrFactorLoss(G4double Z,
                                 G4double kineticEnergy, G4double cut)

//calculates the correction factor for the energy loss due to bremsstrahlung for positrons
//the same correction is in the (discrete) bremsstrahlung

{
  static const G4double K = 132.9416*eV ;
  static const G4double a1=4.15e-1, a3=2.10e-3, a5=54.0e-5 ;

  G4double x   = log(kineticEnergy/(K*Z*Z)), x2 = x*x, x3 = x2*x;
  G4double eta = 0.5+atan(a1*x+a3*x3+a5*x3*x2)/pi;
  G4double e0  = cut/kineticEnergy;

  G4double factor = 0.0;
  if (e0 < 1.0) {
    factor=log(1.-e0)/eta;
    factor=exp(factor);
  }
  factor = eta*(1.-factor)/e0;

  return factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungHEModel::CrossSectionPerVolume(
					      const G4Material* material,
                                              const G4ParticleDefinition* p,
                                                    G4double kineticEnergy,
                                                    G4double cutEnergy,
                                                    G4double maxEnergy)
{
  if(!particle) SetParticle(p);
  G4double cross = 0.0;
  G4double tmax = min(maxEnergy, kineticEnergy);
  G4double cut  = max(cutEnergy, minThreshold);
  if(cut >= tmax) return cross;

  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();
  G4double dum=0.;
  
  for (size_t i=0; i<material->GetNumberOfElements(); i++) {

    cross += theAtomNumDensityVector[i] * ComputeCrossSectionPerAtom(p,
                      kineticEnergy, (*theElementVector)[i]->GetZ(), dum, cut);
    if (tmax < kineticEnergy) {
      cross -= theAtomNumDensityVector[i] * ComputeCrossSectionPerAtom(p,
                     kineticEnergy, (*theElementVector)[i]->GetZ(), dum, tmax);
    }
  }

  // now compute the correction due to the supression(s)

  G4double kmax = tmax;
  G4double kmin = cut;

  G4double totalEnergy = kineticEnergy+electron_mass_c2 ;
  G4double kp2 = MigdalConstant*totalEnergy*totalEnergy
                                             *(material->GetElectronDensity());

  G4double fsig = 0.;
  G4int nmax = 100;
  G4double vmin=log(kmin);
  G4double vmax=log(kmax) ;
  G4int nn = (G4int)(nmax*(vmax-vmin)/(log(highKinEnergy)-vmin));
  G4double u,fac,c,v,dv,y ;
  if(nn > 0) {

      dv = (vmax-vmin)/nn ;
      v  = vmin-dv ;
      for(G4int n=0; n<=nn; n++) {

        v += dv;  
        u = exp(v);              
        fac = SupressionFunction(material, kineticEnergy, u);
        y = u/kmax;
        fac *= (4.-4.*y+3.*y*y)/3.;
        fac *= probsup*(u*u/(u*u+kp2))+1.-probsup;

        if ((n==0)||(n==nn)) c=0.5;
        else                 c=1. ;

        fac  *= c;
        fsig += fac;
      }
      y = kmin/kmax ;
      fsig *=dv/(-4.*log(y)/3.-4.*(1.-y)/3.+0.5*(1.-y*y));

  } else {

      fsig = 1.;
  }
  if (fsig > 1.) fsig = 1.;

  // correct the cross section
  cross *= fsig;

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eBremsstrahlungHEModel::ComputeCrossSectionPerAtom(
                                              const G4ParticleDefinition*,
				                    G4double kineticEnergy, 
                                                    G4double Z,   G4double,
						    G4double cut, G4double)
 
// Calculates the cross section per atom in GEANT4 internal units.
//
 
{
  G4double cross = 0.0 ;
  if ( kineticEnergy < 1*keV || kineticEnergy < cut) return cross;

  static const G4double ksi=2.0, alfa=1.00;
  static const G4double csigh = 0.127, csiglow = 0.25, asiglow = 0.020*MeV ;
  static const G4double Tlim = 10.*MeV ;

  static const G4double xlim = 1.2 ;
  static const G4int NZ = 8 ;
  static const G4int Nsig = 11 ;
  static const G4double ZZ[NZ] =
        {2.,4.,6.,14.,26.,50.,82.,92.} ;
  static const G4double coefsig[NZ][Nsig] = {
  // Z=2
  { 0.4638,        0.37748,        0.32249,      -0.060362,      -0.065004,
   -0.033457,      -0.004583,       0.011954,      0.0030404,     -0.0010077,
   -0.00028131},

  // Z=4
  { 0.50008,        0.33483,        0.34364,      -0.086262,      -0.055361,
   -0.028168,     -0.0056172,       0.011129,      0.0027528,    -0.00092265,
   -0.00024348},

  // Z=6
  { 0.51587,        0.31095,        0.34996,       -0.11623,      -0.056167,
   -0.0087154,     0.00053943,      0.0054092,     0.00077685,    -0.00039635,
   -6.7818e-05},

  // Z=14
  { 0.55058,        0.25629,        0.35854,      -0.080656,      -0.054308,
   -0.049933,    -0.00064246,       0.016597,      0.0021789,      -0.001327,
   -0.00025983},

  // Z=26
  { 0.5791,        0.26152,        0.38953,       -0.17104,      -0.099172,
    0.024596,       0.023718,     -0.0039205,     -0.0036658,     0.00041749,
    0.00023408},

  // Z=50
  { 0.62085,        0.27045,        0.39073,       -0.37916,       -0.18878,
    0.23905,       0.095028,      -0.068744,      -0.023809,      0.0062408,
    0.0020407},

  // Z=82
  { 0.66053,        0.24513,        0.35404,       -0.47275,       -0.22837,
    0.35647,        0.13203,        -0.1049,      -0.034851,      0.0095046,
    0.0030535},

  // Z=92
  { 0.67143,        0.23079,        0.32256,       -0.46248,       -0.20013,
    0.3506,        0.11779,        -0.1024,      -0.032013,      0.0092279,
    0.0028592}

    } ;

  G4int iz = 0 ;
  G4double delz = 1.e6 ;
  for (G4int ii=0; ii<NZ; ii++)
  {
    G4double absdelz = std::abs(Z-ZZ[ii]); 
    if(absdelz < delz)
    {
      iz = ii ;
      delz = absdelz;
    }
  }

  G4double xx = log10(kineticEnergy/MeV) ;
  G4double fs = 1. ;

  if (xx <= xlim) {

    fs = coefsig[iz][Nsig-1] ;
    for (G4int j=Nsig-2; j>=0; j--) {

      fs = fs*xx+coefsig[iz][j] ;
    }
    if(fs < 0.) fs = 0.;
  }

  cross = Z*(Z+ksi)*(1.-csigh*exp(log(Z)/4.))*pow(log(kineticEnergy/cut),alfa);

  if (kineticEnergy <= Tlim)
     cross *= exp(csiglow*log(Tlim/kineticEnergy))
              *(1.+asiglow/(sqrt(Z)*kineticEnergy));

  if (!isElectron)
     cross *= PositronCorrFactorSigma(Z, kineticEnergy, cut);

  cross *= fs/Avogadro ;

  if (cross < 0.) cross = 0.;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
G4double G4eBremsstrahlungHEModel::PositronCorrFactorSigma( G4double Z,
                                 G4double kineticEnergy, G4double cut)

//Calculates the correction factor for the total cross section of the positron
//  bremsstrahl.
// Eta is the ratio of positron to electron energy loss by bremstrahlung. 
// A parametrized formula from L. Urban is used to estimate eta. It is a fit to
// the results of L. Kim & al: Phys Rev. A33,3002 (1986)

{
  static const G4double K = 132.9416*eV;
  static const G4double a1 = 4.15e-1, a3 = 2.10e-3, a5 = 54.0e-5;

  G4double x    = log(kineticEnergy/(K*Z*Z));
  G4double x2 = x*x;
  G4double x3 = x2*x;
  G4double eta  = 0.5 + atan(a1*x + a3*x3 + a5*x3*x2)/pi ;
  G4double alfa = (1. - eta)/eta;
  return eta*pow((1. - cut/kineticEnergy), alfa);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DataVector* G4eBremsstrahlungHEModel::ComputePartialSumSigma(
              const G4Material* material,
                    G4double kineticEnergy,
                    G4double cut)

// Build the table of cross section per element. 
//The table is built for MATERIALS.
// This table is used by DoIt to select randomly an element in the material.
{
  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();
  G4double dum = 0.;
  
  G4DataVector* dv = new G4DataVector();

  G4double cross = 0.0;

  for (G4int i=0; i<nElements; i++ ) {

    cross += theAtomNumDensityVector[i] * ComputeCrossSectionPerAtom( particle,
                       kineticEnergy, (*theElementVector)[i]->GetZ(), dum, cut);
    dv->push_back(cross);
  }
  return dv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4DynamicParticle*>* G4eBremsstrahlungHEModel::SampleSecondaries(const G4MaterialCutsCouple* couple,
					       const G4DynamicParticle* dp,
					       G4double tmin,
					       G4double maxEnergy)
// The emitted gamma energy is sampled using a parametrized formula 
// from L. Urban.
// This parametrization is derived from :
//    cross-section values of Seltzer and Berger for electron energies
//    1 keV - 10 GeV,
//    screened Bethe Heilter differential cross section above 10 GeV,
//    Migdal corrections in both case.
//  Seltzer & Berger: Nim B 12:95 (1985)
//  Nelson, Hirayama & Rogers: Technical report 265 SLAC (1985)
//  Migdal: Phys Rev 103:1811 (1956); Messel & Crawford: Pergamon Press (1970)
//
// A modified version of the random number techniques of Butcher&Messel is used
//    (Nuc Phys 20(1960),15).
{
  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double tmax = min(maxEnergy, kineticEnergy);
  if(tmin >= tmax) return 0;

//
// GEANT4 internal units.
//
  static const G4double
     ah10 = 4.67733E+00, ah11 =-6.19012E-01, ah12 = 2.02225E-02,
     ah20 =-7.34101E+00, ah21 = 1.00462E+00, ah22 =-3.20985E-02,
     ah30 = 2.93119E+00, ah31 =-4.03761E-01, ah32 = 1.25153E-02;

  static const G4double
     bh10 = 4.23071E+00, bh11 =-6.10995E-01, bh12 = 1.95531E-02,
     bh20 =-7.12527E+00, bh21 = 9.69160E-01, bh22 =-2.74255E-02,
     bh30 = 2.69925E+00, bh31 =-3.63283E-01, bh32 = 9.55316E-03;

  static const G4double
     al00 =-2.05398E+00, al01 = 2.38815E-02, al02 = 5.25483E-04,
     al10 =-7.69748E-02, al11 =-6.91499E-02, al12 = 2.22453E-03,
     al20 = 4.06463E-02, al21 =-1.01281E-02, al22 = 3.40919E-04;

  static const G4double
     bl00 = 1.04133E+00, bl01 =-9.43291E-03, bl02 =-4.54758E-04,
     bl10 = 1.19253E-01, bl11 = 4.07467E-02, bl12 =-1.30718E-03,
     bl20 =-1.59391E-02, bl21 = 7.27752E-03, bl22 =-1.94405E-04;

  static const G4double tlow = 1.*MeV;

  G4double gammaEnergy;
  G4bool LPMOK = false;
  const G4Material* material = couple->GetMaterial();

  // select randomly one element constituing the material
  const G4Element* anElement = SelectRandomAtom(couple);

  // Extract Z factors for this Element
  G4double lnZ = 3.*(anElement->GetIonisation()->GetlogZ3());
  G4double FZ = lnZ* (4.- 0.55*lnZ);
  G4double ZZ = anElement->GetIonisation()->GetZZ3();

  // limits of the energy sampling
  G4double totalEnergy = kineticEnergy + electron_mass_c2;
  G4ThreeVector direction = dp->GetMomentumDirection();
  G4double xmin     = tmin/kineticEnergy;
  G4double xmax     = tmax/kineticEnergy;
  G4double kappa    = 0.0;
  if(xmax >= 1.) xmax = 1.;
  else     kappa    = log(xmax)/log(xmin);
  G4double epsilmin = tmin/totalEnergy;
  G4double epsilmax = tmax/totalEnergy;

  // Migdal factor
  G4double MigdalFactor = (material->GetElectronDensity())*MigdalConstant
                        / (epsilmax*epsilmax);

  G4double x, epsil, greject, migdal, grejmax, q;
  G4double U  = log(kineticEnergy/electron_mass_c2);
  G4double U2 = U*U;

  // precalculated parameters
  G4double ah, bh;
  G4double screenfac = 0.0;

  if (kineticEnergy > tlow) {
       
    G4double ah1 = ah10 + ZZ* (ah11 + ZZ* ah12);
    G4double ah2 = ah20 + ZZ* (ah21 + ZZ* ah22);
    G4double ah3 = ah30 + ZZ* (ah31 + ZZ* ah32);

    G4double bh1 = bh10 + ZZ* (bh11 + ZZ* bh12);
    G4double bh2 = bh20 + ZZ* (bh21 + ZZ* bh22);
    G4double bh3 = bh30 + ZZ* (bh31 + ZZ* bh32);

    ah = 1.   + (ah1*U2 + ah2*U + ah3) / (U2*U);
    bh = 0.75 + (bh1*U2 + bh2*U + bh3) / (U2*U);

    // limit of the screening variable
    screenfac =
      136.*electron_mass_c2/((anElement->GetIonisation()->GetZ3())*totalEnergy);
    G4double screenmin = screenfac*epsilmin/(1.-epsilmin);

    // Compute the maximum of the rejection function
    G4double F1 = max(ScreenFunction1(screenmin) - FZ ,0.);
    G4double F2 = max(ScreenFunction2(screenmin) - FZ ,0.);
    grejmax = (F1 - epsilmin* (F1*ah - bh*epsilmin*F2))/(42.392 - FZ);

  } else {  

    G4double al0 = al00 + ZZ* (al01 + ZZ* al02);
    G4double al1 = al10 + ZZ* (al11 + ZZ* al12);
    G4double al2 = al20 + ZZ* (al21 + ZZ* al22);
 
    G4double bl0 = bl00 + ZZ* (bl01 + ZZ* bl02);
    G4double bl1 = bl10 + ZZ* (bl11 + ZZ* bl12);
    G4double bl2 = bl20 + ZZ* (bl21 + ZZ* bl22);
 
    ah = al0 + al1*U + al2*U2;
    bh = bl0 + bl1*U + bl2*U2;

    // Compute the maximum of the rejection function
    grejmax = max(1. + xmin* (ah + bh*xmin), 1.+ah+bh);
    G4double xm = -ah/(2.*bh);
    if ( xmin < xm && xm < xmax) grejmax = max(grejmax, 1.+ xm* (ah + bh*xm));
  }

  //
  //  sample the energy rate of the emitted gamma for e- kin energy > 1 MeV
  //
 
  do {
    if (kineticEnergy > tlow) {
      do {
        q = G4UniformRand();
        x = pow(xmin, q + kappa*(1.0 - q));
        epsil = x*kineticEnergy/totalEnergy;
        G4double screenvar = screenfac*epsil/(1.0-epsil);
        G4double F1 = max(ScreenFunction1(screenvar) - FZ ,0.);
        G4double F2 = max(ScreenFunction2(screenvar) - FZ ,0.);
        migdal = (1. + MigdalFactor)/(1. + MigdalFactor/(x*x));
        greject = migdal*(F1 - epsil* (ah*F1 - bh*epsil*F2))/(42.392 - FZ);
	/*
	if ( greject > grejmax ) {
            G4cout << "### G4eBremsstrahlungHEModel Warning: Majoranta exceeded! "
                   << greject << " > " << grejmax
                   << " x= " << x 
		   << " e= " << kineticEnergy
                   << G4endl;
	}
	*/	
      } while( greject < G4UniformRand()*grejmax );

    } else {  

      do {
        q = G4UniformRand();
        x = pow(xmin, q + kappa*(1.0 - q));
        migdal = (1. + MigdalFactor)/(1. + MigdalFactor/(x*x));  
        greject = migdal*(1. + x* (ah + bh*x));
	/*
	if ( greject > grejmax ) {
	  G4cout << "### G4eBremsstrahlungHEModel Warning: Majoranta exceeded! " 
                 << greject << " > " << grejmax 
                 << " x= " << x 
		 << " e= " << kineticEnergy
                 << G4endl;
	}
	*/
      } while( greject < G4UniformRand()*grejmax );
    }
    /*
    if(x > 0.999) {
      G4cout << "### G4eBremsstrahlungHEModel Warning: e= " << kineticEnergy
	     << " tlow= " << tlow
	     << " x= " << x
	     << " greject= " << greject 
	     << " grejmax= " << grejmax
	     << " migdal= " << migdal
	     << G4endl; 
      //      if(x >= 1.0) G4Exception("X=1");
    }
    */
    gammaEnergy = x*kineticEnergy; 

    if (theLPMflag) {
     // take into account the supression due to the LPM effect
      if (G4UniformRand() <= SupressionFunction(material,kineticEnergy,
                                                           gammaEnergy))
        LPMOK = true;
      }
    else LPMOK = true;

  } while (!LPMOK);

  //
  //  angles of the emitted gamma. ( Z - axis along the parent particle)
  //
  //  universal distribution suggested by L. Urban 
  // (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  G4double u;
  const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  if (9./(9.+d) > G4UniformRand()) u = - log(G4UniformRand()*G4UniformRand())/a1;
     else                          u = - log(G4UniformRand()*G4UniformRand())/a2;

  G4double theta = u*electron_mass_c2/totalEnergy;

  G4double sint = sin(theta);

  G4double phi = twopi * G4UniformRand() ;

  G4ThreeVector gammaDirection(sint*cos(phi),sint*sin(phi), cos(theta));
  gammaDirection.rotateUz(direction);

  // create G4DynamicParticle object for the Gamma

  std::vector<G4DynamicParticle*>* vdp = new std::vector<G4DynamicParticle*>;


  G4DynamicParticle* g = new G4DynamicParticle(theGamma,gammaDirection,
                                                        gammaEnergy);
  vdp->push_back(g);
  
  G4double totMomentum = sqrt(kineticEnergy*(totalEnergy + electron_mass_c2));
  G4ThreeVector dir = totMomentum*direction - gammaEnergy*gammaDirection;
  direction = dir.unit();

  // energy of primary
  G4double finalE = kineticEnergy - gammaEnergy;

  // stop tracking and create new secondary instead of primary
  if(gammaEnergy > highEnergyTh) {
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    G4DynamicParticle* el = 
      new G4DynamicParticle(const_cast<G4ParticleDefinition*>(particle),
			    direction, finalE);
    vdp->push_back(el);

    // continue tracking
  } else {
    fParticleChange->SetProposedMomentumDirection(direction);
    fParticleChange->SetProposedKineticEnergy(finalE);
  }

  return vdp;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4eBremsstrahlungHEModel::SelectRandomAtom(
           const G4MaterialCutsCouple* couple) 
{
  // select randomly 1 element within the material

  const G4Material* material = couple->GetMaterial();
  G4int nElements = material->GetNumberOfElements();
  const G4ElementVector* theElementVector = material->GetElementVector();

  const G4Element* elm = 0;

  if(1 < nElements) {

    G4DataVector* dv = partialSumSigma[couple->GetIndex()];
    G4double rval = G4UniformRand()*((*dv)[nElements-1]);

    for (G4int i=0; i<nElements; i++) {
      if (rval <= (*dv)[i]) elm = (*theElementVector)[i];
    }
    if(!elm) {
      G4cout << "G4eBremsstrahlungHEModel::SelectRandomAtom: Warning -"
	     << " no elements found in "
	     << material->GetName()
	     << G4endl;
      elm = (*theElementVector)[0];
    }
  } else elm = (*theElementVector)[0];
 
  //  SetCurrentElement(elm);
  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4eBremsstrahlungHEModel::SupressionFunction(const G4Material* material,
                                 G4double kineticEnergy, G4double gammaEnergy)
{

  G4double Eele = kineticEnergy+electron_mass_c2 ;
  G4double k = gammaEnergy;

  G4double Elpm = LPMconstant*(material->GetRadlen()) ;

  // WARNING tested for single element materials only!
  //  G4double Z = material->GetZ(); 


  // temprary fix

  G4double Z = material->GetElectronDensity()/material->GetTotNbOfAtomsPerVolume();


  G4double supr = 1.0;

  if (theLPMflag) {
    G4double y = k/Eele;
    G4double y2 = y*y;
    G4double yone2 = 2.*(1+sqr(1.-y));

    // numerical safety factor
    // WARNING may need modification for very high energy electrons 
    if (k>0.9999*Eele) {
      k=0.9999*Eele;
    }

    // WARNING needs Eele and Z dependence!
    if (Elpm<=200.*TeV) {
      //   G4cout<<"Elpm="<<Elpm<<endl;

      // *** calculate lpm variable s & sprime ***
      // Klein eqs. (78) & (79)
      G4double sprime = sqrt(0.125*k*Elpm/(Eele*(Eele-k)));
      G4double s1 = pow(Z,2./3.)/(184.*184.);

      G4double h  = log(sprime)/log(sqrt(2.)*s1);
      G4double xiprime = 2;
      if (sprime>1) 
	xiprime = 1.;
      else if (sprime>sqrt(2.)*s1) 
	xiprime = 1+h-0.08*(1-h)*(1-sqr(1-h))/log(sqrt(2.)*s1);
      
      G4double s = sprime/sqrt(xiprime); 

      // *** merging with density effect ***
      // using Ter-Mikaelian eq. (20.9)
      G4double k2 = gammaEnergy*gammaEnergy ;
      G4double Eele2 = Eele*Eele;
      G4double eDensity = material->GetElectronDensity();
      G4double kp2=MigdalConstant*Eele2*eDensity;
      s = s * (1 + (kp2/k2) );

      // *** recalculate Xi ***
      // Klein eq. (75)
      G4double xi = 1;
      if (s<s1) xi = 2;
      if ( (s1<s) && (s<=1) ) xi = 1 + log(s)/log(s1);

      // Klein eqs. (77)
      G4double s2=s*s;
      G4double s3=s*s2;
      G4double s4=s2*s2;

      G4double phi=0.,G=0.;
      if (s<0.1) {
	// high suppression limit
	phi = 6.*s - 18.84955592153876*s2 + 39.47841760435743*s3 
	  - 57.69873135166053*s4;
	G = 37.69911184307752*s2 - 236.8705056261446*s3 + 807.7822389*s4;
      }
      else if (s<1.9516) {
	// intermediate suppression
	// using eq.77 approxim. valid s<2.      
	phi = 1-exp(-6*s*(1+(3-pi)*s)
		    +s3/(0.623+0.795*s+0.658*s2));
	if (s<0.415827397755) {
	  // using eq.77 approxim. valid 0.07<s<2
	  G4double psi = 1-exp(-4*s-8*s2/(1+3.936*s+4.97*s2-0.05*s3+7.50*s4));
	  G = 3*psi-2*phi;
	}
	else {
	  // using alternative parametrisiation
          G4double par[] = {-0.16072300849123999, 3.7550300067531581, 
			    -1.7981383069010097, 0.67282686077812381, -0.1207722909879257};
	  G4double pre =  par[0] + s*(par[1] + s*(par[2] + s*(par[3] + s*par[4])));
	  G = tanh(pre);
	}
      }
      else {
	// low suppression limit valid s>2.
	phi = 1. - 0.0119048/s4;
	G = 1. - 0.0230655/s4;
      }

      // *** make sure suppression is smaller than 1 ***
      // *** caused by Migdal approximation in xi    ***
      if (xiprime*phi>1. || s>0.57)  xiprime=1./phi;


      // *** calculate suppression factor ***
      G4double mainLPM = xiprime*(y2 * G + yone2*phi);
      G4double main = (y2 + yone2);
      supr = mainLPM/main;
    }
    
  } 
  return supr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


