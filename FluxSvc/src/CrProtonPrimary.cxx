/**
 * CrProtonPrimary:
 * The primary particle generator for primary cosmic-ray protons.
 *
 * First written by M. Ozaki(ISAS) and developed by K. Hirano(Hiroshima Univ.)
 * 2001-04-26 modified by T. Mizuno(Hiroshima Univ.)
 *
 */

#include <cmath>

// Geant4
//#include "Randomize.hh"		// in source/global/HEPRandom/include/
//#include "SystemOfUnits.h"	// in source/global/management/include/
// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrProtonPrimary.h"


typedef  double G4double;


// private function definitions.
namespace {
  const G4double pi = 3.14159265358979323846264339;
  // rest energy of proton in units of GeV
  const G4double restE = 0.938;
  // lower and higher energy limit of primary proton in units of GeV
  const G4double lowE_primary  = 1.0;
  const G4double highE_primary = 100.0;
 
  // The following value were computed from the model in advance.
  // This should be changed when the COR or force-field potential changes.
  // Currently it is computed for 4.44 [GV] and 1100 [MV], respectively.
  // The integral should be done from lowE_* to highE_*.
  // FIX ME!!

  // integrated downward flux in units of [m**-2 s**-1 sr**-1]
  const G4double INTEGRAL_primary = 329.67;



  // gives v/c as a function of kinetic Energy
  inline G4double beta(G4double E /* GeV */){
    return sqrt(1 - pow(E/restE+1, -2));
  }

  // gives rigidity (GV) as a function of kinetic Energy
  inline G4double rigidity(G4double E /* GeV */){
    return sqrt(pow(E + restE, 2) - pow(restE, 2));
  }

  // gives kinetic energy (GeV) as a function of rigidity
  inline G4double energy(G4double rigidity /* GV */){
    return sqrt(pow(rigidity, 2) + pow(restE, 2)) - restE;
  }

  //============================================================
  /**
   *  Generate a random distribution of primary cosmic ray protons
   *  j(E) = mod_spec(E, phi) * beta(E) * geomag_cut(E, CutOff)
   *    mod_spec(E, phi) = org_spec(E+phi*1e-3) * 
   *     ((E+restE)**2 - restE**2)/((E+restE+phi*1e-3)**2-restE**2)
   *    org_spec(E) = A * rigidity(E)**-a
   *      a = 2.79
   *    rigidity(E) = sqrt((E+restE)**2 - restE**2)
   *    beta(E) = sqrt(1 - (E/restE+1)**-2)
   *    geomag_cut(E, CutOff) = 1/(1 + (rigidity(E)/CutOff)**-12.0)
   *      CutOff = 4.44 for Theta_M = 0.735 and altitude = 50km (balloon experiment)
   *      phi = 540, 1100 [MV] for Solar minimum, maximum
   *    E: [GeV]
   *    j: [m**-2 s**-1 sr**-1 MeV**-1 ]
   *
   *  References:
   *  org_spec: AMS data (Alcaraz et al. 2000, Phys. Letter B 472, 215)
   *  geomag_cut formula: an eyeball fitting function of the ratio of BESS98 and AMS98 proton spectra.
   *  CutOff: calculated as (Rc/GV) = 14.9 * (1+h/R)^-2 * (cos(theta_M))^4,
   *          where h gives altitude from earth surface, R means an radius of earth,
   *          and theta_M is geomagnetic lattitude. See references below.
   *          "Handbook of space astronomy and astrophysics" 2nd edition, p225 
   *            (Zombeck, 1990, Cambridge University Press)
   *          "High Energy Astrophysics" 2nd edition, p325-330
   *            (M. S. Longair, 1992, Cambridge University Press)
   *  mod_spec: Gleeson, L. J. and Axford, W. I. 1968, ApJ, 154, 1011-1026 (Eq. 11)
   */

  const G4double A_primary = 16.9; // normalization of incident spectrum
  const G4double a_primary = 2.79; // differential spectral index

  // gives geomagnetic cutoff factor of primary cosmic rays
  // as a function of kinetic energy E(GeV) 
  // and geomagnetic lattitude theta_M(rad)
  inline G4double geomag_cut(G4double E, G4double cor /* GV */){
    return 1./(1 + pow(rigidity(E)/cor, -12.0));
  }

  // Unmodulated cosmic ray spectrum outside the Solar system
  inline G4double org_spec(G4double E /* GeV */)
  {
    return A_primary * pow(rigidity(E), -a_primary);
  }

  // Force-field approximation of the Solar modulation
  inline G4double mod_spec(G4double E /* GeV */, G4double phi /* MV */){
    return org_spec(E + phi*1e-3) * (pow(E+restE, 2) - pow(restE, 2))
      / (pow(E+restE+phi*1e-3,2) - pow(restE,2));
  }

  // The property function that the primary component should obey
  inline G4double primaryCRspec(G4double E /* GeV */, G4double cor /* GV*/, G4double phi /* MV */){
    return mod_spec(E, phi) * geomag_cut(E, cor) * beta(E);
  }

  // envelope function in lower energy
  // If we express rigidity as R, energy as E, cutoff rigidity as Rc,
  // and cutoff energy as Ec, geomagnetic cutoff function
  // is enveloped as follows.
  //=============================================================
  // 1/(1+(R/Rc)^-12.0) < (R/Rc)^12.0 < C * E^(a_envelope) / Rc^12.0
  //=============================================================
  // Here, C and a_envelope are determined das
  //  R(lowE_primary)^12.0 = C*(lowE_primary)^(a_envelope)
  // and 
  //  (Rc/Rc)^12.0 = C * Ec^(a_envelope) / Rc^12.0
  inline G4double primaryCRenvelope1(G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */){
    G4double a_envelope = 
      12.0*log(cor/rigidity(lowE_primary))/log(energy(cor));
    G4double A_envelope = 
      pow(rigidity(lowE_primary),12.0)/pow(lowE_primary,a_envelope);
    return A_primary * A_envelope * pow(cor, -12.0) * 
      pow(E, a_envelope - a_primary);
  }

  // integral of envelope function in lower energy
  inline G4double primaryCRenvelope1_integral
  (G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */){
    G4double a_envelope = 
      12.0*log(cor/rigidity(lowE_primary))/log(energy(cor));
    G4double A_envelope = 
      pow(rigidity(lowE_primary),12.0)/pow(lowE_primary,a_envelope);
    return A_primary * A_envelope * pow(cor, -12.0) * 
      1./(a_envelope - a_primary + 1) * pow(E, a_envelope - a_primary + 1);
  }

  // inverse function of integral of envelope function 
  // this function returns energy obeying envelope function
  inline G4double primaryCRenvelope1_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    G4double a_envelope = 
      12.0*log(cor/rigidity(lowE_primary))/log(energy(cor));
    G4double A_envelope = 
      pow(rigidity(lowE_primary),12.0)/pow(lowE_primary,a_envelope);
    return pow( (a_envelope-a_primary+1) / (A_primary * A_envelope * pow(cor, -12.0))
		* value, 1./(a_envelope-a_primary+1));
  }

  // envelope function in higher energy
  inline G4double primaryCRenvelope2(G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    return A_primary * pow(E, -a_primary);
  }

  // integral of envelope function in higher energy
  inline G4double primaryCRenvelope2_integral(G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    return A_primary/(-a_primary+1) * pow(E, -a_primary+1);
  }

  // inverse function of integral of envelope function 
  // this function returns energy obeying envelope function
  inline G4double primaryCRenvelope2_integral_inv(G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return pow((-a_primary+1)/ A_primary * value , 1./(-a_primary+1));
  }

  // The random number generator for the primary component
  G4double primaryCRenergy(HepRandomEngine* engine, G4double cor, G4double solarPotential){
    G4double rand_min_1 = 
      primaryCRenvelope1_integral(lowE_primary, cor, solarPotential);
    G4double rand_max_1 = 
      primaryCRenvelope1_integral(energy(cor), cor, solarPotential);
    G4double rand_min_2 =
      primaryCRenvelope2_integral(energy(cor), cor, solarPotential);
    G4double rand_max_2 =
      primaryCRenvelope2_integral(highE_primary, cor, solarPotential);
    /*****
    G4cerr << "rand_min_1: " << rand_min_1 << ", rand_max_1: " << rand_max_1
           << G4endl;
    G4cerr << "rand_min_2: " << rand_min_2 << ", rand_max_2: " << rand_max_2
           << G4endl;
    *****/

    G4double envelope1_area = rand_max_1 - rand_min_1;
    G4double envelope2_area = rand_max_2 - rand_min_2;

    double r, E; // E means energy
    while(1){
      if (engine->flat() <= envelope1_area/(envelope1_area + envelope2_area)){
        // envelop in lower energy
        r = engine->flat() * (rand_max_1 - rand_min_1) + rand_min_1;
        E = primaryCRenvelope1_integral_inv(r, cor, solarPotential);
        if (engine->flat() <= primaryCRspec(E, cor, solarPotential) / primaryCRenvelope1(E, cor, solarPotential))
          break;
      }
      else{
        // envelope in higher energy
        r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
        E = primaryCRenvelope2_integral_inv(r, cor, solarPotential);
        if (engine->flat() <= primaryCRspec(E, cor, solarPotential) / primaryCRenvelope2(E, cor, solarPotential))
          break;
      }
    }
    return E; //THB *GeV;
  }
  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrProtonPrimary::CrProtonPrimary()
{
  ;
}


CrProtonPrimary::~CrProtonPrimary()
{
  ;
}


std::pair<double,double> CrProtonPrimary::dir(double energy, HepRandomEngine* engine) const
  // return: cos(zenith_angle) and azimuth [rad]
  // The downward has plus sign in cos(zenith_angle),
  // and azimuth = 0 from the east, +pi/2 from the north.
{
  // Assuming isotropic from the upper (i.e. sky side) hemisphere.
  // With the cosine correction, the zenith-angle distribution should
  // be sin(theta).

  double  zenith_angle = acos(engine->flat());//THB * rad;
  double  azimuth      = engine->flat() * 2 * pi;

  return std::pair<double,double>(cos(zenith_angle), azimuth); //THB / rad);
}


double CrProtonPrimary::energySrc(HepRandomEngine* engine) const
{
  return primaryCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


double CrProtonPrimary::flux(double) const
{
  /*****
  // Integrated over the upper (sky-side) hemisphere.
  return  2 * pi * INTEGRAL_primary;  // [m**-2 s**-1]
  *****/
  // Integrated over the upper (sky-side) hemisphere.
  return  0.5 * INTEGRAL_primary;  // [m**-2 s**-1 Sr**-1]
}

double CrProtonPrimary::solidAngle() const
{
  return  2 * pi;
}


const char* CrProtonPrimary::particleName() const
{
  return "proton";
}


float CrProtonPrimary::operator()(float r)
{
  HepJamesRandom  engine;
  engine.setSeed(r * 900000000);
  // 900000000 comes from HepJamesRandom::setSeed function's comment...

  return (float)energySrc(&engine);
}


double CrProtonPrimary::calculate_rate(double old_rate)
{
  return  old_rate;
}


float CrProtonPrimary::flux(float latitude, float longitude) const
{
  return  flux(0.);
}


float CrProtonPrimary::flux(std::pair<double,double> coords) const
{
  return  flux(0.);
}


std::string CrProtonPrimary::title() const
{
  return  "CrProtonPrimary";
}


float CrProtonPrimary::fraction(float energy)
// This function doesn't work in this class... :-(
{
  return  0;
}


std::pair<float,float> CrProtonPrimary::dir(float energy) const
{
  HepJamesRandom  engine;

  std::pair<double,double>  d = dir(energy, &engine);
  
  return  std::pair<float,float>(d.first, d.second);
}

