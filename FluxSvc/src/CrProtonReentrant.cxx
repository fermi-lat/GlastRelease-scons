/**
 * CrProtonReentrant:
 * The primary particle generator for reentrant secondary
 * cosmic-ray protons.
 *
 * First written by M. Ozaki(ISAS) and developed by K. Hirano(Hiroshima Univ.)
 * 2001-04-26 modified by T. Mizuno(Hiroshima Univ.)
 *   energy spectrum and angular distribution are the same as those of splash proton
 */

#include <math.h>

// Geant4
//#include "Randomize.hh"		// in source/global/HEPRandom/include/
//#include "SystemOfUnits.h"	// in source/global/management/include/
// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrProtonReentrant.h"


typedef  double G4double;


// private function definitions.
namespace {
  const G4double pi = 3.14159265358979323846264339;
  // rest energy of proton in units of GeV
  const G4double restE = 0.938;
  // lower and higher energy limit of secondary proton in units of GeV
  const G4double lowE_reentrant  = 0.1;
  const G4double highE_reentrant = 10.0;

  // The following value were computed from the model in advance.
  // These should be changed when the model is changed.
  // The integral should be done from lowE_* to highE_*.
  // FIX ME!!
  // integrated downward flux in units of [m**-2 s**-1 sr**-1]
  const G4double INTEGRAL_reentrant = 26.449;


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
   *  Generate a random distribution of secondary cosmic ray protons:
   *  j(E) = A * E**-a * exp(-pow(energy/LOW_CUTOFF,-1))
   *    A = 3.0e-03
   *    a = 2.5
   *    LOW_CUTOFF = 0.16
   *  References: 
   *  Energy spectrum above is obtained by eye-ball fitting of 
   *    AMS data (Alcaraz et al. 2000, Phys. Letter B 472, 215)
   */

  // Normalization factor of incident spectrum
  const G4double A_reentrant = 3.0e-3;
  // Differential spectral index
  const G4double a_reentrant = 2.5;
  // cutoff energy
  G4double LOW_CUTOFF = 0.16;


  // spectral model of cosmic-ray reentrant protons
  inline G4double reentrantCRspec(G4double E /* GeV */, G4double cor /* GV*/, G4double phi /* MV */){
    return A_reentrant * pow(E, -a_reentrant) * exp(-pow(E/LOW_CUTOFF,-1));
  }

  // envelope function
  inline G4double reentrantCRenvelope(G4double E /* GeV */, G4double cor /* GV*/, G4double phi /* MV */){
    return 2.04e-3*pow(LOW_CUTOFF,-0.3)*pow(E,-2.2)
      *exp(-pow(E/(LOW_CUTOFF*0.65),-1.2));
  }

  // simplified integral of envelope function
  inline G4double reentrantCRenvelope_integral
  (G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */){
    return exp(-pow(E/(LOW_CUTOFF*0.65), -1.2));
  }

  // inverse function of simplified integral of envelope function
  inline G4double reentrantCRenvelope_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return LOW_CUTOFF * 0.65 * pow(-log(value), 1./-1.2);
  }


  // returns energy obeying cosmic-ray reentrant proton's spectrum
  G4double reentrantCRenergy(HepRandomEngine* engine, G4double cor, G4double solarPotential){
    G4double rand_min =
      reentrantCRenvelope_integral(lowE_reentrant, cor, solarPotential);
    G4double rand_max =
      reentrantCRenvelope_integral(highE_reentrant, cor, solarPotential);

    G4double r, E;

    while(1){
      r = engine->flat() * (rand_max - rand_min) + rand_min;
      E = reentrantCRenvelope_integral_inv(r, cor, solarPotential);
      if (engine->flat() * reentrantCRenvelope(E, cor, solarPotential) 
	  < reentrantCRspec(E, cor, solarPotential)){
        break;
      }
    }
    return E;//THB * GeV;
  }

} // End of noname-namespace: private function definitions.


//
//
//

CrProtonReentrant::CrProtonReentrant()
{
  ;
}


CrProtonReentrant::~CrProtonReentrant()
{
  ;
}


std::pair<double,double> CrProtonReentrant::dir(double energy, HepRandomEngine* engine) const
  // return: cos(zenith_angle) and azimuth [rad]
  // The downward has plus sign in cos(zenith_angle),
  // and azimuth = 0 from the east, +pi/2 from the north.
{
  /***
   * Here we assume that
   *  zenith angle(theta): the flux (per steradian) is proportional to 1 + 0.6*sin(theta)
   *  azimuth angle(phi) : isotropic
   *
   * Reference:
   *   Tylka, A. J. 2000-05-12, GLAST team internal report "A Review of Cosmic-Ray Albedo Studies: 1949-1970" (1 + 0.6 sin(theta))
   */

  G4double zenith_angle;
  while (1){
    zenith_angle = acos(engine->flat());//THB*rad;
    if (engine->flat()*1.6<1+0.6*sin(zenith_angle)){break;}
  }

  G4double azimuth = engine->flat() * 2 * pi;

  return std::pair<double,double>(cos(zenith_angle), azimuth); //THB / rad);
}


double CrProtonReentrant::energySrc(HepRandomEngine* engine) const
{
  return reentrantCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


double CrProtonReentrant::flux(double) const
{
  /*****
  // Integrated over the upper (sky-side) hemisphere.
  // Integral of sin(theta)*(1+0.6sin(theta)) in [0,pi/2].
  return  2*pi * (1 + 0.15*pi) * INTEGRAL_reentrant;  // [m**-2 s**-1]
  *****/
  // Integrated over the upper (sky-side) hemisphere.
  // Integral of sin(theta)*(1+0.6sin(theta)) in [0,pi/2].
  return  0.5 * (1 + 0.15*pi) * INTEGRAL_reentrant;  // [m**-2 s**-1 Sr**-1]
}


double CrProtonReentrant::solidAngle() const
{
  return 2 * pi;
}


const char* CrProtonReentrant::particleName() const
{
  return "proton";
}


float CrProtonReentrant::operator()(float r)
{
  HepJamesRandom  engine;
  engine.setSeed(r * 900000000);
  // 900000000 comes from HepJamesRandom::setSeed function's comment...
  
  return (float)energySrc(&engine);
}


double CrProtonReentrant::calculate_rate(double old_rate)
{
  return  old_rate;
}


float CrProtonReentrant::flux(float latitude, float longitude) const
{
  return  flux(0.);
}


float CrProtonReentrant::flux(std::pair<double,double> coords) const
{
  return  flux(0.);
}


std::string CrProtonReentrant::title() const
{
  return  "CrProtonReentrant";
}


float CrProtonReentrant::fraction(float energy)
// This function doesn't work in this class... :-(
{
  return  0;
}


std::pair<float,float> CrProtonReentrant::dir(float energy) const
{
  HepJamesRandom  engine;

  std::pair<G4double,G4double>  d = dir(energy, &engine);
    
  return  std::pair<float,float>(d.first, d.second);
}

