/**
 * CrElectronReentrant:
 * The primary particle generator for reentrant secondary
 * cosmic-ray electrons.
 *
 * Ver 1.0 on 2001-04-18 by Masanobu Ozaki <ozaki@astro.isas.ac.jp>
 *    separate the reentrant component from the all-in-one code.
 * 2001-04-27 modified by T. Mizuno
 *   model function is based on AMS data)
 *   spectrum and angular distribution are the same as those of splash electron
 *
 * $Header$
 */

#include <math.h>

// Geant4
//#include "Randomize.hh"		// in source/global/HEPRandom/include/
//#include "SystemOfUnits.h"	// in source/global/management/include/
// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrElectronReentrant.h"


typedef  double G4double;


// private function definitions.
namespace {
  const G4double pi    = 3.14159265358979323846264339;
  const G4double restE = 5.11e-4; // rest energy of electron in [GeV]
  // lower and higher energy limit of re-entrant electron in units of GeV
  const G4double lowE_reent  = 0.1;
  const G4double highE_reent = 10.0;

  // The following value were computed from the model in advance.
  // These should be changed when the model is changed.
  // The integral should be done from lowE_* to highE_*.
  // FIX ME!!
  const G4double INTEGRAL_reent = 17.812; // [m**-2 s**-1 Sr**-1]


  // gives v/c as a function of kinetic Energy
  inline G4double beta(G4double E /* GeV */)
  {
#if 0	// if E ~ restE
    return sqrt(1 - pow(E/restE+1, -2));
#else	// if E >> restE
    return 1.0;
#endif
  }


  // gives rigidity (GV) as a function of kinetic Energy
  inline G4double rigidity(G4double E /* GeV */)
  {
#if 0	// if E ~ restE
    return sqrt(pow(E + restE, 2) - pow(restE, 2));
#else	// if E >> restE
    return E;
#endif
  }


  // gives kinetic energy (GeV) as a function of rigidity
  inline G4double energy(G4double rigidity /* GV */)
  {
#if 0	// if energy ~ restE
    return sqrt(pow(rigidity, 2) + pow(restE, 2)) - restE;
#else	// if energy >> restE
    return rigidity;
#endif
  }



  //=====================================================================
  /**
   * Generate a random distribution of re-entrant cosmic ray electrons
   * (same distribution as splash cosmic ray electrons):
   * j(E) = A * E**-a
   *     A = 1.33e-04
   *     a = 3.53
   *   E: [GeV]
   *   j: [MeV**-1 m**-2 s**-1 Sr**-1]
   *
   * References:
   *   Alcaraz, J. et al. 2000, Phys.Let.B, 484, 10-22
   */

  const G4double A_reent = 1.33e-4;
  const G4double a_reent = 3.53;

  const G4double Cut_E = 0.06; // cutoff energy (GeV) used in envelope func.
  const G4double A_Envelope = 1.5e-04;


  // spectral model of reentrant cosmic-ray electrons
  inline G4double reentrantCRspec(G4double E /* GeV */)
  {
    return A_reent * pow(E, -a_reent);
  }


  // envelope function
  inline G4double reentrantCRenvelope(G4double E /* GeV */)
  {
    return A_Envelope * pow(E, -a_reent) * exp(-pow(E/Cut_E, -a_reent+1));
  }


  // simplified integral of envelope function
  inline G4double reentrantCRenvelope_integral(G4double E /* GeV */)
  {
    return exp(-pow(E/Cut_E, -a_reent+1));
  }


  // inverse function of simplified integral of envelope function
  inline G4double reentrantCRenvelope_integral_inv(G4double value)
  {
    return Cut_E * pow(-log(value), -1./(a_reent-1));
  }


  // returns energy obeying re-entrant cosmic-ray electron spectrum
  G4double reentrantCRenergy(HepRandomEngine* engine)
  {
    G4double rand_min = reentrantCRenvelope_integral(lowE_reent);
    G4double rand_max = reentrantCRenvelope_integral(highE_reent);

    double r, E;

    while (1){
      r = engine->flat() * (rand_max - rand_min) + rand_min;
      E = reentrantCRenvelope_integral_inv(r);
      // G4cerr << "E: " << E << ", reentrantCRspec(E): " << reentrantCRspec(E) << ", "
      // G4cerr << "reentrantCRenvelope1(E)" << reentrantCRenvelope(E) << G4endl;
      if (engine->flat() <= reentrantCRspec(E) / reentrantCRenvelope(E)){
	break;
      }
    }
    return E ;//THB* GeV;
  }
} // End of noname-namespace: private function definitions.


//
//
//

CrElectronReentrant::CrElectronReentrant()
{
  ;
}


CrElectronReentrant::~CrElectronReentrant()
{
  ;
}


std::pair<double,double> CrElectronReentrant::dir(double energy, HepRandomEngine* engine) const
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

  double zenith_angle;
  while (1){
    zenith_angle = acos(engine->flat());//THB *rad;
    if (engine->flat()*1.6<1+0.6*sin(zenith_angle)){break;}
  }

  double azimuth = engine->flat() * 2 * pi;

  return std::pair<double,double>(cos(zenith_angle), azimuth ); //THB/ rad);
}


double CrElectronReentrant::energySrc(HepRandomEngine* engine) const
{
  return reentrantCRenergy(engine);
}


double CrElectronReentrant::flux(double) const
{
  /*****
  // Integrated over the upper (sky-side) hemisphere.
  // Integral of sin(theta)*(1+0.6sin(theta)) in [0,pi/2].
  return  2*pi * (1 + 0.15*pi) * INTEGRAL_reent;  // [m**-2 s**-1]
  *****/
  // Integrated over the upper (sky-side) hemisphere.
  // Integral of sin(theta)*(1+0.6sin(theta)) in [0,pi/2].
  return  0.5 * (1 + 0.15*pi) * INTEGRAL_reent;  // [m**-2 s**-1 Sr**-1]
}


double CrElectronReentrant::solidAngle() const
{
  return 2 * pi;
}


const char* CrElectronReentrant::particleName() const
{
  return "e-";
}


float CrElectronReentrant::operator()(float r)
{
  HepJamesRandom  engine;
  engine.setSeed(r * 900000000);
  // 900000000 comes from HepJamesRandom::setSeed function's comment...

  return (float)energySrc(&engine);
}


double CrElectronReentrant::calculate_rate(double old_rate)
{
  return  old_rate;
}


float CrElectronReentrant::flux(float latitude, float longitude) const
{
  return  flux(0.);
}


float CrElectronReentrant::flux(std::pair<double,double> coords) const
{
  return  flux(0.);
}


std::string CrElectronReentrant::title() const
{
  return  "CrElectronReentrant";
}


float CrElectronReentrant::fraction(float energy)
  // This function doesn't work in this class... :-(
{
  return  0;
}


std::pair<float,float> CrElectronReentrant::dir(float energy) const
{
  HepJamesRandom  engine;

  std::pair<double,double>  d = dir(energy, &engine);

  return  std::pair<float,float>(d.first, d.second);
}
