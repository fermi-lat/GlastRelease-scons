/****************************************************************************
 * CrAlphaPrimary.cc:
 ****************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrAlpha, the entry-point class for the cosmic-ray alpha generation.
 ****************************************************************************
 * This program provides the primary component of cosmic ray alphas.
 * Its angular distribution is assumed to be uniform for downward 
 * (theta = pi - zenith angle is 0 to pi/2) and zero for theta > pi/2.
 * The following itemizes other important features.
 * 1) One intrinsic spectrum is assumed common to all locations 
 *    on earth: a power-low with the geomagnetic cut-off at lower energy.  
 *    It is modified by the solar modulation, the geomagnetic cutoff, 
 *    and the east-west asymmetry effect (not yet implemented).
 * 2) A method CrAlphaPrimary::energySrc returns an energy and 
 *    CrAlphaPrimary::dir returns a direction in 
 *    cos(theta) and phi (azimuth angle). 
 ****************************************************************************
 * Definitions:
 * 1) The z-axis points upward (from Calorimeter to Tracker).  
 * 2) Partile of theta=0 points for downward (i.e., comes from zenith)
 *    and that of theta=pi points for upward (comes from nadir).
 * 3) Partile of phi=0 comes along x-axis (from x>0 to x=0) and
 *    that of phi=pi/2 comes along y-axis (from y>0 to y=0).
 * 4) Energy means kinetic energy unless defined otherwise. 
 * 5) Magnetic latitude theta_M is in radian. 
 * 6) Particle direction is defined by cos(theta) and phi (in radian).
 ****************************************************************************
 * 2002-04 written by Y. Fukazawa (Hiroshima)
 * 2003-02 Parameters are updated by T. Mizuno
 ****************************************************************************
 */

//$Header$

#include <math.h>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrAlphaPrimary.hh"


typedef  double G4double;


// private function definitions.
namespace {
  const G4double pi = 3.14159265358979323846264339;
  // rest energy of alpha particle in units of GeV
  const G4double restE = 3.72;

  // The lower and higher (kinetic) energy limits of primary alphas
  // generated in this program in units of GeV.  These values are set 
  // in constructor and when the satellite position is set
  G4double lowE_primary;
  G4double highE_primary;
  // The CutE_primary is the kinetic Energy corresponding to 
  // cutoff-rigidity. The value is set also in constructor and 
  // when the satellite position is set.
  G4double cutE_primary;

  const G4double z_alpha = 2.0;


  // gives back v/c as a function of kinetic Energy
  inline G4double beta(G4double E /* GeV */){
    return sqrt(1 - pow(E/restE+1, -2));
  }

  // gives back the rigidity (p/Ze where p is the momentum, e means 
  // electron charge magnitude, and Z is the atomic number) in units of [GV],
  // as a function of kinetic Energy [GeV].
  inline G4double rigidity(G4double E /* GeV */){
    return sqrt(pow(E + restE, 2) - pow(restE, 2))/z_alpha;
  }

  // gives back the kinetic energy [GeV] as a function of rigidity [GV]
  inline G4double energy(G4double rigidity /* GV */){
    return sqrt(pow(rigidity*z_alpha, 2) + pow(restE, 2)) - restE;
  }

  //============================================================
  /**
   *  Generate a random distribution of primary cosmic ray alphas
   *  j(E) = mod_spec(E, phi) * geomag_cut(E, CutOff)
   *    mod_spec(E, phi) = org_spec(E+2*phi*1e-3) * 
   *     ((E+restE)**2 - restE**2)/((E+restE+2*phi*1e-3)**2-restE**2)
   *    org_spec(E) = A * rigidity(E)**-a
   *      A = 1.50 and a = 2.77
   *    rigidity(E) = sqrt((E+restE)**2 - restE**2) / 2.0
   *    geomag_cut(E, CutOff) = 1/(1 + (rigidity(E)/CutOff)**-12.0)
   *      CutOff = 4.46 for Theta_M = 0.735 and altitude = 35km 
   *                                         (balloon experiment)
   *      phi = 540 and 1100 [MV] for Solar minimum and maximum, respectively
   *    E: [GeV]
   *    j: [c/s/m^2/sr/MeV]
   *
   *  References:
   *  org_spec: Set et al. 1991, ApJ 378, 763 (LEAP)
   *            Sanuki et al. 2000, ApJ 545, 1135 (BESS)
   *  geomag_cut formula: 
   *    an eyeball fitting function to represent the AMS proton data
   *  CutOff: calculated as (Rc/GV) = 14.9 * (1+h/R)^-2 * (cos(theta_M))^4,
   *          where h is the altitude from earth surface, 
   *                R is the mean radius of earth,
   *            and theta_M is geomagnetic lattitude. 
   *          References:
   *          "Handbook of space astronomy and astrophysics" 2nd edition, p225 
   *            (Zombeck, 1990, Cambridge University Press)
   *          "High Energy Astrophysics" 2nd edition, p325-330
   *            (M. S. Longair, 1992, Cambridge University Press)
   *  Solar modulation model (mod_spec) is from:
   *     Gleeson, L. J. and Axford, W. I. 1968, ApJ, 154, 1011-1026 (Eq. 11)
   */

  const G4double A_primary = 1.50; // normalization of incident spectrum
  const G4double a_primary = 2.77; // differential spectral index

  // Gives back the geomagnetic cutoff factor to the intrinsic 
  // primary cosmic ray spectrum for a kinetic energy E(GeV) 
  // and a geomagnetic lattitude theta_M(rad)
  inline G4double geomag_cut(G4double E, G4double cor /* GV */){
    return 1./(1 + pow(rigidity(E)/cor, -12.0));
  }

  // The unmodulated primary alpha spectrum outside the Solar system is 
  // returned.
  inline G4double org_spec(G4double E /* GeV */)
  {
    return A_primary * pow(rigidity(E), -a_primary);
  }

  // The modulated alpha flux for a "phi" value is returned. 
  // Force-field approximation of the Solar modulation is used.
  // The value of phi(potential) is doubled due to the charge of 2.
  inline G4double mod_spec(G4double E /* GeV */, G4double phi /* MV */){
    return org_spec(E + z_alpha*phi*1e-3) * (pow(E+restE, 2) - pow(restE, 2))
      / (pow(E+restE+z_alpha*phi*1e-3,2) - pow(restE,2));
  }

  // The final spectrum for the primary alpha.
  inline G4double primaryCRspec
  (G4double E /* GeV */, G4double cor /* GV*/, G4double phi /* MV */){
    return mod_spec(E, phi) * geomag_cut(E, cor);
  }

  // To speed up generation of the spectrum below the geomagnetic cut 
  // off, random numbers are generated to "an envelope function" for
  // whose integral the inverse function can be found.  The final one 
  // is obtained by throwing away some random picks of energy.

  // Envelope function in the lower energy range 
  // (E<Ec, where Ec corresponds to cutoff rigidity).
  // Primary spectrum of cosmic-ray alpha is enveloped by a linear function
  // between lowE_primary and cutE_primary.
  inline G4double primaryCRenvelope1
  (G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */){
    G4double coeff = 
      ( primaryCRspec(cutE_primary,cor,phi)
        - primaryCRspec(lowE_primary,cor,phi) ) 
      / (cutE_primary-lowE_primary);
    return coeff * (E-lowE_primary) + primaryCRspec(lowE_primary,cor,phi);
  }

  // Integral of the envelope function in the lower energy range
  inline G4double primaryCRenvelope1_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    G4double coeff =
      ( primaryCRspec(cutE_primary,cor,phi)
        - primaryCRspec(lowE_primary,cor,phi) ) 
      / (cutE_primary-lowE_primary);
    return 0.5 * coeff * pow(E-lowE_primary,2) + 
      primaryCRspec(lowE_primary,cor,phi) * (E-lowE_primary);
  }

  // The envelope function in the higher energy range
  // (E>Ec, where Ec corresponds to cutoff rigidity).
  // Primary spectrum of cosmic-ray alpha is enveloped by power-law function
  // between cutE_primary and highE_primary
  inline G4double primaryCRenvelope2
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    return A_primary * pow(E/z_alpha, -a_primary);
  }

  // The integral of the envelope function in the higher energy range
  inline G4double primaryCRenvelope2_integral
  (G4double E /* GeV */, G4double cor /* MV */, G4double phi /* MV */){
    return A_primary*z_alpha/(-a_primary+1) * pow(E/z_alpha, -a_primary+1);
  }

  // The inverse function of the integral of the envelope function
  // in the higher energy range.
  // This function returns energy obeying envelope function.
  inline G4double primaryCRenvelope2_integral_inv
  (G4double value, G4double cor /* MV */, G4double phi /* MV */){
    return z_alpha*pow((-a_primary+1)/ (A_primary*z_alpha) * value , 
		       1./(-a_primary+1));
  }

  // The random number generator for the primary component
  G4double primaryCRenergy(HepRandomEngine* engine, 
                           G4double cor, G4double solarPotential){
    G4double rand_min_1 = 
      primaryCRenvelope1_integral(lowE_primary, cor, solarPotential);
    G4double rand_max_1 = 
      primaryCRenvelope1_integral(cutE_primary, cor, solarPotential);
    G4double rand_min_2 =
      primaryCRenvelope2_integral(cutE_primary, cor, solarPotential);
    G4double rand_max_2 =
      primaryCRenvelope2_integral(highE_primary, cor, solarPotential);

    G4double envelope1_area = rand_max_1 - rand_min_1;
    G4double envelope2_area = rand_max_2 - rand_min_2;

    double r, E; // E means energy in GeV
    while(1){
      if (engine->flat() <= envelope1_area/(envelope1_area + envelope2_area)){
        // Use the envelop function in the lower energy range
        // (E<Ec where Ec corresponds to the cutoff rigidity).
        // We enveloped alpha spectrum by linear function between
        // lowE_primary and cutE_primary, and assume that the flux
        // at lowE_primary is 0.
        G4double E1, E2;
        E1 = engine->flat() * (cutE_primary-lowE_primary) + lowE_primary;
        E2 = engine->flat() * (cutE_primary-lowE_primary) + lowE_primary;
        if (E1>E2){E=E1;} else {E=E2;}
        if (engine->flat() <= 
            primaryCRspec(E, cor, solarPotential) 
            / primaryCRenvelope1(E, cor, solarPotential))
          break;
      }
      else{
        // Use the envelop function in the higher energy range
        // (E>Ec where Ec corresponds to the cutoff rigidity).
        r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
        E = primaryCRenvelope2_integral_inv(r, cor, solarPotential);
        if (engine->flat() <= primaryCRspec(E, cor, solarPotential) 
            / primaryCRenvelope2(E, cor, solarPotential))
          break;
      }
    }
    return E;
  }

  // This array stores vertically downward flux in unit of [c/s/m^s/sr]
  // as a function of COR and phi (integral_array[COR][phi]).
  // The flux is integrated between lowE_primary and highE_primary.
  // COR = 0.5, 1, 2, ..., 15 [GV]
  // phi = 500, 600, ..., 1100 [MV]
  G4double integral_array[16][7] = {
    {321.7, 269.0, 229.5, 199.1, 174.4, 154.5, 138.1}, // COR = 0.5GV
    {294.1, 250.1, 216.0, 188.9, 166.9, 148.7, 133.5}, // COR = 1 GV
    {203.9, 181.2, 162.3, 146.2, 132.4, 120.6, 110.2}, // COR = 2 GV
    {136.9, 125.3, 115.2, 106.2, 98.2, 91.1, 84.7}, // COR = 3 GV
    {96.0, 89.5, 83.7, 78.3, 73.5, 69.1, 65.0}, // COR = 4 GV
    {70.6, 66.6, 63.0, 59.6, 56.5, 53.6, 50.9}, // COR = 5 GV
    {54.0, 51.4, 48.9, 46.7, 44.5, 42.5, 40.7}, // COR = 6 GV
    {42.5, 40.7, 39.0, 37.4, 35.9, 34.5, 33.2}, // COR = 7 GV
    {34.4, 33.1, 31.8, 30.7, 29.6, 28.5, 27.5}, // COR = 8 GV
    {28.3, 27.4, 26.4, 25.5, 24.7, 23.9, 23.1}, // COR = 9 GV
    {23.7, 23.0, 22.3, 21.6, 20.9, 20.3, 19.7}, // COR = 10 GV
    {20.1, 19.5, 19.0, 18.5, 17.9, 17.4, 17.0}, // COR = 11 GV
    {17.3, 16.8, 16.4, 15.9, 15.5, 15.1, 14.7}, // COR = 12 GV
    {15.0, 14.6, 14.2, 13.9, 13.5, 13.2, 12.9}, // COR = 13 GV
    {13.1, 12.8, 12.5, 12.2, 11.9, 11.6, 11.4}, // COR = 14 GV
    {11.5, 11.2, 11.0, 10.7, 10.5, 10.3, 10.1} // COR = 15 GV
  };
  


  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrAlphaPrimary::CrAlphaPrimary()
{
  // Set lower and higher energy limit of the primary alpha (GeV).
  // At lowE_primary, flux of primary alpha can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 400.0; // corresponds to 100 GeV/n
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);
}


CrAlphaPrimary::~CrAlphaPrimary()
{
  ;
}


// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrAlphaPrimary::setPosition(double latitude, double longitude){
  CrSpectrum::setPosition(latitude, longitude);

  // Set lower and higher energy limit of the primary alpha (GeV).
  // At lowE_primary, flux of primary alpha can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 400.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrAlphaPrimary::setPosition
(double latitude, double longitude, double time){
  CrSpectrum::setPosition(latitude, longitude, time);

  // Set lower and higher energy limit of the primary alpha (GeV).
  // At lowE_primary, flux of primary alpha can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 400.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);
  
}

// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrAlphaPrimary::
setPosition(double latitude, double longitude, double time, double altitude){
  CrSpectrum::setPosition(latitude, longitude, time, altitude);

  // Set lower and higher energy limit of the primary alpha (GeV).
  // At lowE_primary, flux of primary alpha can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 400.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set geomagnetic cutoff rigidity and calculate the energies related.
// These energies are used to generate the particle. 
void CrAlphaPrimary::setCutOffRigidity(double cor){
  CrSpectrum::setCutOffRigidity(cor);

  // Set lower and higher energy limit of the primary alpha (GeV).
  // At lowE_primary, flux of primary alpha can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 400.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Gives back particle direction in (cos(theta), phi)
std::pair<double,double> CrAlphaPrimary::dir(double energy, 
                                              HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward direction has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0).
{
  // We assume isotropic distribution from the upper hemisphere.
  // After integration over the azimuth angle (phi), 
  // the theta distribution should be sin(theta) for a constant theta width.

  double theta = acos(engine->flat());
  double phi   = engine->flat() * 2 * pi;

  return std::pair<double,double>(cos(theta), phi);
}


// Gives back particle energy
double CrAlphaPrimary::energySrc(HepRandomEngine* engine) const
{
  return primaryCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


// flux() returns the value integrated over whole energy and direction
// and devided by 4pi sr: then the unit is [c/s/m^2/sr].
// This value is used as relative normalization among
// "primary", "reentrant" and "splash".
double CrAlphaPrimary::flux() const
{
  // Straight downward (theta=0) flux integrated over energy,
  // given by integral_array[16][7]
  G4double energy_integral;

  // cutoff rigidity and solar potential to calculate the energy integral.
  G4double cor, phi; 

  // cor is restricted in 0.5 < cor < 14.9[GV] and 
  // phi is restricted in 500 < phi < 1100[MV]
  cor = m_cutOffRigidity;
  if (m_cutOffRigidity < 0.5){ cor = 0.5; }
  if (m_cutOffRigidity > 14.9){ cor = 14.9; }
  phi = m_solarWindPotential;
  if (m_solarWindPotential < 500.0){ phi = 500.0; }
  if (m_solarWindPotential > 1100.0){ phi = 1100.0; }

  if (cor >=1.0){
    phi = phi/100.0 - 5; // 500 MV corresponds to 0, 600 MV corresponds to 1, et..
    G4double tmp1 = 
      integral_array[int(cor)][int(phi)] +
      (cor - int(cor)) * (integral_array[int(cor)+1][int(phi)]-integral_array[int(cor)][int(phi)]);
    G4double tmp2 = 
      integral_array[int(cor)][int(phi)+1] +
      (cor - int(cor)) * (integral_array[int(cor)+1][int(phi)+1]-integral_array[int(cor)][int(phi)+1]);
    energy_integral = tmp1 + (tmp2-tmp1)*(phi-int(phi));
  }

  if (cor < 1.0){
    // 500 MV corresponds to 0, 600 MV corresponds to 1, etc..
    phi = phi/100.0 - 5; 
    G4double tmp1 = 
      integral_array[0][int(phi)] +
      (2*cor - 1) *(integral_array[1][int(phi)]-integral_array[0][int(phi)]);
    G4double tmp2 = 
      integral_array[0][int(phi)+1] +
      (2*cor - 1) *(integral_array[1][int(phi)+1]-integral_array[0][int(phi)+1]);
    energy_integral = tmp1 + (tmp2-tmp1)*(phi-int(phi));
  }

  // Integrated over the upper (sky-side) hemisphere and divided by 4pi.
  return  0.5 * energy_integral;  // [c/s/m^2/sr]


}

// Gives back solid angle from which particle comes
double CrAlphaPrimary::solidAngle() const
{
  return  2 * pi;
}


// Gives back particle name
const char* CrAlphaPrimary::particleName() const
{
  return "alpha";
}


float CrAlphaPrimary::operator()(float r)
{
  HepJamesRandom  engine;
  engine.setSeed(r * 900000000);
  // 900000000 comes from HepJamesRandom::setSeed function's comment...

  return (float)energySrc(&engine);
}


double CrAlphaPrimary::calculate_rate(double old_rate)
{
  return  old_rate;
}


float CrAlphaPrimary::flux(float latitude, float longitude) const
{
  return  flux();
}


float CrAlphaPrimary::flux(std::pair<double,double> coords) const
{
  return  flux();
}


std::string CrAlphaPrimary::title() const
{
  return  "CrAlphaPrimary";
}


float CrAlphaPrimary::fraction(float energy)
  // This function doesn't work in this class... :-(
{
  return  0;
}


std::pair<float,float> CrAlphaPrimary::dir(float energy) const
{
  HepJamesRandom  engine;

  std::pair<double,double>  d = dir(energy, &engine);
  
  return  std::pair<float,float>(d.first, d.second);
}

