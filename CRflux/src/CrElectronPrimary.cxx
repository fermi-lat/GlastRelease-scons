/****************************************************************************
 * CrElectronPrimary.cc:
 ****************************************************************************
 * Read comments at the head of CosmicRayGeneratorAction.cc
 * and GLAST-LAT Technical Note (LAT-TD-250.1) by T. Mizuno et al. 
 * for overall flow of cosmic ray generation in BalloonTestV13.
 * This program is interfaced to CosmicRayGeneratorAction.cc 
 * via CrElectron, the entry-point class for the cosmic-ray electron 
 * generation.
 ****************************************************************************
 * This program provides the primary component of cosmic ray electrons.
 * Its angular distribution is assumed to be uniform for downward 
 * (theta = pi - zenith angle is 0 to pi/2) and zero for theta > pi/2.
 * The following itemizes other important features.
 * 1) One intrinsic spectrum is assumed common to all locations 
 *    on earth: a power-low with the geomagnetic cut-off at lower energy.  
 *    It is modified by the solar modulation, the geomagnetic cutoff, 
 *    and the east-west asymmetry effect (not yet implemented).
 * 2) Note that there are two other components of cosmic ray electrons. 
 *    "reentrant" and "splash".  They are handled in other programs.
 * 3) A method CrElectronPrimary::energySrc returns an energy and 
 *    CrElectronPrimary::dir returns a direction 
 *    in cos(theta) and phi (azimuth angle). 
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
 * 2001-04 Written by M. Ozaki (ISAS).
 * 2001-05 Spectrum part modified by T. Mizuno (Hiroshima Univ.)
 * 2001-05 Comments added and unused codes removed by T. Kamae (SLAC)
 * 2001-05 Program checked by T. Kamae and H. Mizushima (Hiroshima)
 * 2001-12 Modified by T. Mizuno to construct a `stand-alone' module
 * 2004-04 Model parameters are modified; functions are based on the
 *         review by Webber.
 ****************************************************************************
 */

// $Header$

#include <math.h>

// CLHEP
#include <CLHEP/config/CLHEP.h>
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>
#include <CLHEP/Random/JamesRandom.h>

#include "CrElectronPrimary.hh"


typedef double G4double;


// private function definitions.
namespace {
  // The rest energy (rest mass) of electron in [GeV]
  const G4double restE = 5.11e-4; // rest energy of electron in [GeV]
  // The lower and higher (kinetic) energy limits of primary electrons 
  // generated in this program in units of GeV.  These values are set 
  // in constructor and when the satellite position is set
  G4double lowE_primary;
  G4double highE_primary;
  // The CutE_primary is the kinetic Energy corresponding to 
  // cutoff-rigidity.   The value is set in constructor.
  G4double cutE_primary;

  // Gives back v/c as a function of kinetic Energy
  inline G4double beta(G4double E /* GeV */)
  {
#if 0	// if E ~ restE
    return sqrt(1 - pow(E/restE+1, -2));
#else	// if E >> restE
    return 1.0;
#endif
  }


  // Gives back the rigidity (p/e where p is the momentum and e means
  // electron charge magnitude) in units of [GV]
  // as a function of the kinetic energy [GeV].
  inline G4double rigidity(G4double E /* GeV */)
  {
#if 0	// if E ~ restE
    return sqrt(pow(E + restE, 2) - pow(restE, 2));
#else	// if E >> restE
    return E;
#endif
  }


  // Gives back the kinetic energy [GeV] for a rigidity [GV].
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
   * Generate a random distribution of primary cosmic ray electrons:
   * j(E) = mod_spec(E, phi) * geomag_cut(E, CutOff)
   *   mod_spec(E, phi) = org_spec(E+phi*1e-3) * 
   *    ((E+restE)**2 - restE**2) / ((E+restE+phi*1e-3)**2 - restE**2)
   *   org_spec(E) = A * rigidity(E)**-a
   *     A = 6.5e-1
   *     a = 3.3
   *   rigidity(E) = sqrt((E+restE)**2 - restE**2)
   *   beta(E) = sqrt(1 - (E/restE+1)**-2)
   *   geomag_cut(E, CutOff) = 1/(1 + (rigidity(E)/CutOff)**-6.0)
   *     CutOff = 4.46 for Theta_M = 0.735 and altitude = 35km 
   *                                           (balloon experiment)
   *     phi = 540, 1100 [MV] for Solar minimum, maximum
   *   E: [GeV]
   *   j: [c/s/m^2/sr/MeV]
   *
   * References:
   *  org_spec: Komori, Y. et al. 1999, Proc. of Dai-Kikyu (large balloon) 
   *                                    Sympo. Heisei-11yr, 33-36  (Fig. 2)
   *            R.L.Golden et al. 1996, ApJL 457, 103
   *            Webber 1983, Composition and Origin of Cosmic Rays
   *              (ed. M. M. Shapiro)
   *            Longair, M. S. 1992, High Energy Astrophysics, vol 1
   *              (2nd edition, Cambridge University Press)
   *  mod_spec: Gleeson, L. J. and Axford, W. I. 
   *           1968, ApJ, 154, 1011-1026 (Eq. 11)
   *  geomag_cut formula: 
   *    an eyeball fitting function to represent the AMS proton data.
   *  CutOff: calculated as (Rc/GV) = 14.9 * (1+h/R)^-2 * (cos(theta_M))^4,
   *          where h gives altitude from earth surface, 
   *                R means an radius of earth,
   *            and theta_M is geomagnetic lattitude. See references below.
   *          "Handbook of space astronomy and astrophysics" 2nd edition, p225 
   *            (Zombeck, 1990, Cambridge University Press)
   *          "High Energy Astrophysics" 2nd edition, p325-330
   *            (M. S. Longair, 1992, Cambridge University Press)
   */

  // Normalization factor of the original spec.
  const G4double A_primary = 6.5e-1;	
  // Differential spectral index
  const G4double a_primary = 3.3;


  // Gives back the geomagnetic cutoff factor to the intrinsic 
  // primary cosmic ray spectrum for a kinetic energy E (GeV) 
  // and a cut-off rigidity value (GV).
  inline G4double geomag_cut(G4double E /* GeV */, G4double cor /* GV */)
  {
    return 1./(1 + pow(rigidity(E)/cor, -6.0));
  }


  // The unmodulated primary electron spectrum outside the Solar system is 
  // returned.
  inline G4double org_spec(G4double E /* GeV */)
  {
    return A_primary * pow(rigidity(E), -a_primary);
  }


  // The modulated electron flux for a "phi" value is calculated 
  // according to the force-field approximation of the solar modulation.
  inline G4double mod_spec(G4double E /* GeV */, G4double phi /* MV */)
  {
    return org_spec(E + phi * 1e-3) * (pow(E+restE, 2) - pow(restE, 2))
      / (pow(E+restE+phi*1e-3, 2) - pow(restE, 2));
  }


  // The final spectrum for the primary electron.
  inline G4double primaryCRspec
  (G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */)
  {
    return mod_spec(E, phi) * geomag_cut(E, cor);
  }


  // To generate the spectrum (primaryCRspec) efficiently (ie. minimum 
  // call of random numbers), the spectrum is devided into 2 parts:
  // between lowE_primary and cutE_primary and between cutE_primary
  // and highE_primary.  The higher portion is generated by using
  // the inverse function of the integral of the power-low spectrum.
  // The lower portion takes a complicated formula and the inverse 
  // function of its integral doesn't come easily.  A simpler function
  // that envelopes the true spectrum is generated instead.  The 
  // excess is then removed by comparing with the true spectrum.
  // A linear function is used to envelope the lower part of the 
  // spectrum.

  // Envelope function in the lower energy range
  // (E<Ec, where Ec corresponds to cutoff rigidity) is given.
  // Primary spectrum of cosmic-ray electron is enveloped by a linear function
  // between lowE_primary and cutE_primary
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


  // Envelope function in higher energy
  // (E>Ec, where Ec corresponds to cutoff rigidity) is given.
  // Primary spectrum of cosmic-ray electron is enveloped by power-law function
  // between cutE_primary and highE_primary
  inline G4double primaryCRenvelope2
  (G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */)
  {
    return A_primary * pow(E, -a_primary);
  }


  // Integrated envelope function in the higher energy range.
  inline G4double primaryCRenvelope2_integral
  (G4double E /* GeV */, G4double cor /* GV */, G4double phi /* MV */)
  {
    return A_primary * pow(E, -a_primary+1) / (-a_primary+1);
  }


  // Inverse function of the integrated envelope function
  // in the higher energy range.
  inline G4double primaryCRenvelope2_integral_rev
  (G4double value, G4double cor /* GV */, G4double phi /* MV */)
  {
    return pow(value * (-a_primary+1) / A_primary, 1./(-a_primary+1));
  }


  // The rundam number generator for the primary component.
  G4double primaryCRenergy(HepRandomEngine* engine, 
			   G4double cor, G4double solarPotential)
  {
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

    while (1){
      if (engine->flat() <= 
	  envelope1_area / (envelope1_area + envelope2_area)){
        // The envelope function for the lower energy part:
        // We enveloped the spectrum by a linear function between
        // lowE_primary and cutE_primary, and assume the flux to be
        // zero below lowE_primary.
	G4double E1, E2;
        E1 = engine->flat() * (cutE_primary-lowE_primary) + lowE_primary;
        E2 = engine->flat() * (cutE_primary-lowE_primary) + lowE_primary;
        if (E1>E2){E=E1;} else {E=E2;}
        if (engine->flat() <= 
	    primaryCRspec(E, cor, solarPotential) 
	    / primaryCRenvelope1(E, cor, solarPotential))
          break;
      } else {
        // Envelope in the higher energy range.
        r = engine->flat() * (rand_max_2 - rand_min_2) + rand_min_2;
        E = primaryCRenvelope2_integral_rev(r, cor, solarPotential);
        if (engine->flat() <= 
	    primaryCRspec(E, cor, solarPotential) 
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
    {122.24, 89.316, 67.518, 52.466, 41.709, 33.798, 27.837}, // COR = 0.5GV
    {67.749, 53.615, 43.251, 35.462, 29.486, 24.817, 21.112}, // COR = 1 GV
    {27.631, 23.653, 20.429, 17.784, 15.589, 13.752, 12.201}, // COR = 2 GV
    {14.416, 12.827, 11.473, 10.31, 9.304, 8.43, 7.666}, // COR = 3 GV
    {8.675, 7.899, 7.218, 6.615, 6.081, 5.605, 5.18}, // COR = 4 GV
    {5.721, 5.29, 4.904, 4.556, 4.242, 3.957, 3.698}, // COR = 5 GV
    {4.021, 3.759, 3.521, 3.304, 3.105, 2.922, 2.754}, // COR = 6 GV
    {2.961, 2.791, 2.635, 2.491, 2.358, 2.234, 2.12}, // COR = 7 GV
    {2.26, 2.144, 2.036, 1.936, 1.843, 1.756, 1.675}, // COR = 8 GV
    {1.774, 1.692, 1.615, 1.543, 1.475, 1.412, 1.352}, // COR = 9 GV
    {1.424, 1.364, 1.308, 1.254, 1.204, 1.156, 1.111}, // COR = 10 GV
    {1.166, 1.12, 1.077, 1.037, 0.998, 0.962, 0.927}, // COR = 11 GV
    {0.969, 0.934, 0.901, 0.87, 0.84, 0.811, 0.784}, // COR = 12 GV
    {0.816, 0.789, 0.763, 0.738, 0.714, 0.692, 0.67}, // COR = 13 GV
    {0.696, 0.674, 0.653, 0.633, 0.614, 0.596, 0.578}, // COR = 14 GV
    {0.599, 0.582, 0.565, 0.548, 0.533, 0.518, 0.503} // COR = 15 GV
  };

  //============================================================

} // End of noname-namespace: private function definitions.


//
//
//

CrElectronPrimary::CrElectronPrimary()
{
  // Set lower and higher energy limits of primary electron.
  // At lowE_primary, flux of primary electron can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 100.0;
  // "cutE_primary" is the kinetic energy corresponding to 
  // cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);
}


CrElectronPrimary::~CrElectronPrimary()
{
  ;
}


// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrElectronPrimary::setPosition
(double latitude, double longitude){
  CrSpectrum::setPosition(latitude, longitude);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary proton can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 100.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}


// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrElectronPrimary::setPosition
(double latitude, double longitude, double time){
  CrSpectrum::setPosition(latitude, longitude, time);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary proton can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 100.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}


// Set satellite position and calculate energies related to COR.
// These energies will be used to generate particles.
void CrElectronPrimary::setPosition
(double latitude, double longitude, double time, double altitude){
  CrSpectrum::setPosition(latitude, longitude, time, altitude);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary proton can be 
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 100.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Set geomagnetic cutoff rigidity and calculate the energies related.
// These energies are used to generate the particle. 
void CrElectronPrimary::setCutOffRigidity(double cor){
  CrSpectrum::setCutOffRigidity(cor);

  // Set lower and higher energy limit of the primary proton (GeV).
  // At lowE_primary, flux of primary proton can be
  // assumed to be 0, due to geomagnetic cutoff
  lowE_primary = energy(m_cutOffRigidity/2.5);
  highE_primary = 100.0;
  // energy(GeV) corresponds to cutoff-rigidity(GV)
  cutE_primary = energy(m_cutOffRigidity);

}

// Gives back particle direction in (cos(theta), phi)
std::pair<double,double> CrElectronPrimary::dir(double energy, 
						HepRandomEngine* engine) const
  // return: cos(theta) and phi [rad]
  // The downward has plus sign in cos(theta),
  // and phi = 0 for the particle comming along x-axis (from x>0 to x=0)
  // and phi=pi/2 for that comming along y-axis (from y>0 to y=0).
{
  // Assuming isotropic from the upper (i.e. sky side) hemisphere.
  // After integration over the phi, the theta distribution should
  // be sin(theta) for a constant theta width.

  /// Cos(theta) ranges from 1 to -0.4
  double theta = acos(1.4*engine->flat()-0.4);
  double phi   = engine->flat() * 2 * M_PI;

  return std::pair<double,double>(cos(theta), phi);
}



// Gives back particle energy
double CrElectronPrimary::energySrc(HepRandomEngine* engine) const
{
  return primaryCRenergy(engine, m_cutOffRigidity, m_solarWindPotential);
}


// flux() returns the energy integrated flux averaged over
// the region from which particle is coming from 
// and the unit is [c/s/m^2/sr].
// flux()*solidAngle() is used as relative normalization among
// "primary", "reentrant" and "splash".
double CrElectronPrimary::flux() const
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
    phi = phi/100.0 - 5; // 500 MV corresponds to 0, 600 MV corresponds to 1, etc..
    G4double tmp1 = 
      integral_array[int(cor)][int(phi)] +
      (cor - int(cor)) * (integral_array[int(cor)+1][int(phi)]-integral_array[int(cor)][int(phi)]);
    G4double tmp2 = 
      integral_array[int(cor)][int(phi)+1] +
      (cor - int(cor)) * (integral_array[int(cor)+1][int(phi)+1]-integral_array[int(cor)][int(phi)+1]);
    energy_integral = tmp1 + (tmp2-tmp1)*(phi-int(phi));
  }

  if (cor < 1.0){
    phi = phi/100.0 - 5; // 500 MV corresponds to 0, 600 MV corresponds to 1, etc..
    G4double tmp1 = 
      integral_array[0][int(phi)] +
      (2*cor - 1) *(integral_array[1][int(phi)]-integral_array[0][int(phi)]);
    G4double tmp2 = 
      integral_array[0][int(phi)+1] +
      (2*cor - 1) *(integral_array[1][int(phi)+1]-integral_array[0][int(phi)+1]);
    energy_integral = tmp1 + (tmp2-tmp1)*(phi-int(phi));
  }

  // We assume that the flux is uniform above the earth horizon.
  // Then the average flux is equal to the vertically downward one.
  return  energy_integral;  // [c/s/m^2/sr]
}


// Gives back solid angle from which particle comes
double CrElectronPrimary::solidAngle() const
{
  // * 1.4 since Cos(theta) ranges from 1 to -0.4
  return  2 * M_PI * 1.4;
}


// Gives back particle name
const char* CrElectronPrimary::particleName() const
{
  return "e-";
}


// Gives back the name of the component
std::string CrElectronPrimary::title() const
{
  return  "CrElectronPrimary";
}

