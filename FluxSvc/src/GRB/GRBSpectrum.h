/**
 * GRBSpectrum: Spectrum Class for GRB Source Simulation
 *
 * /author Nicola Omodei nicola.omodei@pi.infn.it 
 * /author Johann Cohen Tanugi johann.cohen@pi.infn.it
 *
 * \b revision:
 *   02/15/2002
 */

#ifndef GRBSpectrum_H
#define GRBSpectrum_H
// $Heading:$
//

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include "FluxSvc/mainpage.h"
#include "FluxSvc/Spectrum.h"
#include "facilities/Observer.h"
#include "src/GPS.h"
#include "CLHEP/Random/RandomEngine.h"
#include "GRBSim.h"

///  Spectrum that reads its differential spectrum from a table
class GRBSpectrum : public Spectrum
{
 public:
  /// params is the filename to read
  GRBSpectrum(const std::string& /*params*/);
  ~GRBSpectrum();
  
  /// calculate the flux, particles/m^2/s, given a time
  double flux(double time)const;
  /// return rate,given time;
  double rate(double time)const;
  
  double solidAngle() const;
  
  std::pair<float,float> dir(float energy) const;
  
  /// sample a single particle energy from the spectrum
  float operator() (float) const ;
  
  double energySrc(HepRandomEngine*, double /*time*/ );
  
  inline std::string title() const {return "GRBSpectrum";}
  inline const char * particleName() const {return "GRBgamma";}
  inline  const char * nameOf() const {return "GRBSpectrum";}
  
  
 private:
  GRBSim*              m_grbsim;
  std::vector<double>  m_spectrum;
};
#endif
