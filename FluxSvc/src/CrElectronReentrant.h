//$Header$
#ifndef CrElectronReentrant_H
#define CrElectronReentrant_H

#include <utility>
#include <string>
#include "CrSpectrum.h"

// Forward declaration:
class HepRandomEngine;

/*!
 * CrElectronReentrant:
 *   The reentrant cosmic-ray electron spectrum (and incident angle) source.
 */

class CrElectronReentrant : public CrSpectrum
{
public:
  CrElectronReentrant();
  ~CrElectronReentrant();

  std::pair<double,double> dir(double energy, HepRandomEngine* engine) const;
  double energySrc(HepRandomEngine* engine) const;
  /// flux returns the value integrated over whole energy and direction
  /// and devided by 4pi Sr: then the unit is [s**-1 m**-2 Sr**-1]
  double flux() const;
  double solidAngle() const;
  const char* particleName() const;

  //
  // "flux" package stuff
  //
  /// r in [0,1)
  float operator()(float r);
  double calculate_rate(double old_rate);
  float flux(float latitude, float longitude) const;
  float flux(std::pair<double,double> coords) const;
  std::string title() const;
  /// fraction function doesn't work in this class... :-(
  float fraction(float energy);
  std::pair<float,float> dir(float energy) const;
};
#endif // CrElectronReentrant_H
