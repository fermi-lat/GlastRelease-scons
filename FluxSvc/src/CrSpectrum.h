//$Header$

#ifndef CrSpectrum_H
#define CrSpectrum_H

#ifndef BALLOONSIM
#define BALLOONSIM
#endif


//! CrSpectrum:
//! The base class for cosmic-ray spectrum generator class.
//! All the CR spectrum generator should be derived from this class.


#include <utility>
#include "FluxSvc/Spectrum.h"
class HepRandomEngine;

class CrSpectrum : public Spectrum
{
public:
  CrSpectrum();
  ~CrSpectrum();

  /// time is the elapsed seconds from 2000-01-01 00:00:00 UT and can be negative value.
  void setTime(double time);
  /// The unit of coordinates is [deg].
  void setPosition(double latitude, double longitude, double time);
  void setPosition(double latitude, double longitude);
  /// [deg]
  double longitude() const;
  /// [deg]
  double latitude() const;
  /// [deg]
  double geomagneticLongitude() const;
  /// [deg]
  double geomagneticLatitude() const;
  /// [GV]
  double cutOffRigidity() const;
  /// [MV]
  double solarWindPotential() const;
  

#ifdef BALLOONSIM
  virtual double energySrc(HepRandomEngine* engine) const = 0;
  virtual std::pair<double,double> dir(double energy, HepRandomEngine* engine)const  = 0;
#endif // BALLOONSIM

protected:
  /// [s]
  double m_time;
  /// [deg]
  double m_geomagneticLatitude;
  /// [deg]
  double m_geomagneticLongitude;
  /// [GV]
  double m_cutOffRigidity;
  /// [MV]
  double m_solarWindPotential;
};

#endif // CrSpectrum_H
