//$Header$

#ifndef CrProton_H
#define CrProton_H


//!  The class that calls each cosmic-ray proton components based on
//!  their flux.

#include <vector>
#include <utility>
#include <string>
#include "FluxSvc/Spectrum.h"

class CrSpectrum;

class CrProton : public Spectrum
{
public:
    // params[0] is bit flag for determining which to include (default 7)
    // 1: primary
    // 2: reentrant
    // 4: splash

  CrProton(const std::string& params);
  ~CrProton();

  CrSpectrum* selectComponent(HepRandomEngine* engine);
  double energySrc(HepRandomEngine* engine);
  std::pair<double,double> dir(double energy, HepRandomEngine* engine);
  CrSpectrum* component() const;

	//! calculate the flux, particles/m^2/sr.
  	virtual double    flux ( ) const;

	virtual const char * particleName()const{ return "p";}
  virtual std::string title()const{return "CrProton";}
  virtual double solidAngle( )const;

private:
  std::vector<CrSpectrum*>  m_subComponents;
  CrSpectrum*               m_component;
};
#endif // CrProton_H

