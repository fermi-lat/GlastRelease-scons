/** @file SampledBackground.h

*/
#ifndef interleave_PseudoBackground_h
#define interleave_PseudoBackground_h

#include "flux/Spectrum.h"

class SampledBackground: public Spectrum {
public:
     SampledBackground(const std::string& params);

     virtual ~SampledBackground(){};
  
  // Gives back energy
  virtual double energy(double time);

  // Gives back paticle direction in cos(theta) and phi[rad]
  virtual std::pair<double,double> dir(double energy);

  // Gives back the total flux (summation of each component's flux)
  virtual double flux(double time) const;  // calculate the flux [c/s/m^2/sr]

  // Gives back the particle kind
  virtual const char* particleName() const{ return "p"; }

  // Gives back the component name
  virtual std::string title() const{ return "SampledBackground"; }

  // Gives back solid angle from which particles come
  virtual double solidAngle() const;
  
  // Gives back the interval to the next event
  //virtual double interval(double time);


private:
};

#endif
