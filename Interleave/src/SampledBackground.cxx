/** @file SampledBackground.cxx

@brief declaration and definition of the class SampledBackground

$Header$

*/
#include "SampledBackground.h"
#include <cmath>
#include "InterleaveAlg.h" // for static rate to use

#include "flux/EventSource.h"


SampledBackground::SampledBackground(const std::string& params)
{}

  // Gives back energy
double SampledBackground::energy(double time){ return 0.;}

  // Gives back paticle direction in cos(theta) and phi[rad]
   std::pair<double,double> SampledBackground::dir(double energy)
   { 
       return std::make_pair(0., 0.);
   }

   // Gives back the total flux into the standard area 
   double SampledBackground::flux(double time) const // calculate the flux [c/s/m^2/sr]
   { 
       return  InterleaveAlg::triggerRate()/(EventSource::totalArea() * solidAngle());
   }


  // Gives back solid angle from which particles come
   double SampledBackground::solidAngle() const
   { return 4*M_PI;}
  
  // Gives back the interval to the next event
  // double SampledBackground::interval(double time)
  // { return 1.0;}
