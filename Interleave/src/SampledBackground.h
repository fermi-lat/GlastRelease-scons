/** @file SampledBackground.h
    @brief declaration of SampledBackground class
$Header$

*/
#ifndef interleave_PseudoBackground_h
#define interleave_PseudoBackground_h

#include "flux/Spectrum.h"

/** @class SampledBackground
    @brief  a Spectrum class that is used to implement the background sampling


*/
class SampledBackground: public Spectrum {
public:
    SampledBackground(const std::string& params);

    virtual ~SampledBackground(){};

    /// Gives back energy
    virtual double energy(double time);

    /// Gives back paticle direction in cos(theta) and phi[rad]
    virtual std::pair<double,double> dir(double energy);

    /// Gives back the total flux (summation of each component's flux)
    virtual double flux(double time) const;  // calculate the flux [c/s/m^2/sr]

    /// Gives back the particle kind
    virtual const char* particleName() const{ return "p"; }

    /// Gives back the component name
    virtual std::string title() const{ return "SampledBackground"; }

    /// Gives back solid angle from which particles come
    virtual double solidAngle() const;


private:
};

#endif
