// $id$
// File: TrappedProtonSpectrum.h
//
//
//!  Purpose: Simulate the flux and spectrum of protons trapped in
//!  the earth's magnetic field.  Implemented only for a 600 km circular
//!  orbit, so only the South Atlantic Anomaly (SAA) is seen.
//
/*!  Interface:
    In the constructor, latitude and longitude can be specified, either
       as two floats or by a pair.  Coordinates are in degrees
    After a TrappedProtonSpectrum is constructed, the position can
       be modified by using the setPosition() member, passing it latitude
       and longitude either as two floats or as a pair.  The pair form
       is especially useful in combination with an Orbit object.
    Flux in protons / (m^2 sec ster) is returned by the flux() member.
       Alternately, flux at a given position can be obtained by 
       specifying the coordinates.  The values are obtained by interpolating
       on a 5-degree grid, so finer detail should not be expected.
    Random samples from the energy spectrum (in GeV) can be obtained by using
       the () operator, supplying a random argument between 0 and 1.
       This behavior is found in all classes which inherit from Spectrum.

  Expected usage:
    To obtain flux as a function of position, use flux(lat, lon) or
    flux(pair).  To obtain spectrum samples, use ().  The setPosition
    members and flux() (with no arguments) are probably of little utility.
  Patrick Nolan, Stanford University, 1998
*/
#ifndef TRAPPED_PROTON_SPECTRUM_H
#define TRAPPED_PROTON_SPECTRUM_H

#include "FluxSvc/Spectrum.h"

#include "FluxSvc/Orbit.h"

#include <utility>


class TrappedProtonSpectrum  : public Spectrum
{
 public:
    TrappedProtonSpectrum(float lat=0., float lon=0.);
    TrappedProtonSpectrum( const std::vector<float>& params);
 private:
    void init(float lat, float lon);//initializes variables

 public:
    virtual double calculate_rate (double old_rate);
    virtual double flux() const;
    virtual float flux(float lat, float lon) const;
    virtual float flux(std::pair<double, double> coords) const;
    virtual float operator() (float)const;
    virtual  const char * particleName() const;
    virtual  std::string title() const;
	void	setPosition ( float lat, float lon ){}

    //   use default destructor, copy constructor, and assignment op.

 private:
    /// table of total flux vs. latitude and longitude
    float m_fluxTbl[73][13];
};

#endif // TRAPPED_PROTON_SPECTRUM_H

