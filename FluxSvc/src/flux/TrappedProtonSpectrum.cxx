// $Id$

#include "FluxSvc/TrappedProtonSpectrum.h"

#include "FluxSvc/GPS.h"


#include <cmath>

#include "FluxSvc/SpectrumFactory.h"

static SpectrumFactory<TrappedProtonSpectrum> factory;
const ISpectrumFactory& TrappedProtonSpectrumFactory = factory;

//-------------------------- constructors

TrappedProtonSpectrum::TrappedProtonSpectrum(float lat, float lon) {
    init(lat,lon);
}

TrappedProtonSpectrum::TrappedProtonSpectrum(const std::vector<float>& params)
{
	init(params[0], params[1]);
}

//-------------------------- flux() (current position)

double TrappedProtonSpectrum::flux() const {
    // calculate flux for the current cutoff
    return flux(m_lat, m_lon);
}

//-------------------------- flux() (specified position)

float TrappedProtonSpectrum::flux(float lat, float lon) const {
    // Flux as a function of latitude and longitude in a 600 km orbit.
    // Linear interpolate in a table with a 5 degree sampling grid.
    int ilat = static_cast<int>(lat/5.+6.);
    float a = fmod(lat+30., 5.)/5.;
    int ilon = static_cast<int>(lon/5.);
    float b = fmod(lon, 5.)/5.;
    float rawflux = m_fluxTbl[ilon][ilat] * (1.-a) * (1.-b) +
	m_fluxTbl[ilon+1][ilat] * (1.-a) * b +
	m_fluxTbl[ilon][ilat+1] * a * (1.-b) +
	m_fluxTbl[ilon+1][ilat+1] * a * b;
    return 10000. * rawflux / (4.*M_PI); // convert cm^2 to m^2 ster
}

//-------------------------- flux() (position given in a pair)

float TrappedProtonSpectrum::flux(std::pair<double,double> coords) const {
    // Flux as a function of latitude and longitude
    return flux(coords.first, coords.second);
}

//-------------------------- operator()  sample an energy value

float TrappedProtonSpectrum::operator() (float x)const {
    // return a random value of energy sampled from the spectrum
    float y = 1. - x;
    if (y > 0.1107) return -0.09088 * log(y);
    else            return -0.1011  * log(y/0.8004);
}

double TrappedProtonSpectrum::calculate_rate(double old_rate) {
  setPosition( GPS::instance()->lat(), GPS::instance()->lon() );
  return flux();
}

void TrappedProtonSpectrum::init(float lat, float lon) {
    m_lat = lat;
    m_lon = lon;
	
#include "TrappedProtonSpectrum.inc"  // numerical data: energies, fluxes

    // table of total flux as a function of latitude and longitude
    for (int ii=0; ii < 73; ii++) {
	for (int jj=0; jj < 13; jj++) {
	    m_fluxTbl[ii][jj] = gfluxes[jj+13*ii];
	}
    }
}

const char * TrappedProtonSpectrum::particleName() const {
  return "p";
}

std::string TrappedProtonSpectrum::title() const {
  return "TrappedProtonSpectrum";
}
