// $Id$

#include "HeSpectrum.h"

#include "HeSpectrum.inc"  // numerical data: energies, fluxes, gfluxes

#include <cmath>
#include <algorithm>
#include <functional>
#include "CLHEP/Random/Random.h"
#include "GPS.h"

#include "FluxSvc/SpectrumFactory.h"

static SpectrumFactory<HeSpectrum> factory;
const ISpectrumFactory& HeSpectrumFactory = factory;



// First deal with subclass HeSpectrum::InterpVec::
//  which is a std::vector<float> with a couple of methods added for
//  searching and interpolation.
//-------------------------- InterpVec constructor

HeSpectrum::InterpVec::InterpVec() : std::vector<float>() {}
// default constructor.

//-------------------------- InterpVec::search

HeSpectrum::Intrp HeSpectrum::InterpVec::search(float x) const {
	using namespace std;
    // Binary search a vector for values that straddle x using STL lower_bound.
    // Deal gracefully with values off the end by extrapolation.
    // Contents of vector must be monotonic.
    // Return integer index for use by another vector

	vector<float>::const_iterator loc;  // search. ascending or descending?
    if (back() > front()) loc=lower_bound(begin(), end(), x, less<float>());
    else                  loc=lower_bound(begin(), end(), x, greater<float>());
    // lower_bound points to the entry after x

    int ind;  // convert STL iterator to integer index
    if (loc == begin())    ind = 1;   // extrapolate below table
    else if (loc == end()) ind = (end() - begin()) - 1;  // above table
    else                   ind = (loc - begin());  // in the table

    return Intrp(ind, (x-(*this)[ind-1]) / ((*this)[ind]-(*this)[ind-1]));
}

//-------------------------- InterpVec::interpolate

float HeSpectrum::InterpVec::interpolate(Intrp y) const {
    // linear interpolation between two std::vector values
    return (*this)[y.first-1] +
	y.second * ((*this)[y.first]-(*this)[y.first-1]);
}

//-------------------------- Beginning of HeSpectrum proper

const float HeSpectrum::m_rearth = 6371.f;  // radius of earth in km
const float HeSpectrum::m_altitude = 600.f;  // altitude of circular orbit

//Initializes parameters during construction
void HeSpectrum::init(float lat, float lon) {

    int nen = sizeof(energies)/sizeof(float);
    m_en.reserve(nen);
    float amus = 4.;  // Mass of 4He; CHIME  energy is in MeV/amu
    int i;
    for (i=0; i< nen; i++) m_en.push_back(amus*energies[i]/1000.);
    // convert MeV to GeV

    m_normfact = .5*(1.+sqrt(m_altitude*(m_altitude+2.*m_rearth)) /
		     (m_rearth+m_altitude));
              //geometrical shadow factor in 600 km orbit
    m_fluxes.reserve(nen);
    m_fl.reserve(nen);

    const float width = 0.115f; // CHIME log spacing between standard energies
    m_etop = m_en.back() * (1.+.5*width);  // boundary between table & tail
    m_expo = -2.65f + 1.0f;  // power law exponent (integral) of high-energy tail

    // Populate table of differential galactic fluxes -- no geomagnetic cutoff
    float tfl = 0.f;
    for (i=74; i >= 0; i--) {
	tfl += width*1000.*m_en[i]*m_normfact*fluxes[i];
	m_fluxes.insert(m_fluxes.begin(), tfl);
    }

    // table of total flux as a function of latitude and longitude
    for (int ii=0; ii < 73; ii++) {
	for (int jj=0; jj < 13; jj++) {
	    m_fluxTbl[ii][jj] = gfluxes[jj+13*ii];
	}
    }
    setPosition(lat, lon);

    // cos(angle between zenith and horizon)
    m_coscutoff = -sqrt(m_altitude*m_altitude+2.*m_altitude*m_rearth)
	/ (m_altitude+m_rearth);

    m_particle_name = "alpha";

    // set callback to be notified when the position changes
    m_observer.setAdapter( new ActionAdapter<HeSpectrum>(this,
	&HeSpectrum::askGPS) );

    GPS::instance()->notification().attach( &m_observer );

}


//-------------------------- constructors---------------------

HeSpectrum::HeSpectrum(const std::string& paramstring) {
   std::vector<float> params;

   parseParamList(paramstring,params);

   if(params.empty())
      init(0.0f,0.0f);
   else
      init(params[0], params[1]);
}

std::string HeSpectrum::title() const {
    return "HeSpectrum";
}

const char * HeSpectrum::particleName() const {
    return m_particle_name.c_str();
}
void HeSpectrum::setParticleName(std::string name)
{
    m_particle_name = name;
}

//-------------------------- flux()  (specified cutoff value)

float HeSpectrum::flux(float cut) const {
    // Total flux in nuclei / m^2 sec ster
    // Interpolate in table if possible, otherwise assume power law
    //  tail at high energy.
    if (cut > m_etop) return m_upper * pow(cut/m_etop, m_expo);
    else              return m_upper + m_fl.interpolate(m_en.search(cut));
}

//-------------------------- flux() (current cutoff value)

double HeSpectrum::flux(double) const {
    // calculate flux for the current cutoff
    return m_flux;
}

//--------------------------

double HeSpectrum::solidAngle() const
{
    // the normalization for the flux calculation
    return 4.*M_PI;
}

//-------------------------- flux() (specified position)

double HeSpectrum::flux(double lat, double lon) const {
    // Flux as a function of latitude and longitude in a 600 km orbit.
    // Linear interpolate in a table with a 5 degree sampling grid.
    int ilat = static_cast<int>(lat/5.+6.);
    double/*float*/ a = fmod(lat+30., 5.)/5.;
    int ilon = static_cast<int>(lon/5.);
    double/*float*/ b = fmod(lon, 5.)/5.;
    return m_fluxTbl[ilon][ilat] * (1.-a) * (1.-b) +
	m_fluxTbl[ilon+1][ilat] * (1.-a) * b +
	m_fluxTbl[ilon][ilat+1] * a * (1.-b) +
	m_fluxTbl[ilon+1][ilat+1] * a * b;
}

//-------------------------- operator()  sample an energy value

float HeSpectrum::operator() (float x)const {
    // return a random value of energy sampled from the spectrum

    float join = (m_tot-m_upper)/m_tot;
    if (x < join) return m_en.interpolate(m_fl.search((1.-x)*m_tot-m_upper));
    else          return m_etop*pow((1.-x)/(1.-join), 1./m_expo);
}

//-------------------------- setPosition (separate coordinates)

int HeSpectrum::askGPS()
{
    setPosition(GPS::instance()->lat(), GPS::instance()->lon());
    return 0; // can't be void in observer pattern
}

void HeSpectrum::setPosition(double lat, double lon) {
    // Do the initialization necessary when moving to a new position:
    // look up cutoff energy, build a new table of integral alpha
    // fluxes

    m_lat = lat;
    m_lon = lon;

    // Integrated flux in the power law tail above the table.
    m_upper = -0.115*1000.*m_en.back()*m_normfact*fluxes[74]
	* pow(1.+.5*0.115,m_expo)/(0.115*m_expo);

    m_cutoff = findCutoff(lat,lon);

    // Populate table of integral fluxes modified by geomagnetic cutoff.
    float tfl = 0.;
    m_fl.erase(m_fl.begin(), m_fl.end());
    for (int i = 74; i >= 0; i--) {
	tfl += 0.115*1000.*m_en[i]*m_normfact*fluxes[i]*exposure(m_en[i]);
	m_fl.insert(m_fl.begin(), tfl);
    }

    m_tot = m_fl.front() + m_upper;
    m_flux = flux(m_cutoff);
}

//-------------------------- findCutoff (from a flux)

float HeSpectrum::findCutoff(float rflux) const {
    // determine the cutoff value which will produce the desired flux
    return m_en.interpolate(m_fluxes.search(rflux-m_upper));
}

//-------------------------- findCutoff (from a position)

float HeSpectrum::findCutoff(float lat, float lon) const {
    // determine the cutoff value at a geographical location
    return findCutoff(flux(lat,lon));
}


//------------------------- calculate_rate()

double HeSpectrum::calculate_rate(double old_rate)
{
    return flux(GPS::instance()->lat(), GPS::instance()->lon());
}

//------------------------- rad2()

float HeSpectrum::rad2() const {
    // square of (distance from center of magnetic dipole / earth radius)
    // Dipole is offset from the earth's center

    /*float*/double d2 =
	pow((m_rearth+m_altitude)*sin(m_lat)            - 145.1, 2) +
	pow((m_rearth+m_altitude)*cos(m_lat)*cos(m_lon) + 371.2, 2) +
	pow((m_rearth+m_altitude)*cos(m_lat)*sin(m_lon) - 233.7, 2);
    return d2 / (m_rearth*m_rearth);
}

//------------------------- cosomega()

float HeSpectrum::cosomega(float E) const {
    // Opening angle of geomagnetic cutoff cone.
    // This is a pretty backward way to do it.  It starts from a cutoff
    // energy which is derived from a desired rate.  Then it calculates
    // the magnetic latitude from that.  Really, the cutoff and latitude
    // should be derived from the position by looking in tables.
    // Also, this is simple Størmer theory, ignoring  the penumbral
    // region and multipole effects.

    const float Mp = 4*0.938f;  // mass of alpha in GeV
    const float charge = 2.;    // electric charge of alpha
    float pcut = sqrt(m_cutoff*m_cutoff + 2.*Mp*m_cutoff);  // cutoff momentum
    const float moment = 59.8f; // magnetic moment of earth in convenient units
    float coslambda4 = 4. * pcut  * rad2() / (charge*moment);
      // magnetic latitude**4
    float p = sqrt(E*E + 2.*Mp*E); // momentum
    float coso = 4. * (sqrt(pcut/p) - pcut/p) * pow(coslambda4, -0.75);
      // opening angle of Størmer cone
    if (coso > 1.) return 1.;
    else if (coso < -1.) return -1.;
    else return coso;
}

//-------------------------- exposure()

float HeSpectrum::exposure(float E) const {
    // Geomagnetic cutoff fraction.  Varies from 0 to 1.
    return 0.5 * (1. + cosomega(E));
}

//-------------------------- dir()

std::pair<float,float> HeSpectrum::dir(float energy)const
{

    // Random particle direction from Størmer cone

    // Rejection method for the direction.  Direction must be inside
    // the Stormer cone and above the horizon.
    float coszenith, sinpolar, cospolar, azi;
    const int try_max = 1000;
    static int max_tried = 0;
    int trial=0;
    do {
	// uniform distribution within Størmer cone
        cospolar = -1. + HepRandom::getTheGenerator()->flat()*
            (cosomega(energy)+1.);
	sinpolar = sqrt(1.-cospolar*cospolar);
        azi = 2.*M_PI* HepRandom::getTheGenerator()->flat();
	coszenith = cos(azi)*sinpolar;
    } while (coszenith < m_coscutoff && trial++ < try_max);

    max_tried = std::max(trial, max_tried); // keep track

    // Transform to local earth-based coordinates.
    /*float*/double earthazi = atan2(sinpolar*sin(azi), cospolar);

    return std::make_pair<float,float>(coszenith, earthazi);
}
