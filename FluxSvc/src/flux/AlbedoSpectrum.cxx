// $Header$

#include "FluxSvc/AlbedoSpectrum.h"

#include "facilities/error.h"
#include <cmath>
#include <string.h>
#include <algorithm>
#include "FluxSvc/SpectrumFactory.h"

static SpectrumFactory<AlbedoSpectrum> factory;
const ISpectrumFactory& AlbedoSpectrumFactory = factory;

// local namespace
namespace {

    // define the fromMeV function which will convert from
    // MeV to GeV
    float   fromMeV ( float e ) { return e / 1000.; }
}

AlbedoSpectrum::AlbedoSpectrum(const std::string& name)
{
    if( name == "gamma" )
	type_ = photon;
    else if( name == "n" || name=="neutron" )
	type_ = neutron;
    else
        std::cerr << " AlbedoSpectrum:  not \"gamma\" or \"n\" " << std::endl;
}


const char * AlbedoSpectrum::particleName()const {
    if (type_==photon) return "gamma";
    else return "neutron";
}

std::string AlbedoSpectrum::title()const
{
    return std::string("albedo-")+std::string(particleName());
}


static double albedo_photon_energy(double);
static double albedo_neutron_energy(double);


float AlbedoSpectrum::operator()(float r)const
{
    r = std::max(1e-6, 1.0-r); // make monotonic increasing, avoid problem at zero.
    return fromMeV((type_==photon)? albedo_photon_energy(r)
	: albedo_neutron_energy(r) );
}


// Routines which return random energies typical of the spectra
// of particles produced by earth albedo.

#include <math.h>

// lower cutoffs
static const double neutron_emin= 100.,
photon_emin=10.;
static const double neutronexp =(1.6 - 1.0);

static double albedo_photon_energy(double r)
/* gamma ray spectrum from near the horizon. */
/* Ref: Thompson et al., 1981, JGR, 86, 1265 */
{ return photon_emin / r; }

static double albedo_neutron_energy(double r)
/* neutron spectrum good up to 10 GeV?  */
/* Ref: Armstrong et al., 1973, JGR, 78, 2715
Gehrels, 1992, NIM, A313, 513
Gehrels,                                */
{ return neutron_emin / pow(r, 1./neutronexp); }


double AlbedoSpectrum::calculate_rate(double old_rate)
{
    return old_rate;
}
