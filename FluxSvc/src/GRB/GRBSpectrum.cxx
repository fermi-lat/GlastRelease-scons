///
///   GRBSpectrum: Spectrum class for the GRB source
///    Authors: Nicola Omodei & Johann Cohen Tanugi 
///
#include "GRBSpectrum.h"
#include "GRBConstants.h"
#include <math.h>

// define a factory for anonomous instantiation
#include "FluxSvc/SpectrumFactory.h"

static SpectrumFactory<GRBSpectrum> factory;
const ISpectrumFactory& GRBSpectrumFactory = factory;

static const GRBConstants cst;

//static const   TROOT root("GUI", "GUI test environement");


GRBSpectrum::GRBSpectrum(const std::string& params) 
{
  m_grbsim = new GRBSim();
  m_grbsim->Start();
}

GRBSpectrum::~GRBSpectrum() 
{
  delete m_grbsim;
}
double GRBSpectrum::solidAngle() const
{
  return 1.0;
}

///return flux, given a time
double GRBSpectrum::flux(double time) const
{
  m_grbsim->ComputeFlux(time);
  cout<<"**** FLUX ***** "<<m_grbsim->IRate()<<endl;
  return m_grbsim->IRate(); // in ph/(m^2 s) 
}

///testing rate return flux, given a time
double GRBSpectrum::rate(double time) const
{
  m_grbsim->ComputeFlux(time);
  cout<<"**** RATE ***** "<<m_grbsim->IRate()<<endl;
  return m_grbsim->IRate(); // in ph/(m^2 s) 
}

double GRBSpectrum::energySrc(HepRandomEngine* engine, double time)
{
  //return the flux in photons/(m^2 sec)
  /// Use time to update the spectrum
  m_spectrum.clear();
  m_grbsim->ComputeFlux(time);
  m_spectrum =m_grbsim->Spectrum();
  /// Then returns a uniform random value, to be used by operator()
  return (*this)(engine->flat());
}

float GRBSpectrum::operator() (float u) const
{
  float energy = m_grbsim->DrawPhotonFromSpectrum(m_spectrum);
  return (float)energy;
}

std::pair<float,float> GRBSpectrum::dir(float energy)const
{
  return std::make_pair(0.5,0.5);
}
