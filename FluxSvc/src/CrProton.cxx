//
// CrProton.cc
//
//   The entry-point class for primary particle generators of primary and
//   secondary (i.e. albedo) cosmic-ray protons.
//

#include <math.h>
#include <map>
#include <vector>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>

#include "CrProton.h"
#include "CrProtonPrimary.h"
#include "CrProtonSplash.h"
#include "CrProtonReentrant.h"

// derived from CrSpectrum that is derived from Spectrum
#include "CrSpectrum.h"

// define a factory for anonomous instantiation
#include "FluxSvc/SpectrumFactory.h"

static SpectrumFactory<CrProton> factory;
const ISpectrumFactory& CrProtonFactory = factory;



CrProton::CrProton(const std::string& paramstring)
: m_component(0)
{
   std::vector<float> params;

   parseParamList(paramstring,params);

    int flag = params.empty() || params[0]==0 ? 7 : params[0];

    // including each component (primary/re-entrant/splash protons)...
   if(flag& 1) m_subComponents.push_back(new CrProtonPrimary);
   if(flag& 2) m_subComponents.push_back(new CrProtonReentrant);
   if(flag& 4) m_subComponents.push_back(new CrProtonSplash);
}


CrProton::~CrProton()
{
  std::vector<CrSpectrum*>::iterator  i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    delete *i;
  }
}


CrSpectrum* CrProton::selectComponent(HepRandomEngine* engine)
{
  std::map<CrSpectrum*,double>       integ_flux;
  double                             total_flux = 0;
  std::vector<CrSpectrum*>::iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
    integ_flux[*i] = total_flux;
  }

  double  rnum = engine->flat() * total_flux;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    if (integ_flux[*i] >= rnum) { break; }
  }

  m_component = *i;

  return *i;
}


double CrProton::energySrc(HepRandomEngine* engine)
{
//THB  if (!m_component){ selectComponent(engine); }
 selectComponent(engine);
  return m_component->energySrc(engine);
}


std::pair<double,double> CrProton::dir(double energy, HepRandomEngine* engine)
// return: cos(zenith_angle) and azimuth [rad]
{
  if (!m_component){ selectComponent(engine); }

  return m_component->dir(energy, engine);
}


CrSpectrum* CrProton::component() const
{
  return m_component;
}

double CrProton::flux ( ) const
{
  double          total_flux = 0;
  std::vector<CrSpectrum*>::const_iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
  }
  return total_flux;
}

double CrProton::solidAngle( )const
{
    return 4*M_PI;
}

