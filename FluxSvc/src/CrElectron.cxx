/**
 * CrElectron:
 * The entry-point class for primary particle generators of primary and
 * secondary (i.e. albedo) cosmic-ray electrons.
 *
 * The included components (primary, re-entrand and splash protons)
 * have to be defined in the constructor.
 *
 * Ver 1.0 on 2001-04-18 by Masanobu Ozaki <ozaki@astro.isas.ac.jp>
 *
 * $Header$
 */

#include <math.h>
#include <map>
#include <vector>

// CLHEP
#include <CLHEP/Random/RandomEngine.h>
#include <CLHEP/Random/RandGeneral.h>

#include "CrElectron.h"
#include "CrElectronPrimary.h"
#include "CrElectronSplash.h"
#include "CrElectronReentrant.h"
  // derived from CrSpectrum that is derived from Spectrum
#include "CrSpectrum.h"

// define a factory for anonomous instantiation
#include "SpectrumFactory.h"

static SpectrumFactory<CrElectron> factory;
const ISpectrumFactory& CrElectronFactory = factory;




CrElectron::CrElectron(const std::string& paramstring)
: m_component(0)
{
   std::vector<float> params;

   parseParamList(paramstring,params);

   int flag = params.empty() || params[0]==0 ? 7 : params[0];

    // including each component (primary/re-entrant/splash protons)...
   if(flag& 1) m_subComponents.push_back(new CrElectronPrimary);
   if(flag& 2) m_subComponents.push_back(new CrElectronReentrant);
   if(flag& 4) m_subComponents.push_back(new CrElectronSplash);

}


CrElectron::~CrElectron()
{
  std::vector<CrSpectrum*>::iterator  i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    delete *i;
  }
}


CrSpectrum* CrElectron::selectComponent(HepRandomEngine* engine)
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


double CrElectron::energySrc(HepRandomEngine* engine)
{
//THB: always select  if (!m_component){ selectComponent(engine); }
	selectComponent(engine); 

  return m_component->energySrc(engine);
}


std::pair<double,double> CrElectron::dir(double energy, HepRandomEngine* engine)
// return: cos(zenith_angle) and azimuth [rad]
{
  if (!m_component){ selectComponent(engine); }

  return m_component->dir(energy, engine);
}


CrSpectrum* CrElectron::component() const
{
  return m_component;
}

double CrElectron::flux ( ) const
{
  double          total_flux = 0;
  std::vector<CrSpectrum*>::const_iterator i;
  for (i = m_subComponents.begin(); i != m_subComponents.end(); i++){
    total_flux += (*i)->flux();
  }
  return total_flux;
}

double CrElectron::solidAngle( )const
{
    return 4	*M_PI;
}

