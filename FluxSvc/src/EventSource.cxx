//  $Header$

#include "../FluxSvc/EventSource.h"

#include "dom/DOM_Element.hpp"
#include "xml/Dom.h"
#include "GPS.h"
#include "Orbit.h"
#include "CLHEP/Random/RandExponential.h"

#include <strstream>

unsigned int  EventSource::s_id = 0;
double  EventSource::s_total_area = 6.; // area in m^2

EventSource::EventSource (double aFlux, unsigned acode)
  :  m_enabled(true), m_flux(aFlux), m_solid_angle(1.), m_code(acode)
{
  std::strstream  s;

  s << "Source_" << (++s_id) << '\0';
  if (acode == 0) code(s_id); // automatically assign event codes...

  m_name = s.str();
  s.freeze(false);
}

EventSource::EventSource (const DOM_Element& xelem)
  :  m_enabled(true), m_flux(1.), m_solid_angle(1.), m_code(0)
{
  m_name = xml::Dom::getAttribute(xelem, "name");
  m_flux = atof (xml::Dom::getAttribute(xelem, "flux").c_str());

  std::string code_str = xml::Dom::getAttribute(xelem, "code");
  if (code_str != std::string("")) {
    m_code = atoi(code_str.c_str());
  }
  else  {
    m_code = ++s_id;
  }

  // this is set by default to be overriden when the solid_angle 
  // element is present...
  DOM_Element   angles = 
    xml::Dom::findFirstChildByName(xelem, "solid_angle");

  if (angles != DOM_Element()) {
    double  mincos = atof(xml::Dom::getAttribute(angles, "mincos").c_str());
    double  maxcos = atof(xml::Dom::getAttribute(angles, "maxcos").c_str());

    m_solid_angle = (maxcos - mincos)*2*M_PI;
  }
  else if (xml::Dom::findFirstChildByName(xelem, "direction") != DOM_Element())
    m_solid_angle = 1.;
}


EventSource::~EventSource()
{}

double EventSource::flux () const
{
  return m_flux;  // default if not overridden
}

void   EventSource::flux (double value) {
  m_flux = value;
}

GPStime EventSource::interval ()
{
  double  r = rate();

  if (r > 0) {
    return  RandExponential::shoot(1./r);
  } else if (r < 0) {
    return 1./r;
  } else return -1.;
}

double  EventSource::rate ( double solid_angle, double fluxval )
{
  double  r = solid_angle*fluxval*s_total_area;   
  // integral of the flux*cross-sectional-area over
  // the solid angle (which is passed in as an integral)
  m_solid_angle = solid_angle;
  flux( fluxval );
  return  r;
}

double  EventSource::rate ( )const
{
  return enabled()? (solidAngle()*flux()*s_total_area) :0;
}

void    EventSource::rate ( double rate )
{
  flux(  rate/(m_solid_angle*s_total_area) );
}

Orbit*  EventSource::makeOrbit () const
{
  return new Orbit;
}



// UI titles
std::string EventSource::fullTitle () const 
{ return std::string("EventSource");   }
std::string EventSource::displayTitle () const  {  return m_name; }
