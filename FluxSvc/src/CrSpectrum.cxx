//$Header$

// FIX ME!!:
//   - computte the geomagnetic coordinates from the input (lon,lat)
#include "CrSpectrum.h"

CrSpectrum::CrSpectrum()
{
  setPosition(0., 0., 0.); // (latitude, longitude, time)
}


CrSpectrum::~CrSpectrum()
{
  ;
}


void CrSpectrum::setTime(double time)
{
  setPosition(m_lat, m_lon, time);
}


void CrSpectrum::setPosition(double latitude, double longitude)
{
  setPosition(latitude, longitude, m_time);
}


void CrSpectrum::setPosition(double latitude, double longitude, double time)
{
  m_lon  = longitude;
  m_lat  = latitude;
  m_time = time;

  // compute the geomagnetic coordinates
  // temporaly returns geographic latitude and longitude: FIX ME!!
  m_geomagneticLongitude = longitude;
  m_geomagneticLatitude  = latitude;

  // compute the Cut-Off-Rigidity (GV)
  // temporaly retruns fixed value: FIX ME!!
  // ------------------------------
  // For proton incoming vertically to the magnetic field,
  // cut off rigidity (Rc) is calculated as
  // (Rc/GV) = 14.9 * (1+h/R)^-2 * (cos(theta_M))^4,
  // where h gives altitude from earth surface, R means an radius of earth,
  // and theta_M is geomagnetic lattitude.
  // For a balloon experiment at Palestine, h is about 50km, 
  // theta_M is about 0.735, and R is about 6380km.
  // We thus obtaine Rc=4.44 GV.
  // Reference: 
  // "Handbook of space astronomy and astrophysics" 2nd edition, p225
  //   Zombeck, 1990, Cambridge University Press 
  // "High Energy Astrophysics" 2nd edition, p325-330
  //   M. S. Longair, 1992, Cambridge University Press
  // ------------------------------
  m_cutOffRigidity = 4.44;  // temporarily fixed to Palestine value.

  // compute the force-field approximation potential (MV)
  // temporaly returns fized value: FIX ME!!
  //  m_solarWindPotential = 540.;  // temporarily fixed to the typical minimum-activity value.
  m_solarWindPotential = 1100.;  // temporarily fixed to the typical minimum-activity value.
}


double CrSpectrum::longitude() const
{
  return m_lon;
}


double CrSpectrum::latitude() const
{
  return m_lat;
}


double CrSpectrum::geomagneticLongitude() const
{
  return m_geomagneticLongitude;
}


double CrSpectrum::geomagneticLatitude() const
{
  return m_geomagneticLatitude;
}


double CrSpectrum::cutOffRigidity() const
{
  return m_cutOffRigidity;
}


double CrSpectrum::solarWindPotential() const
{
  return m_solarWindPotential;
}

