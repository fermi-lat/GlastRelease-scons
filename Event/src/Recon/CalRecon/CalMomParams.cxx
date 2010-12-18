
// Implementation file of CalMomParams class.
//  
// Authors: Luca Baldini, Johan Bregeon.


#include "Event/Recon/CalRecon/CalMomParams.h"


Event::CalMomParams::CalMomParams(double transRms, double longRms, double longRmsAsym,
				  double longSkewness, double coreEnergyFrac) :
                                  m_transRms(transRms),
                                  m_longRms(longRms),
                                  m_longRmsAsym(longRmsAsym),
                                  m_longSkewness(longSkewness),
                                  m_coreEnergyFrac(coreEnergyFrac)
{
  // Do nothing here, as all the variables are already set :-)
}

void Event::CalMomParams::clear()
{
  // Call the base class clear method...
  Event::CalParams::clear();

  // ... then reset the additional class members.
  m_transRms       = -1.;
  m_longRms        = -1.;
  m_longRmsAsym    = -1.;
  m_longSkewness   = -9999.;
  m_coreEnergyFrac = -1.;
}

double Event::CalMomParams::getElongation() const
{
  if (m_transRms > 0.){
    return m_longRms/m_transRms;
  }
  return -1.;
}

std::ostream& Event::CalMomParams::fillStream(std::ostream& s) const
{
  // Call the base class method...
  Event::CalParams::fillStream(s) << "\n";

  // ... then print the additional stuff.
  s <<
    "Transverse RMS = " << m_transRms << " mm\n" <<
    "Longitudinal RMS = " << m_longRms << " mm\n" <<
    "Longitudinal RMS asymmetry = " << m_longRmsAsym << "\n" <<
    "Longitudinal skewness = " << m_longSkewness << "\n" <<
    "Elongation = " << getElongation() << "\n" <<
    "Core energy fraction = " << m_coreEnergyFrac;  

  return s; 
}
