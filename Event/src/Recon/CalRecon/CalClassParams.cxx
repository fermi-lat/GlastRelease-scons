
// Implementation file of CalClassParams class.
//  
// Authors: Luca Baldini, Johan Bregeon.


#include "Event/Recon/CalRecon/CalClassParams.h"


Event::CalClassParams::CalClassParams()
{
  clear();
  setProducerName("Not set");
}

void Event::CalClassParams::clear()
{
  m_probMap.clear();
  setProb("gam", -1.);
}

bool Event::CalClassParams::hasClass(const std::string &className) const
{
  return m_probMap.count(className);
}

double Event::CalClassParams::getProb(const std::string &className) const
{
  if (hasClass(className)){
    return m_probMap.find(className)->second;
  }
  else {
    return -1.;
  }
}

std::ostream& Event::CalClassParams::fillStream(std::ostream& s) const
{
  s << "Producer name: '" << getProducerName() << "'\n";
  std::map <std::string, double>::const_iterator iter;
  for (iter = m_probMap.begin(); iter != m_probMap.end(); iter++)
    {
      s << "Probability for class '" << (*iter).first << "': " << (*iter).second << "\n";
    }

  return s; 
}
