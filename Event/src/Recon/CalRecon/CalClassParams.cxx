
// Implementation file of CalClassParams class.
//  
// Authors: Luca Baldini, Johan Bregeon.


#include "Event/Recon/CalRecon/CalClassParams.h"


Event::CalClassParams::CalClassParams()
{
  setProducerName("Not set");
  clear();
}

Event::CalClassParams::CalClassParams(std::string producerName,
				      std::map <std::string, double> probMap) :
  m_producerName(producerName),
  m_probMap(probMap)
{
  //Nothing to do.
}

void Event::CalClassParams::clear()
{
  m_probMap.clear();
}

bool Event::CalClassParams::hasClass(const std::string &className) const
{
  return ( m_probMap.count(className) > 0 );
}

double Event::CalClassParams::getProb(const std::string &className) const
{
  if ( hasClass(className) ) {
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
