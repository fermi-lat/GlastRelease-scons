
// Implementation file of CalClassParams class.
//  
// Authors: Luca Baldini, Johan Bregeon.


#include "Event/Recon/CalRecon/CalNBCClassParams.h"


Event::CalNBCClassParams::CalNBCClassParams()
{
  setProducerName("NBC");
  clear();
}

Event::CalNBCClassParams::CalNBCClassParams(std::string producerName,
                                            std::map <std::string, double> probMap) :
  Event::CalClassParams(producerName, probMap)
{
  //Nothing to do.
}

void Event::CalNBCClassParams::multiply(const std::string &className, double pdfValue)
{
  // If the base map has not an entry for this class, yet, initialize the
  // probability to the pdf value. 
  if ( !hasClass(className) ) {
    setProb(className, pdfValue);
  }
  // Otherwise multiply the current probability value by the pdf value. 
  else {
    setProb(className, pdfValue*getProb(className));
  }
}

void Event::CalNBCClassParams::normalize()
{
  // First round to evaluate the sum of probabilities.
  double norm = 0.0;
  std::map <std::string, double>::const_iterator iter;
  std::map<std::string, double> probMap = getProbMap();
  for (iter = probMap.begin(); iter != probMap.end(); iter++)
    {
      norm += (*iter).second;
      //std::cout << (*iter).first << " " << (*iter).second << std::endl;
    }
  // Then the actual normalization, if it's the case.
  if ( norm > 0 ) {
    for (iter = probMap.begin(); iter != probMap.end(); iter++)
      {
        setProb( (*iter).first, ((*iter).second) / norm );
      }
  }
}
