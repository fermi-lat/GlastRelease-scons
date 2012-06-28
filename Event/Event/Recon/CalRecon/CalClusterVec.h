#ifndef CalClusterVec_H
#define CalClusterVec_H


#include <vector>
#include <string>
#include "Event/Recon/CalRecon/CalCluster.h"


/** 
* @class CalClusterVec
*
* @brief TDS class representing a std::vector of pointers to Event::CalCluster obects.
*
* @author Luca Baldini.
*/


namespace Event { //Namespace Event

  class CalClusterVec : virtual public std::vector<CalCluster*>
    {
    public:
      CalClusterVec() {};
      virtual ~CalClusterVec() {};
  };


}; //Namespace Event

#endif        

