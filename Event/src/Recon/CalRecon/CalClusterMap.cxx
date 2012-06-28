

// Implementation file of CalClusterMap class.
//  
// Authors: Luca Baldini.


#include "Event/Recon/CalRecon/CalClusterMap.h"
#include <stdexcept>


const Event::CalClusterVec Event::CalClusterMap::get(std::string key)
{
  CalClusterMap::const_iterator item = find(key);
  if (item == end()){
    std::string msg =  "Invalid CalClusterVec (" + key + ") requested";
    throw std::invalid_argument(msg);
  }
  return item->second;
}

