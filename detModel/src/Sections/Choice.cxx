#include <string>
#include <map>
#include "detModel/Sections/Choice.h"

namespace detModel{

  void Choice::addCase(std::string name, std::string volume)
  {
    typedef std::map<std::string, std::string> M;
    
    cases.insert(M::value_type(name, volume));
  }

  std::string Choice::getVolumeNameByMode(std::string mode)
  {
    typedef std::map<std::string, std::string> M;
    M::const_iterator i; 
    
    i = cases.find(mode);
    if(i == cases.end()) return "";
    else return i->second;
  }
  
}
