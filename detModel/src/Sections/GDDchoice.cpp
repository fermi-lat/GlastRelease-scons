#include <string>
#include <map>
#include "detModel/Sections/GDDchoice.h"


void GDDchoice::addCase(std::string name, std::string volume)
{
  typedef std::map<std::string, std::string> M;

  cases.insert(M::value_type(name, volume));
}

string GDDchoice::getVolumeNameByMode(string mode)
{
  typedef std::map<std::string, std::string> M;
  M::const_iterator i; 

  i = cases.find(mode);
  if(i == cases.end()) return "";
  else return i->second;
}
