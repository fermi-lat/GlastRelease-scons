#include <string>
#include <map>
#include "detModel/Sections/GDDchoice.h"


void GDDchoice::addCase(string name, string volume)
{
  typedef map<string, string> M;

  cases.insert(M::value_type(name, volume));
}

string GDDchoice::getVolumeNameByMode(string mode)
{
  typedef map<string, string> M;
  M::const_iterator i; 

  i = cases.find(mode);
  if(i == cases.end()) return "";
  else return i->second;
}
